"""Kommandozeilen-Einstieg der ROS Research Suite.

Beispiele:
    python main.py --dataset GSE12220
    python main.py --dataset GSE12220 --analysis log2fc
    python main.py --dataset GSE12220 --inspect-annotation
"""
import argparse
import logging
import sys
import time
from pathlib import Path
from typing import NamedTuple

import pandas as pd

import config
from src import (annotation, differential_expression, metadata, provenance, robustness,
                 ros_analysis, ros_evidence, ros_genes, statistics, timecourse)
from src import design as design_module
from src.geo_client import GeoClient, GeoError, normalize_geo_id
from src.utils import setup_logging

logger = logging.getLogger("main")

AVAILABLE_ANALYSES = ("log2fc", "differential", "timecourse", "ros", "robustness", "validation")
PLANNED_ANALYSES = ("full",)


class ImportedDataset(NamedTuple):
    """Alles, was nach dem Import einer Series-Matrix-Datei vorliegt."""
    label: str
    samples: pd.DataFrame      # rohe Sample-Metadaten
    design: pd.DataFrame       # bestätigte Zuordnung (Behandlung, Zeit, ...)
    expression: pd.DataFrame   # Sonden x Samples, Spalten in Design-Reihenfolge
    series_info: dict          # Serien-Metadaten aus GEO (Titel, PubMed-ID, ...)
    series_file: Path          # lokale Series-Matrix-Datei


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ROS Research Suite – Analyse öffentlicher Genexpressionsdaten."
    )
    parser.add_argument("--dataset", required=True, help="GEO-ID, z. B. GSE12220")
    parser.add_argument("--analysis", choices=AVAILABLE_ANALYSES + PLANNED_ANALYSES,
                        help="Analyse; verfügbar: log2fc, differential, timecourse, ros, robustness, validation")
    parser.add_argument("--inspect-annotation", action="store_true",
                        help="Spalten der Plattform-Annotation anzeigen")
    parser.add_argument("--report", action="store_true", help="(ab Phase 9 verfügbar)")
    parser.add_argument("--force-download", action="store_true",
                        help="Cache ignorieren und neu herunterladen")
    parser.add_argument("--reconfirm", action="store_true",
                        help="Zuordnung (Kontrolle/Behandlung/Zeit) erneut bestätigen")
    parser.add_argument("--yes", action="store_true",
                        help="Zuordnung ohne Nachfrage bestätigen (nur wenn sie eindeutig ist)")
    parser.add_argument("--verbose", action="store_true", help="Ausführlichere Logs")
    return parser


def dataset_label(path) -> str:
    """'GSE12220_series_matrix.txt.gz' -> 'GSE12220' (bei mehreren Plattformen 'GSE1-GPL2')."""
    return path.name.replace("_series_matrix.txt.gz", "")


def run_import(client: GeoClient, geo_id: str, args: argparse.Namespace) -> list[ImportedDataset]:
    """Phase 1+2: Daten holen, Zuordnung bestätigen, Expressionsdaten speichern."""
    files = client.fetch_series_matrix(geo_id, force=args.force_download)
    datasets = []

    for path in files:
        label = dataset_label(path)
        series_info, sample_rows = metadata.read_header(path)
        samples = metadata.build_sample_table(sample_rows)

        print(f"\n=== {label} ===")
        print(f"Titel:   {series_info.get('Series_title', '(unbekannt)')}")
        print(f"Samples: {len(samples)}")

        meta_dir = config.METADATA_DIR / geo_id
        meta_dir.mkdir(parents=True, exist_ok=True)
        samples.to_csv(meta_dir / f"{label}_samples.csv")

        draft = design_module.build_design(samples)
        design = design_module.confirm_design(
            draft, meta_dir, label,
            interactive=sys.stdin.isatty(), assume_yes=args.yes, reconfirm=args.reconfirm)

        expression = metadata.read_expression_table(path)
        expression = metadata.align_expression(expression, design)
        info = metadata.describe_expression(expression)

        out_dir = config.PROCESSED_DIR / geo_id
        out_dir.mkdir(parents=True, exist_ok=True)
        expression.to_csv(out_dir / f"{label}_expression.csv")

        print(f"\nExpression: {info['n_probes']} Sonden x {info['n_samples']} Samples, "
              f"fehlende Werte: {info['n_missing']}")
        print(f"Wertebereich: {info['min']:.2f} bis {info['max']:.2f} "
              f"(Median {info['median']:.2f}) -> "
              + ("sieht nach log2-Skala aus (Faustregel)" if info["likely_log_scale"]
                 else "KEINE log2-Skala erkennbar, Normalisierung nötig"))
        datasets.append(ImportedDataset(label, samples, design, expression, series_info, path))
    return datasets


def platform_id_of(samples: pd.DataFrame) -> str:
    """Die (eine) Plattform-ID der Samples, z. B. 'GPL2529'."""
    if "platform_id" not in samples.columns:
        raise annotation.AnnotationError("In den Metadaten steht keine Plattform-ID.")
    platforms = samples["platform_id"].unique()
    if len(platforms) != 1:
        raise annotation.AnnotationError(f"Mehrere Plattformen in einer Datei: {list(platforms)}")
    return str(platforms[0])


def load_annotation(client: GeoClient, samples: pd.DataFrame, force: bool) -> pd.DataFrame:
    path = client.fetch_platform_annotation(platform_id_of(samples), force=force)
    return annotation.parse_platform_table(path)


def inspect_annotation(client: GeoClient, dataset: ImportedDataset, force: bool) -> None:
    """Zeigt, welche Spalten die Annotation hat und welche davon erkannt wurden."""
    table = load_annotation(client, dataset.samples, force)
    rate = annotation.probe_match_rate(dataset.expression.index, table)
    print(f"\n=== Annotation {platform_id_of(dataset.samples)} ===")
    print(f"Einträge: {len(table)} | Spalten: {list(table.columns)}")
    print(f"Anteil unserer Sonden in der Annotation: {rate:.1%}")
    print(f"Erkannt: {annotation.detect_gene_columns(table)}")
    print("\nDie ersten 3 Einträge:")
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(table.head(3).T.to_string(max_colwidth=60))


def build_log2fc_table(client: GeoClient, dataset: ImportedDataset, force: bool) -> pd.DataFrame:
    """log2FC jedes Zeitpunkts gegen t=0 (je Behandlung), wenn möglich mit Gennamen."""
    table = statistics.log2fc_table(dataset.expression, dataset.design)

    try:  # Ohne Annotation geht es mit Sonden-IDs weiter (ein fehlender Dienst stoppt nichts)
        annotation_table = load_annotation(client, dataset.samples, force)
        rate = annotation.probe_match_rate(table.index, annotation_table)
        if rate < config.MIN_ANNOTATION_MATCH_RATE:
            raise annotation.AnnotationError(
                f"Nur {rate:.1%} der Sonden passen zur Annotation (Minimum "
                f"{config.MIN_ANNOTATION_MATCH_RATE:.0%}). Falsche Plattform?")
        columns = annotation.detect_gene_columns(annotation_table)
        table = annotation.attach_annotation(table, annotation_table, columns)
        print(f"\nGen-Spalten aus der Annotation: {columns}")
        if "organism" in table.columns:
            print("\nOrganismen der Sonden laut Annotation:")
            print(table["organism"].value_counts(dropna=False).to_string())
            if config.ORGANISM_FILTER:
                table, removed = annotation.filter_by_organism(table, config.ORGANISM_FILTER)
                print(f"Filter '{config.ORGANISM_FILTER}': {removed} Sonden entfernt, "
                      f"{len(table)} bleiben.")
    except (annotation.AnnotationError, GeoError) as error:
        logger.warning("Weiter OHNE Gennamen: %s", error)

    table.index.name = "probe_id"
    return table


def run_log2fc(client: GeoClient, geo_id: str, dataset: ImportedDataset, force: bool) -> None:
    """Berechnet log2FC jedes Zeitpunkts gegen t=0 (je Behandlung) und speichert die Tabelle."""
    table = build_log2fc_table(client, dataset, force)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out_path = config.TABLES_DIR / f"{dataset.label}_log2fc_vs_t0.csv"
    table.to_csv(out_path)

    fc_columns = [c for c in table.columns if c.startswith("log2fc_")]
    values = table[fc_columns].to_numpy(dtype=float)
    print(f"\nlog2FC-Tabelle: {table.shape[0]} Sonden, {len(fc_columns)} Zeitpunkt-Spalten")
    print(f"log2FC-Bereich: {values.min():.2f} bis {values.max():.2f}")
    print(f"Gespeichert: {out_path}")


def print_top(result: pd.DataFrame, treatment: str, direction: str) -> pd.DataFrame:
    """Zeigt (Top 10) und gibt die Top-N-Tabelle zurück."""
    top = differential_expression.top_genes(result, treatment, direction, config.TOP_N_GENES)
    value = f"{treatment}_max_log2fc" if direction == "up" else f"{treatment}_min_log2fc"
    show = [c for c in ("display_name", "systematic_name") if c in top.columns]
    show += [value, f"{treatment}_peak_time_min", f"{treatment}_n_timepoints_above"]
    word = "hochreguliert" if direction == "up" else "herunterreguliert"
    print(f"\nTop 10 {word} ({treatment}):")
    print(top[show].head(10).round(2).to_string())
    return top


class Prepared(NamedTuple):
    """Gemeinsame Vorbereitung für Phase 4-6b."""
    table: pd.DataFrame        # log2FC je Zeitpunkt, mit Gennamen, ohne Kontrollsonden
    treatments: list[str]
    result: pd.DataFrame       # responsive Gene, Richtung, Spezifität, statistische Evidenz
    counts: dict               # Anzahl Sonden nach jedem Filterschritt


def prepare_responses(client: GeoClient, dataset: ImportedDataset, force: bool) -> Prepared:
    """Tabelle bauen, Kontrollsonden abtrennen, responsive Gene bestimmen."""
    table = build_log2fc_table(client, dataset, force)
    n_after_organism_filter = len(table)
    table, controls = differential_expression.split_control_probes(table)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    controls.to_csv(config.TABLES_DIR / f"{dataset.label}_control_probes_log2fc.csv")
    print(f"\n{len(controls)} Kontrollsonden ({config.CONTROL_PROBE_PREFIX}…) abgetrennt, "
          f"{len(table)} Sonden werden analysiert.")
    treatments = [t for t in dataset.design["treatment"].unique()
                  if differential_expression.timepoint_columns(table, t)]
    result = differential_expression.classify_all(table, treatments)
    for treatment in treatments:   # Standard: statistische Evidenz immer ausweisen
        result[f"{treatment}_statistical_evidence"] = statistics.replication_status(
            dataset.design, treatment)
    counts = {"probes_in_expression_file": len(dataset.expression),
              "probes_after_organism_filter": n_after_organism_filter,
              "control_probes_removed": len(controls), "probes_analyzed": len(table)}
    return Prepared(table, treatments, result, counts)


def files_changed_since(start: float) -> list[str]:
    """Ergebnisdateien (ohne JSON), die seit `start` geschrieben wurden."""
    return sorted(str(p.relative_to(config.PROJECT_ROOT)) for p in config.RESULTS_DIR.rglob("*")
                  if p.is_file() and p.suffix != ".json" and p.stat().st_mtime >= start)


def finish_run(client: GeoClient, geo_id: str, dataset: ImportedDataset, analysis: str,
               counts: dict | None, start: float) -> None:
    """Schreibt/ergänzt die Metadaten-Datei dieses Datensatzes (Reproduzierbarkeit)."""
    annotation_file, annotation_table = None, None
    try:
        annotation_file = client.fetch_platform_annotation(platform_id_of(dataset.samples))
        annotation_table = annotation.parse_platform_table(annotation_file)
    except (annotation.AnnotationError, GeoError):
        logger.warning("Annotation für die Metadaten nicht verfügbar.")
    info = provenance.build_dataset_info(geo_id, dataset.series_info, dataset.samples, dataset.design,
                                         dataset.series_file, annotation_file, annotation_table)
    path = provenance.update_metadata(
        config.RESULTS_DIR / f"{dataset.label}_analysis_metadata.json", info, analysis, counts,
        files_changed_since(start))
    print(f"Metadaten gespeichert: {path}")


def run_differential(client: GeoClient, dataset: ImportedDataset, force: bool) -> dict:
    """Phase 4: responsive Gene, Top-Listen, H2O2-vs-MMS-Vergleich, MA-/Antwort-Plots."""
    prepared = prepare_responses(client, dataset, force)
    table, treatments, result = prepared.table, prepared.treatments, prepared.result
    label = dataset.label
    result.to_csv(config.TABLES_DIR / f"{label}_differential_expression.csv")

    print(f"\nKriterium: |log2FC| >= {config.MIN_ABS_LOG2FC} in mindestens "
          f"{config.MIN_CONSECUTIVE_TIMEPOINTS} aufeinanderfolgenden Zeitpunkten "
          "(keine p-Werte: keine Replikate).")
    for treatment in treatments:
        counts = result[f"{treatment}_direction"].value_counts()
        print(f"\n{treatment}: " + ", ".join(f"{name}: {counts.get(name, 0)}"
                                              for name in ("up", "down", "both", "none")))
        for direction in ("up", "down"):
            top = print_top(result, treatment, direction)
            top.to_csv(config.TABLES_DIR / f"{label}_top_{direction}_{treatment}.csv")

    print("\nStatistische Evidenz je Gen: "
          + ", ".join(f"{t}: {statistics.replication_status(dataset.design, t)}" for t in treatments)
          + " -> 'responsiv' beschreibt Effektstärke + Konsistenz, nicht Signifikanz.")
    specificity = [c for c in result.columns if c.startswith("specificity_")]
    if specificity:
        low = config.SPECIFICITY_LOW_FRACTION * config.MIN_ABS_LOG2FC
        print(f"\n{specificity[0]} (spezifisch = andere Behandlung an allen Zeitpunkten |log2FC| < {low:g}; "
              "dazwischen 'ambiguous'):")
        print(result[specificity[0]].value_counts().to_string())

    try:  # Grafiken sind optional: ohne matplotlib läuft der Rest weiter
        from src import visualization
        config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        for treatment in treatments:
            visualization.ma_plots(
                dataset.expression, dataset.design, result, treatment, config.MIN_ABS_LOG2FC,
                config.FIGURES_DIR / f"{label}_MA_{treatment}.png")
            visualization.response_scatter(
                result, treatment, config.MIN_ABS_LOG2FC,
                config.FIGURES_DIR / f"{label}_response_{treatment}.png")
        print(f"\nGrafiken gespeichert in {config.FIGURES_DIR}")
    except ImportError:
        logger.warning("matplotlib fehlt – Grafiken übersprungen (pip install -r requirements.txt).")
    print(f"Tabellen gespeichert in {config.TABLES_DIR}")
    return prepared.counts


def run_timecourse(client: GeoClient, dataset: ImportedDataset, force: bool) -> dict:
    """Phase 5: Zeitverlaufs-Muster, Sensitivitätsanalyse, Verlaufsplots."""
    prepared = prepare_responses(client, dataset, force)
    table, treatments, result = prepared.table, prepared.treatments, prepared.result
    label = dataset.label
    patterns = timecourse.classify_patterns_all(table, treatments)
    patterns.to_csv(config.TABLES_DIR / f"{label}_timecourse_patterns.csv")

    print(f"\nKriterien: Schwelle |log2FC| >= {config.MIN_ABS_LOG2FC}, "
          f"{config.MIN_CONSECUTIVE_TIMEPOINTS} Zeitpunkte in Folge; 'early' = Beginn <= "
          f"{config.EARLY_ONSET_MAX_MIN} min; 'transient' = am Ende <= "
          f"{config.TRANSIENT_RETURN_FRACTION:.0%} des Maximums und unter der Schwelle.")
    for treatment in treatments:
        column = f"{treatment}_pattern"
        counts = patterns[column].value_counts()
        print(f"\n{treatment} – Muster (Anzahl Gene):")
        print(counts.to_string())
        print("Beispiele (stärkste Antwort zuerst):")
        for pattern in counts.index:
            if pattern == "none":
                continue
            members = patterns.index[patterns[column] == pattern]
            strongest = result.loc[members, f"{treatment}_peak_log2fc"].abs().nlargest(5).index
            names = ", ".join(str(patterns.loc[probe, "display_name"])
                              if "display_name" in patterns.columns else probe
                              for probe in strongest)
            print(f"  {pattern}: {names}")

    sensitivity = differential_expression.threshold_sensitivity(table, treatments)
    sensitivity.to_csv(config.TABLES_DIR / f"{label}_threshold_sensitivity.csv", index=False)
    print("\nSensitivität: Anzahl responsiver Gene bei anderen Kriterien")
    print(sensitivity.pivot_table(index=["threshold", "min_consecutive"], columns="treatment",
                                  values="n_responsive").to_string())

    try:
        from src import visualization
        config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        lead = "H2O2" if "H2O2" in treatments else treatments[0]
        up = differential_expression.top_genes(result, lead, "up", 8).index
        down = differential_expression.top_genes(result, lead, "down", 4).index
        visualization.gene_timecourse_grid(
            table, list(up) + list(down), treatments, config.MIN_ABS_LOG2FC,
            config.FIGURES_DIR / f"{label}_timecourse_top_genes_{lead}.png")
        for treatment in treatments:
            visualization.pattern_profiles(
                table, patterns[f"{treatment}_pattern"], treatment, config.MIN_ABS_LOG2FC,
                config.FIGURES_DIR / f"{label}_pattern_profiles_{treatment}.png")
        print(f"\nGrafiken gespeichert in {config.FIGURES_DIR}")
    except ImportError:
        logger.warning("matplotlib fehlt – Grafiken übersprungen.")
    print(f"Tabellen gespeichert in {config.TABLES_DIR}")
    return prepared.counts


def run_ros(client: GeoClient, dataset: ImportedDataset, force: bool) -> dict:
    """Phase 6: bekannte ROS-Gene auswerten und datengetrieben vergleichen."""
    prepared = prepare_responses(client, dataset, force)
    table, treatments, result = prepared.table, prepared.treatments, prepared.result
    label = dataset.label
    patterns = timecourse.classify_patterns_all(table, treatments)
    ros_table, missing = ros_genes.build_ros_table(result, patterns)

    # --- 1. Zuordnung der Liste zur Messung (Qualitätskontrolle) ---
    n_listed = len(ros_genes.ROS_GENES)
    best = ros_table[ros_table["is_best_probe"]]
    print(f"\nROS-Liste: {n_listed} Gene, {len(best)} davon auf dem Chip gefunden.")
    if missing:
        print("Nicht auf dem Chip gefunden: " + ", ".join(f"{g.symbol} ({g.orf})" for g in missing))
    mismatches = ros_table[ros_table["symbol_check"] == "symbol_mismatch"]
    for probe, row in mismatches.iterrows():
        print(f"ACHTUNG Symbol-Abweichung: Liste sagt {row['symbol']} ({row['orf']}), "
              f"Annotation sagt '{row['annotation_symbol']}' (Sonde {probe}) – bitte in der SGD prüfen.")
    multi = int((best["n_probes"] > 1).sum())
    if multi:
        print(f"{multi} Gene haben mehrere Sonden; gezählt wird jeweils die mit dem stärksten Signal.")

    ros_table.to_csv(config.TABLES_DIR / f"{label}_ros_genes_response.csv")
    summary = ros_analysis.category_summary(ros_table, treatments)
    summary.to_csv(config.TABLES_DIR / f"{label}_ros_category_summary.csv")
    print("\nBekannte ROS-Gene je Kategorie (Anzahl responsiver Gene, median |peak log2FC|):")
    print(summary.to_string())

    for treatment in treatments:
        in_list, rest = ros_analysis.responsive_share(result, ros_table, treatment)
        print(f"\n{treatment}: responsiv sind {in_list:.0%} der ROS-Listen-Gene, "
              f"{rest:.0%} aller übrigen Gene.")
        known = best.reindex(best[f"{treatment}_peak_log2fc"].abs().sort_values(ascending=False).index)
        shown = known[known[f"{treatment}_responsive"]].head(10)
        print(f"Stärkste bekannte ROS-Gene ({treatment}):" + ("" if len(shown) else " keine responsiv."))
        if len(shown):
            print(shown[["symbol", "category", f"{treatment}_peak_log2fc", f"{treatment}_peak_time_min",
                     f"{treatment}_pattern"]].round(2).to_string(index=False))

    # --- 2. Unbiased Discovery: stärkste Antworter, unabhängig von der Liste ---
    top_frames = []
    for treatment in treatments:
        top = ros_analysis.top_responsive(result, treatment, config.TOP_N_DISCOVERY, patterns)
        if top.empty:
            print(f"\nTOP RESPONSIVE GENES ({treatment}): keine responsiven Gene nach dem Kriterium.")
            continue
        info = ros_analysis.discovery_summary(top, result)
        print(f"\nTOP {info['n_top']} RESPONSIVE GENES ({treatment}): {info['n_in_list']} stehen in der "
              f"ROS-Liste ({info['share_in_top']:.0%}); Liste macht {info['base_rate']:.1%} aller Sonden aus.")
        unexpected = top[~top["in_ros_list"]]
        print(f"Nicht in der Liste (Top 15 von {len(unexpected)}):")
        print(unexpected.head(15)[["rank", "display_name", "systematic_name", f"{treatment}_peak_log2fc",
                                   f"{treatment}_pattern"]].round(2).to_string(index=False))
        top_frames.append(top.assign(treatment=treatment).rename(columns={
            f"{treatment}_direction": "direction", f"{treatment}_peak_log2fc": "peak_log2fc",
            f"{treatment}_peak_time_min": "peak_time_min", f"{treatment}_n_timepoints_above": "n_timepoints_above",
            f"{treatment}_pattern": "pattern"}))
    pd.concat(top_frames).to_csv(config.TABLES_DIR / f"{label}_top_responsive_genes.csv")

    try:
        from src import visualization
        config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        visualization.ros_gene_barplot(ros_table, treatments, ros_genes.CATEGORY_LABELS,
                                       config.MIN_ABS_LOG2FC, config.FIGURES_DIR / f"{label}_ros_genes_barplot.png")
        print(f"\nGrafik gespeichert in {config.FIGURES_DIR}")
    except ImportError:
        logger.warning("matplotlib fehlt – Grafik übersprungen.")
    print(f"Tabellen gespeichert in {config.TABLES_DIR}")
    return prepared.counts


def run_robustness(client: GeoClient, dataset: ImportedDataset, force: bool) -> dict:
    """Phase 6b: Robustheit der responsiven Gene, Stabilität der Muster, Spezifität mit Grauzone."""
    prepared = prepare_responses(client, dataset, force)
    table, treatments, result = prepared.table, prepared.treatments, prepared.result
    label = dataset.label

    scores = robustness.responsive_robustness(table, treatments)
    pattern_scores, pattern_overview = robustness.pattern_stability(table, treatments)
    parts = [result[[c for c in ("display_name", "systematic_name") if c in result.columns]], scores,
             pattern_scores]
    specificity_counts = None
    if len(treatments) == 2:
        specificity_scores, specificity_counts = robustness.specificity_sensitivity(table, *treatments)
        parts += [result[[c for c in result.columns if c.startswith("specificity_")]], specificity_scores]
    combined = pd.concat(parts, axis=1)
    combined.to_csv(config.TABLES_DIR / f"{label}_robustness.csv")
    pattern_overview.to_csv(config.TABLES_DIR / f"{label}_pattern_sensitivity.csv", index=False)

    n_combos = len(config.SENSITIVITY_THRESHOLDS) * len(config.SENSITIVITY_CONSECUTIVE)
    print(f"\nRobustheit responsiver Gene (Anteil von {n_combos} Kombinationen Schwelle x Zeitpunkte-in-Folge; "
          f"robust >= {config.ROBUSTNESS_ROBUST_MIN}, fragil <= {config.ROBUSTNESS_FRAGILE_MAX}):")
    for treatment in treatments:
        counts = combined[f"{treatment}_robustness_class"].value_counts()
        print(f"  {treatment}: " + ", ".join(f"{name}: {n}" for name, n in counts.items()))

    print("\nStabilität der Zeitverlaufs-Muster (Anteil der Kombinationen aus Beginn-Grenze x "
          "Rückkehr-Anteil mit gleichem Muster, nur responsive Gene):")
    for treatment in treatments:
        values = combined[f"{treatment}_pattern_stability"].dropna()
        shares = pattern_overview[pattern_overview["treatment"] == treatment]["share_transient"]
        print(f"  {treatment}: Median {values.median():.2f}; {(values >= 0.75).mean():.0%} der Gene >= 0.75 stabil. "
              f"Anteil 'transient' schwankt je nach Parametern zwischen {shares.min():.0%} und {shares.max():.0%}.")

    if specificity_counts is not None:
        default_labels = result[[c for c in result.columns if c.startswith("specificity_")][0]]
        print("\nSpezifität: Gene je Label und wie stabil das Label über Schwelle x Grauzone ist:")
        summary = pd.DataFrame({"n_genes": default_labels.value_counts(),
                                "median_stability": specificity_scores.groupby(default_labels).median(),
                                "share_stable_>=0.75": specificity_scores.groupby(default_labels)
                                .apply(lambda s: (s >= 0.75).mean())}).round(2)
        print(summary.to_string())
        specificity_counts.to_csv(config.TABLES_DIR / f"{label}_specificity_sensitivity.csv", index=False)

    print(f"\nTabellen gespeichert in {config.TABLES_DIR}")
    return prepared.counts


def run_validation(client: GeoClient, dataset: ImportedDataset, force: bool) -> dict:
    """Phase 6b: Beispielrechnung, GO-Abgleich der Genliste, Positivkontrollen."""
    prepared = prepare_responses(client, dataset, force)
    table, treatments, result = prepared.table, prepared.treatments, prepared.result
    label = dataset.label
    lead = "H2O2" if "H2O2" in treatments else treatments[0]

    # 1. Beispielrechnung für das stärkste Gen
    probe = result[f"{lead}_peak_log2fc"].abs().idxmax()
    time_min = float(result.loc[probe, f"{lead}_peak_time_min"])
    example = provenance.worked_example(dataset.expression, dataset.design, lead, probe, time_min,
                                        str(result.loc[probe, "display_name"]))
    print("\n" + "\n".join(example["lines"]))
    table_value = float(table.loc[probe, f"log2fc_{lead}_{time_min:g}min"])
    print(f"  Wert in der log2FC-Tabelle: {table_value:.4f} -> "
          + ("identisch mit M" if abs(table_value - example["M"]) < 1e-9 else "ABWEICHUNG!"))

    # 2. GO-Abgleich der Genliste
    annotation_table = load_annotation(client, dataset.samples, force)
    support = ros_evidence.go_support(annotation_table)
    support.to_csv(config.TABLES_DIR / f"{label}_ros_list_go_support.csv", index=False)
    print("\nGO-Abgleich der kuratierten Genliste (Stichwortsuche in den GO-Begriffen der Annotation):")
    print(support["go_support"].value_counts().to_string())
    unsupported = support[support["go_support"] != "supported"]
    if len(unsupported):
        print("Bitte manuell prüfen (SGD): " + ", ".join(
            f"{r.symbol} ({r.go_support})" for r in unsupported.itertuples()))

    # 3. Positivkontrollen
    patterns = timecourse.classify_patterns_all(table, treatments)
    ros_table, _ = ros_genes.build_ros_table(result, patterns)
    controls = ros_evidence.validate_positive_controls(ros_table, lead)
    controls.to_csv(config.TABLES_DIR / f"{label}_positive_controls.csv", index=False)
    print(f"\nPositivkontrollen ({lead}, erwartete Richtung vs. beobachtet):")
    print(controls.to_string(index=False))
    agree = int((controls["status"] == "agrees").sum())
    print(f"{agree} von {len(controls)} Positivkontrollen stimmen überein. "
          "Die Erwartungen sind unbelegt, bis du in der Spalte 'source' Quellen einträgst "
          "(ros_evidence.POSITIVE_CONTROLS).")
    print(f"\nTabellen gespeichert in {config.TABLES_DIR}")
    return prepared.counts


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        geo_id = normalize_geo_id(args.dataset)
        client = GeoClient(config.RAW_DIR)
        datasets = run_import(client, geo_id, args)
        runners = {"differential": run_differential, "timecourse": run_timecourse, "ros": run_ros,
                   "robustness": run_robustness, "validation": run_validation}
        for dataset in datasets:
            if args.inspect_annotation:
                try:
                    inspect_annotation(client, dataset, args.force_download)
                except (annotation.AnnotationError, GeoError) as error:
                    logger.error("Annotation nicht verfügbar: %s", error)
            if args.analysis in AVAILABLE_ANALYSES:
                started = time.time()
                if args.analysis == "log2fc":
                    run_log2fc(client, geo_id, dataset, args.force_download)
                    counts = {}
                else:
                    counts = runners[args.analysis](client, dataset, args.force_download)
                finish_run(client, geo_id, dataset, args.analysis, counts, started)
    except (GeoError, metadata.MetadataError, statistics.StatisticsError,
            differential_expression.DifferentialExpressionError, ros_genes.RosGeneError,
            annotation.AnnotationError) as error:
        logger.error("%s", error)
        return 1

    if args.analysis in PLANNED_ANALYSES or args.report:
        print("\nHinweis: Diese Analyse bzw. der Report ist noch nicht implementiert "
              "(kommt in einer späteren Phase).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
