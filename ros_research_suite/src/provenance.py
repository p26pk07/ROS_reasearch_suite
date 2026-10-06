"""Herkunfts- und Reproduzierbarkeitsangaben (wissenschaftlicher Standard 2/3/4).

- analysis_metadata: JSON mit Datensatz, Samples, Datenebene, Filtern, Parametern, Methoden, Versionen.
  Wird bei jedem Analyselauf ergänzt (Liste 'analyses_run').
- worked_example: vollständige Beispielrechnung für ein Gen (M, A, log2FC) mit Selbstprüfung.
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import config
from src import statistics

PARAMETER_NAMES = (
    "MIN_ABS_LOG2FC", "MIN_CONSECUTIVE_TIMEPOINTS", "EARLY_ONSET_MAX_MIN", "TRANSIENT_RETURN_FRACTION",
    "SPECIFICITY_LOW_FRACTION", "SENSITIVITY_THRESHOLDS", "SENSITIVITY_CONSECUTIVE",
    "ROBUSTNESS_ROBUST_MIN", "ROBUSTNESS_FRAGILE_MAX", "PATTERN_EARLY_GRID", "PATTERN_RETURN_GRID",
    "ORGANISM_FILTER", "CONTROL_PROBE_PREFIX", "MIN_ANNOTATION_MATCH_RATE", "TOP_N_GENES", "TOP_N_DISCOVERY",
)

METHODS = {
    "log2fc": "log2FC = E_behandelt - E_Baseline (E = log2-RMA-Wert, Baseline = t=0 derselben Behandlung)",
    "differential": "responsive = |log2FC| >= Schwelle in >= k aufeinanderfolgenden Zeitpunkten; "
                    "kein Signifikanztest (keine Replikate)",
    "timecourse": "Muster aus Beginn (early/late), Dauer (sustained/transient/partial_decline/biphasic)",
    "ros": "kuratierte Genliste (Vorwissen) gegen unvoreingenommene Top-Antworter; deskriptiver Vergleich",
    "robustness": "Anteil der Parameterkombinationen, in denen Gen/Muster/Spezifität gleich bleibt",
    "validation": "Beispielrechnung, GO-Abgleich der Genliste, Positivkontrollen",
}


def package_versions() -> dict[str, str]:
    versions = {"software": config.SOFTWARE_VERSION, "python": sys.version.split()[0],
                "platform": platform.platform()}
    for name in ("pandas", "numpy", "scipy", "matplotlib", "requests"):
        try:
            versions[name] = __import__(name).__version__
        except ImportError:
            versions[name] = "not installed"
    return versions


def _json_safe(value):
    if isinstance(value, (tuple, list)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def build_dataset_info(geo_id: str, series_info: dict, samples: pd.DataFrame, design: pd.DataFrame,
                       series_file: Path, annotation_file: Path | None,
                       annotation_table: pd.DataFrame | None) -> dict:
    """Beschreibung des Datensatzes aus den GEO-Metadaten (nichts aus dem Gedächtnis)."""
    def unique(column):
        return sorted(samples[column].dropna().unique().tolist()) if column in samples.columns else []

    treatments = {}
    for treatment, rows in design.groupby("treatment"):
        treatments[treatment] = {
            "concentration": sorted(rows["concentration"].unique().tolist()),
            "time_points_min": sorted(float(t) for t in rows["time_min"].unique()),
            "samples": {sample: float(row["time_min"]) for sample, row in rows.iterrows()},
            "replicates_per_timepoint": int(rows.groupby("time_min").size().max())}
    annotation = {}
    if annotation_file is not None:
        annotation = {"file": str(annotation_file.name),
                      "retrieved": datetime.fromtimestamp(annotation_file.stat().st_mtime,
                                                          timezone.utc).date().isoformat()}
        if annotation_table is not None and "Annotation Date" in annotation_table.columns:
            annotation["annotation_date_in_file"] = sorted(
                annotation_table["Annotation Date"].replace("", np.nan).dropna().unique().tolist())
    return {
        "geo_accession": geo_id,
        "title": series_info.get("Series_title", ""),
        "pubmed_id_from_geo": series_info.get("Series_pubmed_id", "nicht in GEO-Metadaten angegeben"),
        "geo_submission_date": series_info.get("Series_submission_date", ""),
        "organism": unique("organism_ch1"), "strain": sorted(design["strain"].unique().tolist()),
        "platform": unique("platform_id"), "n_samples": int(len(samples)),
        "treatments": treatments,
        "data_level": "processed (von den Autoren bereitgestellte, RMA-berechnete log2-Werte; "
                      "CEL-Rohdateien werden NICHT verwendet; keine eigene Normalisierung)",
        "data_processing_reported_by_authors": unique("data_processing"),
        "series_matrix_file": series_file.name, "annotation": annotation,
        "data_source": {"series": "NCBI GEO (ftp.ncbi.nlm.nih.gov/geo/series)",
                        "platform_annotation": "NCBI GEO query (www.ncbi.nlm.nih.gov/geo/query/acc.cgi)"},
    }


def update_metadata(path: Path, dataset_info: dict, analysis: str, counts: dict | None = None,
                    outputs: list[str] | None = None) -> Path:
    """Erzeugt oder ergänzt analysis_metadata.json (Liste der bisherigen Analyseläufe bleibt erhalten)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"analyses_run": []}
    metadata["dataset"] = dataset_info
    metadata["parameters"] = {name: _json_safe(getattr(config, name)) for name in PARAMETER_NAMES}
    metadata["statistical_evidence"] = (
        "Gene: nicht testbar (keine Replikate). Keine p-Werte und keine FDR für einzelne Gene. "
        "Sprache: 'responsiv' bezeichnet Effektstärke + Konsistenz, nicht Signifikanz.")
    metadata["software_versions"] = package_versions()
    metadata["analysis_methods"] = METHODS
    metadata["probe_counts"] = {**metadata.get("probe_counts", {}), **(counts or {})}
    metadata["analyses_run"].append({
        "analysis": analysis, "date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "outputs": outputs or []})
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def worked_example(expression: pd.DataFrame, design: pd.DataFrame, treatment: str, probe: str,
                   time_min: float, gene_name: str = "") -> dict:
    """Vollständige Beispielrechnung für eine Sonde, mit Selbstprüfung der Formeln.

    Die Werte kommen aus statistics.ma_values (dieselbe Funktion wie bei MA-Plots); zusätzlich
    wird das Ergebnis über die lineare Skala (2^E) unabhängig nachgerechnet.
    """
    rows = design[design["treatment"] == treatment]
    baseline_ids = list(rows.index[rows["is_baseline"].to_numpy(dtype=bool)])
    treated_ids = list(rows.index[(rows["time_min"] == time_min).to_numpy()
                                  & ~rows["is_baseline"].to_numpy(dtype=bool)])
    e_base = float(expression.loc[probe, baseline_ids].mean())
    e_treated = float(expression.loc[probe, treated_ids].mean())
    ma = statistics.ma_values(expression.loc[[probe]], design, treatment, time_min).loc[probe]
    linear_ratio = (2 ** e_treated) / (2 ** e_base)
    check_m = bool(np.isclose(np.log2(linear_ratio), ma["M"]))
    check_a = bool(np.isclose(ma["A"], (e_treated + e_base) / 2))
    lines = [
        f"Beispielrechnung: {gene_name or probe} (Sonde {probe}), {treatment}, {time_min:g} min vs. t=0",
        f"  E_Baseline  (log2, {', '.join(baseline_ids)}) = {e_base:.4f}",
        f"  E_behandelt (log2, {', '.join(treated_ids)}) = {e_treated:.4f}",
        f"  M = E_behandelt - E_Baseline = {e_treated:.4f} - {e_base:.4f} = {ma['M']:.4f}   (= log2FC)",
        f"  A = (E_behandelt + E_Baseline) / 2 = {ma['A']:.4f}",
        f"  Gegenprobe (lineare Skala): 2^{e_treated:.4f} / 2^{e_base:.4f} = {linear_ratio:.2f}"
        f"  -> log2 = {np.log2(linear_ratio):.4f}   ({'stimmt' if check_m else 'ABWEICHUNG'})",
        f"  Das ist eine {linear_ratio:.1f}-fache Änderung der RMA-Signalintensität "
        "(Genexpression, keine ROS-Konzentration).",
    ]
    return {"M": float(ma["M"]), "A": float(ma["A"]), "linear_ratio": float(linear_ratio),
            "check_m": check_m, "check_a": check_a, "lines": lines}
