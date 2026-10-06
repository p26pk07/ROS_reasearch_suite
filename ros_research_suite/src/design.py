"""Experimentelles Design: Welches Sample ist welche Behandlung zu welchem Zeitpunkt?

Aus der rohen Sample-Tabelle (src/metadata.py) wird eine saubere Design-Tabelle:

    treatment | concentration | time_min | is_baseline | replicate | needs_review ...

Grundprinzip: Das Programm RÄT nicht stillschweigend. Es zeigt die erkannte
Zuordnung an, und du bestätigst sie (oder korrigierst sie über eine CSV-Datei).
Erst die bestätigte Zuordnung wird gespeichert und später für Analysen genutzt.
"""
import logging
import re
from pathlib import Path
from typing import Callable

import pandas as pd

from src.metadata import MetadataError

logger = logging.getLogger(__name__)


class DesignNotConfirmedError(MetadataError):
    """Die Zuordnung ist unklar oder wurde vom Benutzer nicht bestätigt."""


# --- Mustererkennung ----------------------------------------------------------
# Zeit mit Einheit, z. B. "30 min", "t=30min", "2 h". Der Lookbehind verhindert,
# dass Ziffern mitten in Wörtern/Zahlen (H2O2, 0.3mM) getroffen werden.
_TIME_WITH_UNIT = re.compile(
    r"(?<![A-Za-z0-9.])(\d+(?:[.,]\d+)?)\s*(min(?:ute)?s?|h(?:rs?|ours?)?)\b", re.I)
# Null ist in jeder Einheit gleich, deshalb darf "t=0" / "time 0" auch ohne Einheit gelten.
_TIME_ZERO = re.compile(r"\b(?:t|time)\s*=?\s*0\b", re.I)

_TREATMENT_PATTERNS = {
    "H2O2": re.compile(r"h2o2|hydrogen\s*peroxide", re.I),
    "MMS": re.compile(r"\bmms\b|methyl\s*methane\s*sulfonate", re.I),
    "control": re.compile(r"\bcontrol\b|\bctrl\b|untreated|unstressed|\bmock\b", re.I),
}
_CONCENTRATION = re.compile(r"(\d+(?:[.,]\d+)?)\s*(mM|µM|uM|nM|%)")
_STRAIN = re.compile(r"strain\s*[:=]?\s*([\w\-]+)", re.I)

# Reihenfolge = Priorität: spezifischere Felder zuerst. 'treatment_protocol_ch1'
# fehlt absichtlich, weil dort bei GSE12220 in JEDEM Sample beide Behandlungen stehen.
_TREATMENT_FIELDS = ("source_name_ch1", "title", "description")
_TIME_FIELDS = ("title", "source_name_ch1")


def extract_time_minutes(text: str) -> float | None:
    """Findet eine Zeitangabe und rechnet in Minuten um.

    '30 min' -> 30.0 | '2 h' -> 120.0 | 't=0' -> 0.0 | kein Treffer -> None
    """
    match = _TIME_WITH_UNIT.search(text)
    if match:
        value = float(match.group(1).replace(",", "."))
        return value * 60 if match.group(2).lower().startswith("h") else value
    if _TIME_ZERO.search(text):
        return 0.0
    return None


def detect_treatment(row: pd.Series) -> tuple[str, str]:
    """Bestimmt die Behandlung eines Samples.

    Das erste Feld (nach Priorität), in dem GENAU EINE Behandlung vorkommt, gewinnt.
    Kommen in einem Feld mehrere vor, bleibt es 'unknown' (-> manuelle Prüfung).

    Returns:
        (Behandlung, Feld, in dem sie gefunden wurde)
    """
    for field in _TREATMENT_FIELDS:
        text = str(row.get(field, ""))
        found = {name for name, pattern in _TREATMENT_PATTERNS.items() if pattern.search(text)}
        if len(found) == 1:
            return found.pop(), field
        if len(found) > 1:
            return "unknown", field
    return "unknown", "-"


def detect_time(row: pd.Series) -> float | None:
    """Zeitpunkt in Minuten. Widersprechen sich Titel und Quelle, wird nichts geraten."""
    times = {extract_time_minutes(str(row.get(field, ""))) for field in _TIME_FIELDS}
    times.discard(None)
    return times.pop() if len(times) == 1 else None


def detect_concentration(row: pd.Series) -> str:
    match = _CONCENTRATION.search(str(row.get("description", "")))
    return f"{match.group(1)} {match.group(2)}" if match else ""


def detect_strain(row: pd.Series) -> str:
    if "strain" in row.index and str(row["strain"]) not in ("", "nan"):
        return str(row["strain"])
    for column in row.index:
        if "characteristics" in str(column):
            match = _STRAIN.search(str(row[column]))
            if match:
                return match.group(1)
    return ""


def _finalize(design: pd.DataFrame) -> pd.DataFrame:
    """Berechnet abgeleitete Spalten (Baseline, Replikatnummer, Prüfbedarf)."""
    design["is_baseline"] = design["time_min"] == 0
    design["replicate"] = design.groupby(
        ["treatment", "time_min"], dropna=False).cumcount() + 1
    design["needs_review"] = (design["treatment"] == "unknown") | design["time_min"].isna()
    return design


def build_design(sample_table: pd.DataFrame) -> pd.DataFrame:
    """Erzeugt die (noch unbestätigte) Design-Tabelle aus der Sample-Tabelle."""
    records = []
    for sample_id, row in sample_table.iterrows():
        treatment, evidence = detect_treatment(row)
        records.append({
            "sample_id": sample_id,
            "organism": str(row.get("organism_ch1", "")),
            "strain": detect_strain(row),
            "treatment": treatment,
            "concentration": detect_concentration(row),
            "time_min": detect_time(row),
            "evidence": evidence,
        })
    design = pd.DataFrame(records).set_index("sample_id")
    design["time_min"] = pd.to_numeric(design["time_min"], errors="coerce")
    return _finalize(design)


def apply_overrides(design: pd.DataFrame, overrides: pd.DataFrame) -> pd.DataFrame:
    """Überschreibt Behandlung/Zeit mit manuellen Angaben (Spalten: sample_id, treatment, time_min)."""
    design = design.copy()
    overrides = overrides.set_index("sample_id")
    unknown_ids = [s for s in overrides.index if s not in design.index]
    if unknown_ids:
        raise MetadataError(f"Override enthält unbekannte Samples: {unknown_ids}")
    for sample_id, row in overrides.iterrows():
        if "treatment" in row and pd.notna(row["treatment"]):
            design.loc[sample_id, "treatment"] = str(row["treatment"])
        if "time_min" in row and pd.notna(row["time_min"]):
            design.loc[sample_id, "time_min"] = float(row["time_min"])
        design.loc[sample_id, "evidence"] = "manual"
    return _finalize(design)


def print_design_tree(design: pd.DataFrame) -> None:
    """Zeigt die erkannte Zuordnung als Baum, z. B.

        H2O2 (0.3 mM, strain Y262)
        ├── 0 min [Baseline]: GSM307217
        └── 30 min: GSM307218
    """
    for treatment, group in design.groupby("treatment", sort=False):
        first = group.iloc[0]
        details = ", ".join(x for x in (first["concentration"],
                                        f"strain {first['strain']}" if first["strain"] else "") if x)
        title = "UNKLAR – bitte prüfen" if treatment == "unknown" else treatment
        print(f"\n{title}" + (f" ({details})" if details else ""))
        ordered = group.sort_values("time_min", na_position="last")
        for i, (sample_id, row) in enumerate(ordered.iterrows()):
            branch = "└──" if i == len(ordered) - 1 else "├──"
            time = "Zeit unklar" if pd.isna(row["time_min"]) else f"{row['time_min']:g} min"
            baseline = " [Baseline]" if row["is_baseline"] else ""
            print(f"{branch} {time}{baseline}: {sample_id}")
    if design.groupby(["treatment", "time_min"], dropna=False).size().max() == 1:
        print("\nHinweis: Pro Behandlung und Zeitpunkt gibt es nur 1 Sample (keine Replikate). "
              "Klassische p-Werte sind damit nicht berechenbar.")


def write_override_template(design: pd.DataFrame, path: Path) -> None:
    """Schreibt die aktuelle Vermutung als editierbare CSV (nur, wenn sie noch nicht existiert)."""
    if not path.exists():
        design.reset_index()[["sample_id", "treatment", "time_min"]].to_csv(path, index=False)


def confirm_design(
    design: pd.DataFrame,
    folder: Path,
    label: str,
    interactive: bool = True,
    assume_yes: bool = False,
    reconfirm: bool = False,
    input_func: Callable[[str], str] = input,
) -> pd.DataFrame:
    """Gibt die vom Benutzer bestätigte Design-Tabelle zurück.

    Ablauf:
      1. Existiert schon eine bestätigte Datei -> wird wiederverwendet.
      2. Existiert <label>_design_overrides.csv -> manuelle Korrekturen anwenden.
      3. Zuordnung anzeigen. Ist etwas unklar -> Abbruch mit Anleitung.
      4. Benutzer bestätigt -> Datei <label>_design_confirmed.csv wird gespeichert.
    """
    folder.mkdir(parents=True, exist_ok=True)
    confirmed_path = folder / f"{label}_design_confirmed.csv"
    overrides_path = folder / f"{label}_design_overrides.csv"

    if confirmed_path.exists() and not reconfirm:
        logger.info("Bestätigte Zuordnung gefunden: %s", confirmed_path.name)
        return pd.read_csv(confirmed_path, index_col="sample_id", keep_default_na=False,
                           na_values=[""], dtype={"strain": str, "concentration": str}
                           ).fillna({"strain": "", "concentration": ""})

    if overrides_path.exists():
        design = apply_overrides(design, pd.read_csv(overrides_path))
        logger.info("Manuelle Korrekturen angewendet: %s", overrides_path.name)

    print_design_tree(design)

    def refuse(reason: str) -> DesignNotConfirmedError:
        write_override_template(design, overrides_path)
        return DesignNotConfirmedError(
            f"{reason}\nKorrigiere die Datei:\n  {overrides_path}\n"
            "(Spalten treatment und time_min anpassen, speichern) und starte das Programm erneut.")

    if design["needs_review"].any():
        raise refuse("Bei manchen Samples ist Behandlung oder Zeit unklar.")

    if not assume_yes:
        if not interactive:
            raise DesignNotConfirmedError(
                "Zuordnung nicht bestätigt. Starte interaktiv im Terminal oder nutze --yes.")
        answer = input_func("\nIst diese Zuordnung korrekt? [j/n]: ").strip().lower()
        if answer not in ("j", "ja", "y", "yes"):
            raise refuse("Zuordnung abgelehnt.")

    design.to_csv(confirmed_path)
    logger.info("Zuordnung gespeichert: %s", confirmed_path.name)
    return design
