"""Metadaten- und Datenimport aus GEO-Series-Matrix-Dateien.

Aufbau einer Series-Matrix-Datei (Textdatei, gzip-komprimiert):

    !Series_title        "..."          <- Infos zum ganzen Datensatz
    !Sample_title        "s1"  "s2" ... <- eine Zeile pro Eigenschaft, eine Spalte pro Sample
    !Sample_characteristics_ch1 "time: 30 min" ...
    !series_matrix_table_begin
    "ID_REF"  "GSM1"  "GSM2" ...        <- Expressionstabelle
    ...
    !series_matrix_table_end

Dieses Modul liest nur Dateien ein. Die Interpretation (Behandlung, Zeitpunkt,
Kontrolle) passiert in src/design.py.
"""
import csv
import gzip
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

TABLE_BEGIN = "!series_matrix_table_begin"
TABLE_END = "!series_matrix_table_end"


class MetadataError(Exception):
    """Datei ist keine lesbare Series-Matrix-Datei."""


def _open_text(path: Path):
    """Öffnet .gz und normale Textdateien gleich."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _split_line(line: str) -> list[str]:
    """Zerlegt eine tab-getrennte Zeile und entfernt die Anführungszeichen."""
    return [cell.strip().strip('"') for cell in line.rstrip("\n").split("\t")]


def read_header(path: Path) -> tuple[dict[str, str], list[tuple[str, list[str]]]]:
    """Liest den Kopf der Datei (alles vor der Expressionstabelle).

    Returns:
        series_info: z. B. {'Series_title': '...'}
        sample_rows: Liste (Schlüssel, Werte pro Sample), Reihenfolge bleibt erhalten.
            Eine Liste statt dict, weil Schlüssel wie 'Sample_characteristics_ch1'
            mehrfach vorkommen.
    """
    series_info: dict[str, str] = {}
    sample_rows: list[tuple[str, list[str]]] = []
    found_marker = False
    with _open_text(path) as file:
        for line in file:
            if line.startswith(TABLE_BEGIN):
                found_marker = True
                break
            if not line.startswith("!"):
                continue
            cells = _split_line(line)
            key = cells[0].lstrip("!")
            if key.startswith("Series_"):
                series_info.setdefault(key, " | ".join(cells[1:]))
            elif key.startswith("Sample_"):
                sample_rows.append((key, cells[1:]))
    if not found_marker:
        raise MetadataError(
            f"{path.name}: Keine Expressionstabelle gefunden "
            f"('{TABLE_BEGIN}' fehlt). Falsches Dateiformat?"
        )
    return series_info, sample_rows


def build_sample_table(sample_rows: list[tuple[str, list[str]]]) -> pd.DataFrame:
    """Baut aus den Sample-Zeilen eine Tabelle mit EINER Zeile pro Sample.

    'Sample_characteristics_ch1'-Zeilen der Form 'key: value' werden zu
    eigenen Spalten (z. B. 'time', 'strain'). Zeilen ohne ':' bekommen
    eine laufende Nummer als Spaltennamen.
    """
    if not sample_rows:
        raise MetadataError("Keine Sample-Metadaten (!Sample_...) gefunden.")

    columns: dict[str, list[str]] = {}
    unnamed_counter = 0
    for key, values in sample_rows:
        if key == "Sample_characteristics_ch1":
            key_parts = [v.split(":", 1)[0].strip() for v in values if ":" in v]
            if key_parts and len(set(key_parts)) == 1:
                column = key_parts[0]
                columns[column] = [v.split(":", 1)[1].strip() if ":" in v else v
                                   for v in values]
                continue
            unnamed_counter += 1
            columns[f"characteristics_{unnamed_counter}"] = values
        else:
            name = key.removeprefix("Sample_")
            while name in columns:  # doppelte Schlüssel nicht überschreiben
                name += "_2"
            columns[name] = values

    table = pd.DataFrame(columns)
    if "geo_accession" in table.columns:
        table = table.set_index("geo_accession")
    return table


def read_expression_table(path: Path) -> pd.DataFrame:
    """Liest die Expressionstabelle (Zeilen = Sonden/Gene, Spalten = Samples).

    Zeilen, die mit '!' beginnen (Kopf und Tabellenende), werden übersprungen.
    Werte werden zu Zahlen umgewandelt; 'null' und leere Felder werden NaN.
    """
    try:
        table = pd.read_csv(
            path, sep="\t", comment="!", index_col=0,
            quoting=csv.QUOTE_MINIMAL, na_values=["null", "NULL", ""],
            compression="infer",
        )
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as error:
        raise MetadataError(f"{path.name}: Expressionstabelle nicht lesbar ({error})") from error
    if table.empty:
        raise MetadataError(f"{path.name}: Expressionstabelle ist leer.")
    table.index.name = "probe_id"
    return table.apply(pd.to_numeric, errors="coerce")


def align_expression(expression: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Bringt die Spalten der Expressionstabelle in die Reihenfolge der Design-Tabelle.

    Prüft dabei, dass jedes Sample aus dem Design auch eine Expressionsspalte hat.
    Überzählige Spalten werden entfernt (mit Warnung).
    """
    missing = [s for s in design.index if s not in expression.columns]
    if missing:
        raise MetadataError(
            f"Für diese Samples fehlen Expressionsdaten: {', '.join(missing)}")
    extra = [c for c in expression.columns if c not in design.index]
    if extra:
        logger.warning("Spalten ohne Sample-Metadaten werden ignoriert: %s", extra)
    return expression[list(design.index)]


def describe_expression(expression: pd.DataFrame) -> dict[str, float | int | bool]:
    """Kurze Kennzahlen zur Plausibilitätsprüfung der Werte.

    'likely_log_scale' ist nur eine Faustregel: Mikroarray-Werte nach RMA liegen
    typischerweise zwischen ca. 2 und 16 (log2-Skala), Rohintensitäten gehen in
    die Tausende. Sie ersetzt NICHT die Angabe der Autoren ('data_processing').
    """
    values = expression.to_numpy(dtype=float)
    finite = values[~pd.isna(values)]
    maximum = float(finite.max())
    return {
        "n_probes": int(expression.shape[0]),
        "n_samples": int(expression.shape[1]),
        "n_missing": int(pd.isna(values).sum()),
        "min": float(finite.min()),
        "median": float(pd.Series(finite).median()),
        "max": maximum,
        "likely_log_scale": maximum <= 50,
    }
