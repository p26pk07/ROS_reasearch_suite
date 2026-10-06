"""Plattform-Annotation: Sonden-IDs (z. B. '1769308_at') -> Gennamen.

Die Annotationstabelle kommt von GEO (Plattform, z. B. GPL2529). Welche Spalten
sie enthält, hängt von der Plattform ab. Deshalb:
  1. Spalten werden anhand typischer Namen erkannt (oder in config.py festgelegt),
  2. das Programm zeigt, welche Spalten es benutzt,
  3. es prüft, dass die Sonden-IDs wirklich zur Annotation passen.
Es werden keine Gennamen erraten: Was nicht in der Tabelle steht, bleibt leer.
"""
import csv
import gzip
import io
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import config

logger = logging.getLogger(__name__)

TABLE_BEGIN = "!platform_table_begin"
TABLE_END = "!platform_table_end"

SYMBOL_HINTS = ("gene symbol", "gene_symbol", "symbol", "gene name", "gene_name")
SYSTEMATIC_HINTS = ("systematic", "orf", "locus_tag", "locus tag", "sgd")
ORGANISM_HINTS = ("organism",)


class AnnotationError(Exception):
    """Annotation fehlt, ist unlesbar oder passt nicht zu den Daten."""


def parse_platform_table(path: Path) -> pd.DataFrame:
    """Liest die Annotationstabelle (zwischen !platform_table_begin und _end)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    lines: list[str] = []
    inside = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as file:
        for line in file:
            if line.startswith(TABLE_BEGIN):
                inside = True
            elif line.startswith(TABLE_END):
                break
            elif inside:
                lines.append(line)
    if not lines:
        raise AnnotationError(
            f"{path.name}: Keine Annotationstabelle gefunden. "
            "Möglicherweise liefert GEO für diese Plattform keine.")
    table = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str,
                        quoting=csv.QUOTE_NONE, keep_default_na=False,
                        on_bad_lines="warn")
    if "ID" not in table.columns:
        raise AnnotationError(
            f"Annotationstabelle ohne 'ID'-Spalte. Spalten: {list(table.columns)}")
    table = table.set_index("ID")
    return table[~table.index.duplicated(keep="first")]


def _find_column(columns: list[str], hints: tuple[str, ...]) -> str | None:
    """Erst exakte, dann enthaltene Treffer (Groß-/Kleinschreibung egal)."""
    lowered = {column.lower(): column for column in columns}
    for hint in hints:
        if hint in lowered:
            return lowered[hint]
    for hint in hints:
        for low, original in lowered.items():
            if hint in low:
                return original
    return None


def detect_gene_columns(annotation: pd.DataFrame) -> dict[str, str | None]:
    """Findet die Spalten für Gen-Symbol, systematischen Namen (z. B. YGR088W) und Organismus."""
    columns = list(annotation.columns)
    chosen = {
        "gene_symbol": config.ANNOTATION_SYMBOL_COLUMN or _find_column(columns, SYMBOL_HINTS),
        "systematic_name": (config.ANNOTATION_SYSTEMATIC_COLUMN
                            or _find_column(columns, SYSTEMATIC_HINTS)),
        "organism": _find_column(columns, ORGANISM_HINTS),
    }
    for source in chosen.values():
        if source is not None and source not in columns:
            raise AnnotationError(
                f"Spalte '{source}' (aus config.py) gibt es nicht. Vorhanden: {columns}")
    return chosen


def probe_match_rate(probe_ids: pd.Index, annotation: pd.DataFrame) -> float:
    """Anteil der Sonden, die in der Annotation vorkommen (0 bis 1)."""
    return float(np.mean([probe in annotation.index for probe in probe_ids]))


def attach_annotation(table: pd.DataFrame, annotation: pd.DataFrame,
                      columns: dict[str, str | None]) -> pd.DataFrame:
    """Stellt Gen-Spalten vor die Ergebnistabelle. Fehlende Einträge bleiben leer (NaN)."""
    selected = {}
    for new_name, source in columns.items():
        if source is not None:
            values = annotation[source].reindex(table.index)
            selected[new_name] = values.replace(["", "---"], np.nan)
    if not selected:
        return table
    annotated = pd.DataFrame(selected, index=table.index)
    if "gene_symbol" in annotated and "systematic_name" in annotated:
        # Anzeigename: Symbol, sonst der systematische Name (viele ORFs haben kein Symbol)
        annotated.insert(0, "display_name",
                         annotated["gene_symbol"].fillna(annotated["systematic_name"]))
    return pd.concat([annotated, table], axis=1)


def filter_by_organism(table: pd.DataFrame, organism: str) -> tuple[pd.DataFrame, int]:
    """Entfernt Sonden, deren Organismus bekannt ist und nicht dem gewünschten entspricht.

    Sonden ohne Organismus-Angabe bleiben erhalten (es wird nichts geraten).
    Returns: (gefilterte Tabelle, Anzahl entfernter Sonden)
    """
    if "organism" not in table.columns:
        return table, 0
    keep = table["organism"].isna() | (table["organism"] == organism)
    return table[keep], int((~keep).sum())
