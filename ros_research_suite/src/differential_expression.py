"""Differential Expression für Zeitreihen OHNE Replikate (z. B. GSE12220).

Ohne Replikate lässt sich die Streuung nicht schätzen, deshalb gibt es hier keine
p-Werte. Stattdessen gilt ein Gen als "responsiv", wenn seine Veränderung gegenüber
t=0 sowohl GROSS (|log2FC| >= Schwelle) als auch REPRODUZIERT über aufeinanderfolgende
Zeitpunkte ist. Das ist ein Kriterium für Effektstärke + Konsistenz, kein
Signifikanztest. Die Schwellen stehen in config.py.
"""
import logging
import re

import numpy as np
import pandas as pd

import config

logger = logging.getLogger(__name__)

_COLUMN_PATTERN = re.compile(r"^log2fc_(?P<treatment>.+)_(?P<time>\d+(?:\.\d+)?)min$")
ANNOTATION_COLUMNS = ("display_name", "gene_symbol", "systematic_name")


class DifferentialExpressionError(Exception):
    """Analyse nicht möglich (z. B. keine Zeitpunkt-Spalten für die Behandlung)."""


def split_control_probes(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Trennt Affymetrix-Kontrollsonden ('AFFX-...') von den echten Sonden ab."""
    is_control = table.index.astype(str).str.startswith(config.CONTROL_PROBE_PREFIX)
    return table[~is_control], table[is_control]


def timepoint_columns(table: pd.DataFrame, treatment: str) -> list[tuple[float, str]]:
    """Alle log2FC-Spalten einer Behandlung als (Zeit in min, Spaltenname), nach ZEIT sortiert.

    Numerisch sortiert, nicht alphabetisch: sonst käme '100min' vor '30min'.
    """
    found = []
    for column in table.columns:
        match = _COLUMN_PATTERN.match(str(column))
        if match and match.group("treatment") == treatment:
            found.append((float(match.group("time")), column))
    return sorted(found)


def longest_true_run(mask: np.ndarray) -> np.ndarray:
    """Längste Folge von True pro Zeile.

    mask: Gene x Zeitpunkte (bool). Beispiel [F, T, T, F, T] -> 2.
    Wir laufen Zeitpunkt für Zeitpunkt durch und zählen mit: bei True +1,
    bei False zurück auf 0; das bisherige Maximum merken wir uns.
    """
    current = np.zeros(mask.shape[0], dtype=int)
    best = np.zeros(mask.shape[0], dtype=int)
    for column in range(mask.shape[1]):
        current = np.where(mask[:, column], current + 1, 0)
        best = np.maximum(best, current)
    return best


def classify_treatment(table: pd.DataFrame, treatment: str,
                       min_abs_log2fc: float | None = None,
                       min_consecutive: int | None = None) -> pd.DataFrame:
    """Bewertet jedes Gen für eine Behandlung.

    Spalten (Präfix '<Behandlung>_'): responsive, up, down, direction ('up'/'down'/
    'both'/'none'), n_timepoints_above, max_log2fc, min_log2fc, peak_log2fc, peak_time_min.
    """
    threshold = config.MIN_ABS_LOG2FC if min_abs_log2fc is None else min_abs_log2fc
    needed = config.MIN_CONSECUTIVE_TIMEPOINTS if min_consecutive is None else min_consecutive

    columns = timepoint_columns(table, treatment)
    if not columns:
        raise DifferentialExpressionError(f"Keine log2FC-Spalten für '{treatment}' gefunden.")
    times = np.array([time for time, _ in columns])
    values = table[[name for _, name in columns]].to_numpy(dtype=float)

    is_up = longest_true_run(values >= threshold) >= needed
    is_down = longest_true_run(values <= -threshold) >= needed
    direction = np.select([is_up & is_down, is_up, is_down], ["both", "up", "down"],
                          default="none")

    abs_values = np.where(np.isnan(values), -np.inf, np.abs(values))
    peak_index = abs_values.argmax(axis=1)
    rows = np.arange(len(table))

    prefix = f"{treatment}_"
    return pd.DataFrame({
        prefix + "responsive": is_up | is_down,
        prefix + "up": is_up,
        prefix + "down": is_down,
        prefix + "direction": direction,
        prefix + "n_timepoints_above": (np.abs(values) >= threshold).sum(axis=1),
        prefix + "max_log2fc": np.nanmax(values, axis=1),
        prefix + "min_log2fc": np.nanmin(values, axis=1),
        prefix + "peak_log2fc": values[rows, peak_index],
        prefix + "peak_time_min": times[peak_index],
    }, index=table.index)


def compare_treatments(result: pd.DataFrame, first: str, second: str) -> pd.Series:
    """Vergleicht zwei Behandlungen pro Gen anhand der Richtung.

    Kategorien: 'both_up', 'both_down', '<first>_only_up', '<second>_only_down', ...,
    'different_direction', 'none'. 'only' heißt: NUR in diesem Datensatz und mit
    diesen Schwellen; es beweist keine Spezifität.
    """
    dir_first = result[f"{first}_direction"]
    dir_second = result[f"{second}_direction"]
    labels = []
    for a, b in zip(dir_first, dir_second):
        if a == "none" and b == "none":
            labels.append("none")
        elif b == "none":
            labels.append(f"{first}_only_{a}")
        elif a == "none":
            labels.append(f"{second}_only_{b}")
        elif a == b and a in ("up", "down"):
            labels.append(f"both_{a}")
        else:
            labels.append("different_direction")
    return pd.Series(labels, index=result.index, name=f"comparison_{first}_vs_{second}")


def classify_all(table: pd.DataFrame, treatments: list[str],
                 min_abs_log2fc: float | None = None,
                 min_consecutive: int | None = None) -> pd.DataFrame:
    """Alle Behandlungen bewerten (+ Vergleich, wenn es genau zwei sind)."""
    annotation = [c for c in ANNOTATION_COLUMNS if c in table.columns]
    parts = [table[annotation]]
    for treatment in treatments:
        parts.append(classify_treatment(table, treatment, min_abs_log2fc, min_consecutive))
    result = pd.concat(parts, axis=1)
    if len(treatments) == 2:
        result = pd.concat([result, compare_treatments(result, *treatments),
                            classify_specificity(table, result, *treatments,
                                                 min_abs_log2fc, min_consecutive)], axis=1)
    return result


def top_genes(result: pd.DataFrame, treatment: str, direction: str, n: int) -> pd.DataFrame:
    """Die n stärksten hoch- ('up') bzw. herunterregulierten ('down') responsiven Gene."""
    if direction == "up":
        return result[result[f"{treatment}_up"]].nlargest(n, f"{treatment}_max_log2fc")
    if direction == "down":
        return result[result[f"{treatment}_down"]].nsmallest(n, f"{treatment}_min_log2fc")
    raise ValueError("direction muss 'up' oder 'down' sein")


def threshold_sensitivity(table: pd.DataFrame, treatments: list[str],
                          thresholds=None, consecutive_values=None) -> pd.DataFrame:
    """Wie stabil ist die Zahl responsiver Gene, wenn man die Kriterien ändert?

    Biologisch: Ohne Replikate hängt jede Trefferliste an der gewählten Schwelle.
    Diese Tabelle macht sichtbar, wie stark. Spalten: threshold, min_consecutive,
    treatment, n_responsive, n_up, n_down.
    """
    thresholds = config.SENSITIVITY_THRESHOLDS if thresholds is None else thresholds
    consecutive_values = (config.SENSITIVITY_CONSECUTIVE if consecutive_values is None
                          else consecutive_values)
    rows = []
    for threshold in thresholds:
        for consecutive in consecutive_values:
            for treatment in treatments:
                result = classify_treatment(table, treatment, threshold, consecutive)
                rows.append({
                    "threshold": threshold, "min_consecutive": consecutive,
                    "treatment": treatment,
                    "n_responsive": int(result[f"{treatment}_responsive"].sum()),
                    "n_up": int(result[f"{treatment}_up"].sum()),
                    "n_down": int(result[f"{treatment}_down"].sum())})
    return pd.DataFrame(rows)


def classify_specificity(table: pd.DataFrame, result: pd.DataFrame, first: str, second: str,
                         min_abs_log2fc: float | None = None, min_consecutive: int | None = None,
                         low_fraction: float | None = None) -> pd.Series:
    """Spezifität zweier Behandlungen mit Grauzone.

    'spezifisch' = in der einen Behandlung responsiv UND in der anderen an allen Zeitpunkten
    |log2FC| < low_fraction * Schwelle. Antwortet die andere Behandlung schwächer, aber nicht
    "gering", ist das Gen 'ambiguous' (wird nicht klassifiziert).
    Labels: shared_up, shared_down, opposite_or_mixed, <first>_specific_<up/down/both>,
    <second>_specific_..., ambiguous, none.
    """
    threshold = config.MIN_ABS_LOG2FC if min_abs_log2fc is None else min_abs_log2fc
    fraction = config.SPECIFICITY_LOW_FRACTION if low_fraction is None else low_fraction
    low = fraction * threshold

    a = classify_treatment(table, first, threshold, min_consecutive)
    b = classify_treatment(table, second, threshold, min_consecutive)
    resp_a, resp_b = a[f"{first}_responsive"], b[f"{second}_responsive"]
    dir_a, dir_b = a[f"{first}_direction"], b[f"{second}_direction"]
    weak_a = np.maximum(a[f"{first}_max_log2fc"].abs(), a[f"{first}_min_log2fc"].abs()) < low
    weak_b = np.maximum(b[f"{second}_max_log2fc"].abs(), b[f"{second}_min_log2fc"].abs()) < low

    labels = []
    for ra, rb, da, db, wa, wb in zip(resp_a, resp_b, dir_a, dir_b, weak_a, weak_b):
        if ra and rb:
            labels.append(f"shared_{da}" if da == db and da in ("up", "down") else "opposite_or_mixed")
        elif ra and wb:
            labels.append(f"{first}_specific_{da}")
        elif rb and wa:
            labels.append(f"{second}_specific_{db}")
        elif ra or rb:
            labels.append("ambiguous")
        else:
            labels.append("none")
    return pd.Series(labels, index=table.index, name=f"specificity_{first}_vs_{second}")
