"""Zeitverlauf-Muster der Genantwort (Phase 5).

Aus den log2FC-Werten pro Zeitpunkt (gegen t=0) wird je Gen und Behandlung bestimmt:

  Richtung   up / down / both (erst hoch, dann runter oder umgekehrt)
  Beginn     'early' oder 'late'  (Beginn der Antwort <= EARLY_ONSET_MAX_MIN?)
  Dauer      'sustained'        am LETZTEN Zeitpunkt noch über der Schwelle
             'transient'        am Ende unter der Schwelle UND <= TRANSIENT_RETURN_FRACTION
                                des Maximums (= Rückkehr Richtung Ausgangswert)
             'partial_decline'  am Ende unter der Schwelle, aber noch über diesem Bruchteil
             'biphasic'         Gen ist sowohl hoch- als auch herunterreguliert

Muster-Name: '<Beginn>_<Dauer>_<Richtung>', z. B. 'early_transient_up'.
Das sind Kriterien für Beschreibung, keine Kategorien mit statistischer Absicherung
(keine Replikate). Alle Schwellen stehen in config.py.
"""
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

import config
from src import differential_expression as de


def _first_run(mask: np.ndarray, length: int) -> tuple[np.ndarray, np.ndarray]:
    """Findet pro Gen die erste Folge von `length` aufeinanderfolgenden True-Werten.

    Returns: (hat eine Folge?, Index des Zeitpunkts, an dem die erste Folge beginnt).
    Trick: sliding_window_view legt ein Fenster der Breite `length` über jede Zeile;
    ein Fenster ist 'voll', wenn alle Werte darin True sind.
    """
    n_genes, n_points = mask.shape
    if length > n_points:
        return np.zeros(n_genes, dtype=bool), np.zeros(n_genes, dtype=int)
    full_windows = sliding_window_view(mask, length, axis=1).all(axis=-1)
    return full_windows.any(axis=1), full_windows.argmax(axis=1)  # argmax = erstes True


def classify_patterns(table: pd.DataFrame, treatment: str,
                      min_abs_log2fc: float | None = None, min_consecutive: int | None = None,
                      early_max_min: float | None = None,
                      return_fraction: float | None = None) -> pd.DataFrame:
    """Bestimmt das Zeitverlaufsmuster jedes Gens für eine Behandlung."""
    threshold = config.MIN_ABS_LOG2FC if min_abs_log2fc is None else min_abs_log2fc
    needed = config.MIN_CONSECUTIVE_TIMEPOINTS if min_consecutive is None else min_consecutive
    early_max = config.EARLY_ONSET_MAX_MIN if early_max_min is None else early_max_min
    fraction = config.TRANSIENT_RETURN_FRACTION if return_fraction is None else return_fraction

    columns = de.timepoint_columns(table, treatment)
    if not columns:
        raise de.DifferentialExpressionError(f"Keine log2FC-Spalten für '{treatment}' gefunden.")
    times = np.array([time for time, _ in columns])
    values = table[[name for _, name in columns]].to_numpy(dtype=float)

    has_up, onset_up = _first_run(values >= threshold, needed)
    has_down, onset_down = _first_run(values <= -threshold, needed)
    responsive = has_up | has_down
    both = has_up & has_down
    only_down = has_down & ~has_up

    # Für Gene mit EINER Richtung wird so gespiegelt, dass "Antwort" immer positiv ist.
    sign = np.where(only_down, -1.0, 1.0)
    mirrored = values * sign[:, None]
    final = mirrored[:, -1]
    peak = np.nanmax(mirrored, axis=1)

    onset_index = np.where(both, np.minimum(onset_up, onset_down),
                           np.where(only_down, onset_down, onset_up))
    onset_time = np.where(responsive, times[onset_index], np.nan)

    persistence = np.select(
        [~responsive, both, final >= threshold, final <= fraction * peak],
        ["none", "biphasic", "sustained", "transient"], default="partial_decline")
    onset_class = np.where(~responsive, "none", np.where(onset_time <= early_max, "early", "late"))
    direction = np.select([both, has_up, has_down], ["both", "up", "down"], default="none")

    pattern = [
        "none" if not r else "biphasic" if b else f"{o}_{p}_{d}"
        for r, b, o, p, d in zip(responsive, both, onset_class, persistence, direction)]

    prefix = f"{treatment}_"
    return pd.DataFrame({
        prefix + "pattern": pattern,
        prefix + "onset_class": onset_class,
        prefix + "persistence": persistence,
        prefix + "onset_time_min": onset_time,
        prefix + "final_log2fc": values[:, -1],
    }, index=table.index)


def classify_patterns_all(table: pd.DataFrame, treatments: list[str]) -> pd.DataFrame:
    """Muster für alle Behandlungen, mit Gennamen-Spalten vorne."""
    annotation = [c for c in de.ANNOTATION_COLUMNS if c in table.columns]
    parts = [table[annotation]] + [classify_patterns(table, t) for t in treatments]
    return pd.concat(parts, axis=1)
