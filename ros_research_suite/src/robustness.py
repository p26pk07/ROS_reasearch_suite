"""Robustheits- und Sensitivitätsanalyse (wissenschaftlicher Standard 14/15).

Frage: Bleiben Ergebnisse bestehen, wenn man die (teils konventionellen) Parameter ändert?
Wenn ein Gen oder Muster nur bei EINER Einstellung auftaucht, wird es als fragil gekennzeichnet.
"""
import numpy as np
import pandas as pd

import config
from src import differential_expression as de
from src import timecourse as tc


def responsive_robustness(table: pd.DataFrame, treatments: list[str]) -> pd.DataFrame:
    """Score je Gen und Behandlung: Anteil der Kombinationen (Schwelle x Zeitpunkte in Folge),
    in denen das Gen responsiv ist.

    Klassen (nur für Gene, die mit den Standardwerten responsiv sind):
    robust (>= ROBUSTNESS_ROBUST_MIN), fragile (<= ROBUSTNESS_FRAGILE_MAX), sonst intermediate;
    nicht responsiv bei den Standardwerten -> 'not_responsive'.
    """
    parts = {}
    for treatment in treatments:
        default = de.classify_treatment(table, treatment)[f"{treatment}_responsive"].to_numpy()
        hits = np.zeros(len(table))
        combos = 0
        for threshold in config.SENSITIVITY_THRESHOLDS:
            for consecutive in config.SENSITIVITY_CONSECUTIVE:
                result = de.classify_treatment(table, treatment, threshold, consecutive)
                hits += result[f"{treatment}_responsive"].to_numpy()
                combos += 1
        score = hits / combos
        classes = np.select(
            [~default, score >= config.ROBUSTNESS_ROBUST_MIN, score <= config.ROBUSTNESS_FRAGILE_MAX],
            ["not_responsive", "robust", "fragile"], default="intermediate")
        parts[f"{treatment}_robustness_score"] = np.round(score, 3)
        parts[f"{treatment}_robustness_class"] = classes
    return pd.DataFrame(parts, index=table.index)


def pattern_stability(table: pd.DataFrame, treatments: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wie stabil sind die Zeitverlaufs-Muster gegenüber den zwei Konventions-Parametern?

    Returns:
        (Stabilität je Gen: Anteil der Kombinationen, in denen das Muster dem Standard-Muster
         entspricht; Übersicht: Anzahl Gene je Dauerklasse für jede Parameterkombination)
    """
    stability, overview = {}, []
    for treatment in treatments:
        default = tc.classify_patterns(table, treatment)[f"{treatment}_pattern"].to_numpy()
        same = np.zeros(len(table))
        combos = 0
        for early in config.PATTERN_EARLY_GRID:
            for fraction in config.PATTERN_RETURN_GRID:
                patterns = tc.classify_patterns(table, treatment, early_max_min=early,
                                                return_fraction=fraction)
                same += (patterns[f"{treatment}_pattern"].to_numpy() == default)
                combos += 1
                persistence = patterns[f"{treatment}_persistence"]
                n_responsive = int((persistence != "none").sum())
                overview.append({
                    "treatment": treatment, "early_max_min": early, "return_fraction": fraction,
                    "n_responsive": n_responsive,
                    "n_transient": int((persistence == "transient").sum()),
                    "n_sustained": int((persistence == "sustained").sum()),
                    "share_transient": round((persistence == "transient").sum() / n_responsive, 3)
                    if n_responsive else float("nan")})
        score = same / combos
        stability[f"{treatment}_pattern_stability"] = np.where(default == "none", np.nan, np.round(score, 3))
    return pd.DataFrame(stability, index=table.index), pd.DataFrame(overview)


def specificity_sensitivity(table: pd.DataFrame, first: str, second: str
                            ) -> tuple[pd.Series, pd.DataFrame]:
    """Wie stabil ist die Spezifitäts-Einteilung (Schwelle x Grauzonen-Breite)?

    Returns:
        (Anteil der Kombinationen mit demselben Label wie bei den Standardwerten je Gen,
         Anzahl Gene je Label und Kombination)
    """
    result = de.classify_all(table, [first, second])
    default = result[f"specificity_{first}_vs_{second}"].to_numpy()
    same = np.zeros(len(table))
    combos, counts = 0, []
    for threshold in config.SPECIFICITY_THRESHOLD_GRID:
        for fraction in config.SPECIFICITY_LOW_FRACTION_GRID:
            labels = de.classify_specificity(table, result, first, second, threshold, None, fraction)
            same += (labels.to_numpy() == default)
            combos += 1
            for label, n in labels.value_counts().items():
                counts.append({"threshold": threshold, "low_fraction": fraction,
                               "label": label, "n_genes": int(n)})
    stability = pd.Series(np.round(same / combos, 3), index=table.index, name="specificity_stability")
    return stability, pd.DataFrame(counts)
