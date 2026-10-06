"""Statistische Grundfunktionen der ROS Research Suite.

Wichtig für GSE12220: Die Expressionswerte sind bereits log2-transformiert (RMA).
Dann gilt   log2FC = Wert(behandelt) - Wert(Baseline)   (Differenz, kein Quotient!).
Der Datensatz hat KEINE Replikate, deshalb werden hier zwar p-Wert-Funktionen
bereitgestellt, aber für diesen Datensatz nicht verwendet.
"""
import numpy as np
import pandas as pd


class StatisticsError(Exception):
    """Berechnung nicht möglich (z. B. keine Baseline vorhanden)."""


# --- Beschreibende Statistik --------------------------------------------------
def _clean(values) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array[~np.isnan(array)]


def mean(values) -> float:
    """Arithmetischer Mittelwert (NaN werden ignoriert)."""
    return float(np.mean(_clean(values)))


def median(values) -> float:
    """Median: der mittlere Wert, robust gegen Ausreißer."""
    return float(np.median(_clean(values)))


def standard_deviation(values) -> float:
    """Stichproben-Standardabweichung (n-1 im Nenner). Bei n<2 nicht definiert -> NaN."""
    array = _clean(values)
    return float(np.std(array, ddof=1)) if array.size >= 2 else float("nan")


def standard_error(values) -> float:
    """Standardfehler des Mittelwerts = Standardabweichung / Wurzel(n)."""
    array = _clean(values)
    return standard_deviation(array) / float(np.sqrt(array.size)) if array.size else float("nan")


# --- Fold Change ----------------------------------------------------------------
def fold_change(treated, baseline):
    """Fold Change = behandelt / Baseline, für Werte auf LINEARER Skala (z. B. Rohintensitäten)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(treated, baseline)


def log2_fold_change(treated, baseline):
    """log2 des Fold Change, für Werte auf LINEARER Skala."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(np.divide(treated, baseline))


def log2_fold_change_from_log2(treated_log2, baseline_log2):
    """log2 Fold Change für Werte, die schon log2-transformiert sind: einfach die Differenz.

    Hintergrund: log2(a/b) = log2(a) - log2(b).
    """
    return np.subtract(treated_log2, baseline_log2)


# --- Tests und Mehrfachtest-Korrektur --------------------------------------------
def welch_t_test_p_values(group_a: pd.DataFrame, group_b: pd.DataFrame) -> pd.Series:
    """Welch-t-Test pro Gen (Zeilen = Gene, Spalten = Replikate).

    Biologisch: Unterscheiden sich die Mittelwerte zweier Gruppen mehr, als es
    das Rauschen zwischen Replikaten erklärt?
    Mathematisch: t = Mittelwertdifferenz / Standardfehler der Differenz, ohne
    Annahme gleicher Varianzen.
    Limitation: Braucht mind. 2 Replikate pro Gruppe; sonst kommt NaN zurück
    (kein erfundener p-Wert).
    """
    if group_a.shape[1] < 2 or group_b.shape[1] < 2:
        return pd.Series(np.nan, index=group_a.index, name="p_value")
    from scipy import stats  # erst hier importiert: nur für Datensätze mit Replikaten nötig
    result = stats.ttest_ind(group_a.to_numpy(dtype=float), group_b.to_numpy(dtype=float),
                             axis=1, equal_var=False, nan_policy="omit")
    return pd.Series(result.pvalue, index=group_a.index, name="p_value")


def benjamini_hochberg(p_values) -> np.ndarray:
    """Benjamini-Hochberg-Korrektur (kontrolliert die False Discovery Rate, FDR).

    Biologisch: Bei ~6000 Genen wären bei p<0,05 allein durch Zufall ~300 "Treffer".
    Mathematisch: p sortieren (Rang i von m), p_adj = p * m / i, danach von hinten
    das laufende Minimum bilden (damit die Reihenfolge erhalten bleibt), max. 1.
    Limitation: Gilt für unabhängige/positiv korrelierte Tests; NaN bleibt NaN.
    """
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan)
    valid = ~np.isnan(p)
    valid_p = p[valid]
    m = valid_p.size
    if m == 0:
        return adjusted
    order = np.argsort(valid_p)
    scaled = valid_p[order] * m / np.arange(1, m + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    result = np.empty(m)
    result[order] = np.minimum(scaled, 1.0)
    adjusted[valid] = result
    return adjusted


# --- Log2FC gegen Baseline (Zeitreihe) -------------------------------------------
def log2fc_vs_baseline(expression: pd.DataFrame, design: pd.DataFrame,
                       treatment: str) -> pd.DataFrame:
    """log2FC jedes Zeitpunkts gegen die Baseline (t=0) einer Behandlung.

    Mehrere Samples pro Zeitpunkt/Baseline werden gemittelt (log-Skala).
    Spalten: 'log2fc_<Behandlung>_<Zeit>min'.
    """
    rows = design[design["treatment"] == treatment]
    baseline_mask = rows["is_baseline"].to_numpy(dtype=bool)
    baseline_ids = list(rows.index[baseline_mask])
    if not baseline_ids:
        raise StatisticsError(f"Für '{treatment}' gibt es kein Baseline-Sample (t=0).")
    baseline = expression[baseline_ids].mean(axis=1)

    later = rows[~baseline_mask]
    columns = {}
    for time_min in sorted(later["time_min"].unique()):
        ids = list(later.index[later["time_min"] == time_min])
        columns[f"log2fc_{treatment}_{time_min:g}min"] = (
            log2_fold_change_from_log2(expression[ids].mean(axis=1), baseline))
    return pd.DataFrame(columns, index=expression.index)


def log2fc_table(expression: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """log2FC-Tabelle für alle Behandlungen + je Behandlung max. |log2FC| über die Zeit."""
    parts = []
    for treatment in design["treatment"].unique():
        part = log2fc_vs_baseline(expression, design, treatment)
        if part.empty:
            continue
        part[f"max_abs_log2fc_{treatment}"] = part.abs().max(axis=1)
        parts.append(part)
    if not parts:
        raise StatisticsError("Keine Behandlung mit Zeitpunkten nach der Baseline gefunden.")
    return pd.concat(parts, axis=1)


def ma_values(expression: pd.DataFrame, design: pd.DataFrame, treatment: str,
              time_min: float) -> pd.DataFrame:
    """M und A eines Zeitpunkts gegen die Baseline (t=0) derselben Behandlung.

    Eingabe sind log2-skalierte Werte (RMA). Definitionen:
        M = E_behandelt - E_Baseline          (= log2 des Verhältnisses = log2FC)
        A = (E_behandelt + E_Baseline) / 2     (mittlere log2-Expression)
    Das ist die EINZIGE Stelle im Projekt, an der M und A berechnet werden (MA-Plots,
    Beispielrechnung und Tests nutzen sie). Mehrere Samples pro Zeitpunkt werden gemittelt.
    """
    rows = design[design["treatment"] == treatment]
    is_baseline = rows["is_baseline"].to_numpy(dtype=bool)
    baseline_ids = list(rows.index[is_baseline])
    if not baseline_ids:
        raise StatisticsError(f"Für '{treatment}' gibt es kein Baseline-Sample (t=0).")
    at_time = (rows["time_min"] == time_min).to_numpy()
    treated_ids = list(rows.index[at_time & ~is_baseline])
    if not treated_ids:
        raise StatisticsError(f"Für '{treatment}' gibt es keinen Zeitpunkt {time_min:g} min.")
    baseline = expression[baseline_ids].mean(axis=1)
    treated = expression[treated_ids].mean(axis=1)
    return pd.DataFrame({"M": treated - baseline, "A": (treated + baseline) / 2})


def replication_status(design: pd.DataFrame, treatment: str) -> str:
    """Gibt an, ob für eine Behandlung statistische Tests überhaupt möglich sind.

    'not_testable_no_replicates': höchstens 1 Sample je Zeitpunkt -> keine Varianzschätzung,
    keine p-Werte. 'replicates_present': Tests wären möglich (noch nicht ausgeführt).
    """
    rows = design[design["treatment"] == treatment]
    most_replicates = rows.groupby("time_min").size().max()
    return "not_testable_no_replicates" if most_replicates < 2 else "replicates_present"
