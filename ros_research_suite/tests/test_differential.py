"""Tests für die Differential Expression ohne Replikate (Werte von Hand gewählt)."""
import numpy as np
import pandas as pd
import pytest

from src import differential_expression as de


def make_table() -> pd.DataFrame:
    # Spaltenreihenfolge absichtlich NICHT nach Zeit sortiert (100min vor 60min)
    data = {
        "log2fc_H2O2_30min":  [0.0, 1.2, -1.5, 1.2, 0.2],
        "log2fc_H2O2_100min": [1.5, 1.3, 0.0, -1.2, 0.0],
        "log2fc_H2O2_60min":  [1.2, 0.0, -1.1, 1.2, 0.1],
        "log2fc_H2O2_140min": [1.0, 0.0, 0.0, -1.2, 0.3],
    }
    return pd.DataFrame(data, index=["A_up", "B_gap", "C_down", "D_both", "E_none"])


def test_timepoint_columns_sorted_numerically():
    columns = de.timepoint_columns(make_table(), "H2O2")
    assert [time for time, _ in columns] == [30.0, 60.0, 100.0, 140.0]
    assert de.timepoint_columns(make_table(), "MMS") == []


def test_longest_true_run():
    mask = np.array([[False, True, True, False, True],
                     [True, True, True, True, True],
                     [False, False, False, False, False]])
    assert list(de.longest_true_run(mask)) == [2, 5, 0]


def test_classification_with_consecutive_rule():
    result = de.classify_treatment(make_table(), "H2O2", min_abs_log2fc=1.0, min_consecutive=2)
    # Zeitreihen (30, 60, 100, 140 min):
    # A_up   0.0, 1.2, 1.5, 1.0  -> 3 in Folge hoch
    # B_gap  1.2, 0.0, 1.3, 0.0  -> zwei hohe Werte, aber nie nacheinander
    # C_down -1.5,-1.1, 0.0, 0.0 -> 2 in Folge runter
    # D_both 1.2, 1.2,-1.2,-1.2  -> erst hoch, dann runter
    # E_none klein
    assert dict(result["H2O2_direction"]) == {
        "A_up": "up", "B_gap": "none", "C_down": "down", "D_both": "both", "E_none": "none"}
    assert result.loc["B_gap", "H2O2_n_timepoints_above"] == 2   # Schwelle 2x erreicht, aber nicht in Folge
    assert not result.loc["B_gap", "H2O2_responsive"]


def test_min_consecutive_one_accepts_single_timepoints():
    result = de.classify_treatment(make_table(), "H2O2", min_abs_log2fc=1.0, min_consecutive=1)
    assert result.loc["B_gap", "H2O2_direction"] == "up"


def test_threshold_is_respected():
    result = de.classify_treatment(make_table(), "H2O2", min_abs_log2fc=1.4, min_consecutive=1)
    assert result.loc["A_up", "H2O2_direction"] == "up"       # 1.5 >= 1.4
    assert result.loc["D_both", "H2O2_direction"] == "none"   # 1.2 < 1.4


def test_peak_values_and_time():
    result = de.classify_treatment(make_table(), "H2O2", 1.0, 2)
    assert result.loc["A_up", "H2O2_peak_log2fc"] == 1.5
    assert result.loc["A_up", "H2O2_peak_time_min"] == 100.0
    assert result.loc["C_down", "H2O2_peak_log2fc"] == -1.5
    assert result.loc["C_down", "H2O2_peak_time_min"] == 30.0


def test_unknown_treatment_raises():
    with pytest.raises(de.DifferentialExpressionError):
        de.classify_treatment(make_table(), "MMS", 1.0, 2)


def test_split_control_probes():
    table = pd.DataFrame({"x": [1, 2, 3]}, index=["AFFX-r2-TagA_at", "1769308_at", "AFFX-CreX-3_at"])
    real, controls = de.split_control_probes(table)
    assert list(real.index) == ["1769308_at"] and len(controls) == 2


def _two_treatment_result():
    table = make_table()
    for column in list(table.columns):
        table[column.replace("H2O2", "MMS")] = table[column]
    # MMS: A bleibt up, B wird 'up' (Lücke aufheben), C bleibt, E bleibt
    table.loc["B_gap", [c for c in table.columns if "MMS" in c]] = [1.2, 1.2, 1.2, 1.2]
    table.loc["C_down", [c for c in table.columns if "MMS" in c]] = [0, 0, 0, 0]
    table.loc["D_both", [c for c in table.columns if "MMS" in c]] = [-1.5, -1.5, -1.5, -1.5]
    return de.classify_all(table, ["H2O2", "MMS"], 1.0, 2)


def test_compare_treatments_categories():
    comparison = _two_treatment_result()["comparison_H2O2_vs_MMS"]
    assert comparison["A_up"] == "both_up"
    assert comparison["B_gap"] == "MMS_only_up"
    assert comparison["C_down"] == "H2O2_only_down"
    assert comparison["D_both"] == "different_direction"   # H2O2 'both' vs MMS 'down'
    assert comparison["E_none"] == "none"


def test_top_genes_ordering():
    result = _two_treatment_result()
    up = de.top_genes(result, "H2O2", "up", 2)
    assert list(up.index) == ["A_up", "D_both"]          # max_log2fc: 1.5, 1.2
    down = de.top_genes(result, "H2O2", "down", 5)
    assert down.index[0] == "C_down"                     # stärkster Abfall (-1.5)
    with pytest.raises(ValueError):
        de.top_genes(result, "H2O2", "sideways", 1)


def test_ma_and_response_plots_create_files(tmp_path):
    pytest.importorskip("matplotlib")
    from src import visualization
    result = _two_treatment_result()
    expression = pd.DataFrame(np.random.default_rng(1).uniform(4, 10, (5, 5)),
                              index=result.index, columns=["b", "t30", "t60", "t100", "t140"])
    design = pd.DataFrame({"treatment": ["H2O2"] * 5, "time_min": [0.0, 30, 60, 100, 140],
                           "is_baseline": [True, False, False, False, False]}, index=expression.columns)
    ma = visualization.ma_plots(expression, design, result, "H2O2", 1.0, tmp_path / "ma.png")
    sc = visualization.response_scatter(result, "H2O2", 1.0, tmp_path / "sc.png")
    assert ma.exists() and ma.stat().st_size > 0 and sc.exists() and sc.stat().st_size > 0


def test_response_scatter_handles_many_identical_coordinates(tmp_path):
    pytest.importorskip("matplotlib")
    from src import visualization
    n = 40   # viele Gene mit (fast) gleichen Koordinaten: früher überlagerten sich die Namen
    rng = np.random.default_rng(2)
    result = pd.DataFrame({
        "display_name": [f"GENE{i}" for i in range(n)],
        "H2O2_peak_log2fc": np.concatenate([rng.uniform(2, 3, 20), rng.uniform(-3, -2, 20)]),
        "H2O2_n_timepoints_above": [5] * n,
        "H2O2_responsive": [True] * n}, index=[f"p{i}" for i in range(n)])
    path = visualization.response_scatter(result, "H2O2", 1.0, tmp_path / "scatter.png")
    assert path.exists() and path.stat().st_size > 0
