"""Tests für die statistischen Grundfunktionen (alle Werte von Hand nachgerechnet)."""
import math

import numpy as np
import pandas as pd
import pytest

from src import statistics as st


def test_descriptive_statistics():
    values = [1, 2, 3, 4]
    assert st.mean(values) == 2.5
    assert st.median(values) == 2.5
    assert st.standard_deviation(values) == pytest.approx(1.2909944)   # n-1 im Nenner
    assert st.standard_error(values) == pytest.approx(1.2909944 / 2)   # sd / sqrt(4)


def test_nan_is_ignored_and_single_value_has_no_sd():
    assert st.mean([1, np.nan, 3]) == 2.0
    assert math.isnan(st.standard_deviation([5.0]))


def test_fold_change_and_log2_fold_change():
    assert st.fold_change(8, 2) == 4
    assert st.log2_fold_change(8, 2) == 2
    assert st.log2_fold_change(1, 4) == -2          # Herunterregulation = negativ
    assert st.log2_fold_change_from_log2(5.0, 3.0) == 2.0


def test_fold_change_division_by_zero_is_inf_not_crash():
    assert math.isinf(st.fold_change(1.0, 0.0))


def test_benjamini_hochberg_known_example():
    adjusted = st.benjamini_hochberg([0.01, 0.04, 0.03, 0.005])
    assert adjusted == pytest.approx([0.02, 0.04, 0.04, 0.02])


def test_benjamini_hochberg_keeps_nan_and_caps_at_one():
    adjusted = st.benjamini_hochberg([0.01, np.nan, 0.04])
    assert adjusted[0] == pytest.approx(0.02) and np.isnan(adjusted[1])
    assert adjusted[2] == pytest.approx(0.04)
    assert st.benjamini_hochberg([0.9])[0] == 0.9
    assert (st.benjamini_hochberg([0.5, 0.6, 0.99]) <= 1).all()


def test_welch_t_test_matches_reference():
    pytest.importorskip("scipy")
    a = pd.DataFrame([[1, 2, 3]], index=["g1"])
    b = pd.DataFrame([[4, 5, 6]], index=["g1"])
    assert st.welch_t_test_p_values(a, b)["g1"] == pytest.approx(0.02131, abs=1e-4)


def test_welch_t_test_without_replicates_gives_nan_not_invented_p():
    a = pd.DataFrame([[1.0]], index=["g1"])
    b = pd.DataFrame([[4.0]], index=["g1"])
    assert st.welch_t_test_p_values(a, b).isna().all()


def _example():
    expression = pd.DataFrame(
        {"s0": [5.0, 8.0], "s30": [7.0, 8.0], "s60": [4.0, 6.0], "m0": [1.0, 1.0], "m30": [1.0, 3.0]},
        index=["geneA", "geneB"])
    design = pd.DataFrame({
        "treatment": ["H2O2", "H2O2", "H2O2", "MMS", "MMS"],
        "time_min": [0.0, 30.0, 60.0, 0.0, 30.0],
        "is_baseline": [True, False, False, True, False]}, index=expression.columns)
    return expression, design


def test_log2fc_vs_baseline():
    expression, design = _example()
    result = st.log2fc_vs_baseline(expression, design, "H2O2")
    assert list(result.columns) == ["log2fc_H2O2_30min", "log2fc_H2O2_60min"]
    assert result.loc["geneA", "log2fc_H2O2_30min"] == 2.0     # 7 - 5
    assert result.loc["geneB", "log2fc_H2O2_60min"] == -2.0    # 6 - 8


def test_log2fc_averages_replicate_baselines():
    expression, design = _example()
    expression["s0b"] = [7.0, 10.0]
    design.loc["s0b"] = ["H2O2", 0.0, True]
    result = st.log2fc_vs_baseline(expression, design, "H2O2")
    assert result.loc["geneA", "log2fc_H2O2_30min"] == 1.0     # 7 - mean(5, 7)


def test_missing_baseline_raises():
    expression, design = _example()
    design["is_baseline"] = False
    with pytest.raises(st.StatisticsError):
        st.log2fc_vs_baseline(expression, design, "H2O2")


def test_log2fc_table_has_both_treatments_and_max_abs():
    expression, design = _example()
    table = st.log2fc_table(expression, design)
    assert "log2fc_MMS_30min" in table.columns
    assert table.loc["geneB", "max_abs_log2fc_H2O2"] == 2.0
    assert table.loc["geneB", "log2fc_MMS_30min"] == 2.0
