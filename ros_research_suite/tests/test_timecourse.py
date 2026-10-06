"""Tests für die Zeitverlaufs-Muster (Zeitpunkte 30, 60, 100, 140, 180 min; Schwelle 1,0; 2 in Folge)."""
import numpy as np
import pandas as pd
import pytest

from src import differential_expression as de
from src import timecourse as tc

TIMES = [30, 60, 100, 140, 180]


def make_table() -> pd.DataFrame:
    genes = {
        "early_sustained_up":    [1.5, 2.0, 2.2, 2.0, 1.8],
        "early_transient_up":    [2.0, 3.0, 1.5, 0.8, 0.4],    # Ende 0,4 <= 0,5 * 3,0
        "early_partial_up":      [1.2, 1.6, 1.3, 1.1, 0.9],    # Ende 0,9: unter Schwelle, aber > 0,8
        "late_sustained_down":   [0.0, -0.2, -1.2, -1.5, -1.6],
        "onset60_is_early":      [0.0, 1.2, 1.5, 1.4, 1.3],
        "onset100_is_late":      [0.0, 0.0, 1.2, 1.3, 1.4],
        "biphasic":              [1.5, 1.5, -1.5, -1.5, -1.5],
        "noise":                 [0.2, 0.3, 0.1, 0.5, 0.4],
        "gap_not_consecutive":   [1.5, 0.0, 1.5, 0.0, 0.0],
    }
    table = pd.DataFrame(genes, index=[f"log2fc_H2O2_{t}min" for t in TIMES]).T
    # Spalten absichtlich mischen: die Sortierung nach Zeit muss die Klasse selbst herstellen
    return table[[f"log2fc_H2O2_{t}min" for t in (100, 30, 180, 60, 140)]]


def test_patterns():
    result = tc.classify_patterns(make_table(), "H2O2", 1.0, 2, 60, 0.5)
    assert dict(result["H2O2_pattern"]) == {
        "early_sustained_up": "early_sustained_up",
        "early_transient_up": "early_transient_up",
        "early_partial_up": "early_partial_decline_up",
        "late_sustained_down": "late_sustained_down",
        "onset60_is_early": "early_sustained_up",
        "onset100_is_late": "late_sustained_up",
        "biphasic": "biphasic",
        "noise": "none",
        "gap_not_consecutive": "none",
    }


def test_onset_time_and_final_value():
    result = tc.classify_patterns(make_table(), "H2O2", 1.0, 2, 60, 0.5)
    assert result.loc["late_sustained_down", "H2O2_onset_time_min"] == 100
    assert result.loc["early_sustained_up", "H2O2_onset_time_min"] == 30
    assert result.loc["early_transient_up", "H2O2_final_log2fc"] == 0.4
    assert np.isnan(result.loc["noise", "H2O2_onset_time_min"])


def test_early_boundary_is_configurable():
    result = tc.classify_patterns(make_table(), "H2O2", 1.0, 2, early_max_min=30, return_fraction=0.5)
    assert result.loc["onset60_is_early", "H2O2_pattern"] == "late_sustained_up"   # Beginn 60 > 30


def test_return_fraction_is_configurable():
    strict = tc.classify_patterns(make_table(), "H2O2", 1.0, 2, 60, return_fraction=0.1)
    assert strict.loc["early_transient_up", "H2O2_pattern"] == "early_partial_decline_up"


def test_single_timepoint_rule_changes_result():
    result = tc.classify_patterns(make_table(), "H2O2", 1.0, 1, 60, 0.5)
    assert result.loc["gap_not_consecutive", "H2O2_pattern"] != "none"


def test_run_longer_than_series_gives_no_response():
    result = tc.classify_patterns(make_table(), "H2O2", 1.0, 9, 60, 0.5)
    assert (result["H2O2_pattern"] == "none").all()


def test_unknown_treatment_raises():
    with pytest.raises(de.DifferentialExpressionError):
        tc.classify_patterns(make_table(), "MMS")


def test_classify_patterns_all_keeps_gene_names():
    table = make_table()
    table.insert(0, "display_name", [f"G{i}" for i in range(len(table))])
    result = tc.classify_patterns_all(table, ["H2O2"])
    assert result.columns[0] == "display_name" and "H2O2_pattern" in result.columns


def test_threshold_sensitivity_decreases_with_stricter_criteria():
    table = make_table()
    sensitivity = de.threshold_sensitivity(table, ["H2O2"], thresholds=(0.5, 1.0, 2.0),
                                           consecutive_values=(1, 2))
    counts = sensitivity.set_index(["threshold", "min_consecutive"])["n_responsive"]
    assert counts[(0.5, 1)] >= counts[(1.0, 1)] >= counts[(2.0, 1)]
    assert counts[(1.0, 1)] >= counts[(1.0, 2)]
    assert counts[(1.0, 2)] == 7        # alle außer 'noise' und 'gap_not_consecutive'


def test_timecourse_plots_create_files(tmp_path):
    pytest.importorskip("matplotlib")
    from src import visualization
    table = make_table()
    for column in list(table.columns):
        table[column.replace("H2O2", "MMS")] = table[column] * 0.5
    patterns = tc.classify_patterns(table, "H2O2", 1.0, 2, 60, 0.5)["H2O2_pattern"]
    grid = visualization.gene_timecourse_grid(table, list(table.index[:6]), ["H2O2", "MMS"], 1.0,
                                              tmp_path / "grid.png")
    profiles = visualization.pattern_profiles(table, patterns, "H2O2", 1.0, tmp_path / "prof.png",
                                              min_genes=1)
    assert grid.stat().st_size > 0 and profiles.stat().st_size > 0
