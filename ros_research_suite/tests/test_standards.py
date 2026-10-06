"""Tests der Standards-Anpassung (Phase 6b): M/A, Spezifität, Robustheit, Belege, Herkunftsangaben."""
import json

import numpy as np
import pandas as pd
import pytest

from src import differential_expression as de
from src import provenance as pv
from src import robustness as rb
from src import ros_evidence as ev
from src import statistics as st

TIMES = (30, 60, 100, 140, 180)


# ---------- M / A ----------
def make_expression_and_design():
    expression = pd.DataFrame({"s0": [5.0, 8.0], "s30": [7.0, 8.0], "s60": [4.0, 6.0]},
                              index=["geneA", "geneB"])
    design = pd.DataFrame({"treatment": ["H2O2"] * 3, "time_min": [0.0, 30.0, 60.0],
                           "is_baseline": [True, False, False]}, index=expression.columns)
    return expression, design


def test_ma_values_follow_the_definitions():
    expression, design = make_expression_and_design()
    ma = st.ma_values(expression, design, "H2O2", 30.0)
    assert ma.loc["geneA", "M"] == 2.0                  # 7 - 5
    assert ma.loc["geneA", "A"] == 6.0                  # (7 + 5) / 2
    assert ma.loc["geneB", "M"] == 0.0


def test_m_equals_log2_of_linear_ratio():
    expression, design = make_expression_and_design()
    ma = st.ma_values(expression, design, "H2O2", 60.0)
    linear_ratio = (2 ** expression["s60"]) / (2 ** expression["s0"])
    assert np.allclose(np.log2(linear_ratio), ma["M"])


def test_ma_and_log2fc_table_use_identical_definition():
    expression, design = make_expression_and_design()
    table = st.log2fc_vs_baseline(expression, design, "H2O2")
    assert np.allclose(table["log2fc_H2O2_30min"], st.ma_values(expression, design, "H2O2", 30.0)["M"])
    assert np.allclose(table["log2fc_H2O2_60min"], st.ma_values(expression, design, "H2O2", 60.0)["M"])


def test_ma_values_errors():
    expression, design = make_expression_and_design()
    with pytest.raises(st.StatisticsError):
        st.ma_values(expression, design, "H2O2", 999.0)
    with pytest.raises(st.StatisticsError):
        st.ma_values(expression, design.assign(is_baseline=False), "H2O2", 30.0)


def test_replication_status():
    _, design = make_expression_and_design()
    assert st.replication_status(design, "H2O2") == "not_testable_no_replicates"
    doubled = pd.concat([design, design.rename(index=lambda s: s + "b")])
    assert st.replication_status(doubled, "H2O2") == "replicates_present"


# ---------- Spezifität ----------
def make_two_treatment_table() -> pd.DataFrame:
    h2o2 = {  # Gen: (H2O2-Verlauf, MMS-Verlauf) bei 30, 60, 100 min
        "h2o2_specific": ([1.5, 1.5, 1.5], [0.1, 0.2, 0.1]),
        "mms_specific": ([0.1, 0.2, 0.1], [1.5, 1.5, 1.5]),
        "shared_up": ([1.5, 1.5, 1.5], [1.4, 1.4, 1.4]),
        "shared_down": ([-1.5, -1.5, -1.5], [-1.2, -1.2, -1.2]),
        "opposite": ([1.5, 1.5, 1.5], [-1.5, -1.5, -1.5]),
        "gray_zone": ([1.5, 1.5, 1.5], [0.7, 0.7, 0.0]),     # MMS 0,7: nicht responsiv, aber nicht gering
        "none": ([0.1, 0.1, 0.1], [0.2, 0.1, 0.1]),
    }
    table = pd.DataFrame(index=list(h2o2))
    for time_index, time in enumerate((30, 60, 100)):
        table[f"log2fc_H2O2_{time}min"] = [v[0][time_index] for v in h2o2.values()]
        table[f"log2fc_MMS_{time}min"] = [v[1][time_index] for v in h2o2.values()]
    return table


def test_specificity_with_gray_zone():
    table = make_two_treatment_table()
    result = de.classify_all(table, ["H2O2", "MMS"], 1.0, 2)
    labels = dict(result["specificity_H2O2_vs_MMS"])
    assert labels == {"h2o2_specific": "H2O2_specific_up", "mms_specific": "MMS_specific_up",
                      "shared_up": "shared_up", "shared_down": "shared_down",
                      "opposite": "opposite_or_mixed", "gray_zone": "ambiguous", "none": "none"}


def test_specificity_grayzone_width_is_configurable():
    table = make_two_treatment_table()
    result = de.classify_all(table, ["H2O2", "MMS"], 1.0, 2)
    wide = de.classify_specificity(table, result, "H2O2", "MMS", 1.0, 2, low_fraction=0.8)
    assert wide["gray_zone"] == "H2O2_specific_up"    # 0,7 < 0,8 gilt jetzt als "gering"


# ---------- Robustheit ----------
def make_robustness_table() -> pd.DataFrame:
    genes = {"strong": [3, 3, 3, 3, 3], "borderline": [0, 1.2, 1.2, 0, 0], "noise": [0.1, 0.2, 0.1, 0.0, 0.1],
             "onset60": [0, 1.2, 1.5, 1.4, 1.3]}
    table = pd.DataFrame(index=list(genes))
    for index, time in enumerate(TIMES):
        table[f"log2fc_H2O2_{time}min"] = [v[index] for v in genes.values()]
    return table


def test_responsive_robustness_scores_and_classes():
    scores = rb.responsive_robustness(make_robustness_table(), ["H2O2"])
    assert scores.loc["strong", "H2O2_robustness_score"] == 1.0
    assert scores.loc["strong", "H2O2_robustness_class"] == "robust"
    assert scores.loc["borderline", "H2O2_robustness_score"] == pytest.approx(4 / 12, abs=1e-3)
    assert scores.loc["borderline", "H2O2_robustness_class"] == "fragile"
    assert scores.loc["noise", "H2O2_robustness_class"] == "not_responsive"


def test_pattern_stability_and_share_overview():
    stability, overview = rb.pattern_stability(make_robustness_table(), ["H2O2"])
    assert stability.loc["strong", "H2O2_pattern_stability"] == 1.0
    # onset60: bei early_max=30 'late', bei 60 und 100 'early' -> 6 von 9 Kombinationen gleich
    assert stability.loc["onset60", "H2O2_pattern_stability"] == pytest.approx(6 / 9, abs=1e-3)
    assert np.isnan(stability.loc["noise", "H2O2_pattern_stability"])
    assert len(overview) == 9 and overview["share_transient"].between(0, 1).all()


def test_specificity_sensitivity_returns_stability_in_unit_range():
    table = make_two_treatment_table()
    stability, counts = rb.specificity_sensitivity(table, "H2O2", "MMS")
    assert stability.between(0, 1).all()
    assert stability["none"] == 1.0
    assert {"threshold", "low_fraction", "label", "n_genes"} <= set(counts.columns)


# ---------- Belege ----------
def test_go_term_parsing():
    terms = ev._go_terms("0006979 // response to oxidative stress // inferred from mutant phenotype "
                         "/// 0055114 // oxidation reduction // inferred from electronic annotation")
    assert terms[0] == ("response to oxidative stress", "inferred from mutant phenotype")
    assert len(terms) == 2 and ev._go_terms("") == []


def test_go_support_classification():
    from src.ros_genes import RosGene
    annotation = pd.DataFrame({
        "ORF": ["YGR088W", "YJR104C"],
        "Gene Ontology Biological Process": ["0006979 // response to oxidative stress // IEA",
                                             "0008150 // biological process unknown // ND"],
        "Gene Ontology Molecular Function": ["", ""]}, index=["p1", "p2"])
    genes = (RosGene("CTT1", "YGR088W", "antioxidant_defense", "catalase"),
             RosGene("SOD1", "YJR104C", "antioxidant_defense", "sod"),
             RosGene("GONE", "YAL001C", "antioxidant_defense", "x"))
    support = ev.go_support(annotation, genes).set_index("symbol")
    assert support.loc["CTT1", "go_support"] == "supported"
    assert support.loc["SOD1", "go_support"] == "no_matching_term"
    assert support.loc["GONE", "go_support"] == "not_in_annotation"


def test_positive_controls_compare_expected_and_observed():
    ros_table = pd.DataFrame({
        "symbol": ["A", "B"], "is_best_probe": [True, True],
        "H2O2_direction": ["up", "down"], "H2O2_peak_log2fc": [3.0, -2.0],
        "H2O2_responsive": [True, True]})
    controls = {"A": ("up", "Beleg 1"), "B": ("up", ""), "C": ("up", "")}
    result = ev.validate_positive_controls(ros_table, "H2O2", controls).set_index("symbol")
    assert result.loc["A", "status"] == "agrees"
    assert result.loc["B", "status"] == "does_not_agree" and result.loc["B", "source"] == "Quelle fehlt"
    assert result.loc["C", "status"] == "not_on_chip"


# ---------- Herkunftsangaben ----------
def test_worked_example_self_check():
    expression = pd.DataFrame({"s0": [5.0], "s30": [8.0]}, index=["p1"])
    design = pd.DataFrame({"treatment": ["H2O2", "H2O2"], "time_min": [0.0, 30.0],
                           "is_baseline": [True, False]}, index=expression.columns)
    example = pv.worked_example(expression, design, "H2O2", "p1", 30.0, "GENX")
    assert example["M"] == 3.0 and example["A"] == 6.5
    assert example["linear_ratio"] == pytest.approx(8.0)
    assert example["check_m"] and example["check_a"]
    assert any("M = " in line for line in example["lines"])


def test_update_metadata_accumulates_runs(tmp_path):
    info = {"geo_accession": "GSE1"}
    path = tmp_path / "meta.json"
    pv.update_metadata(path, info, "differential", {"probes_analyzed": 10}, ["results/a.csv"])
    pv.update_metadata(path, info, "timecourse", {"control_probes_removed": 2})
    data = json.loads(path.read_text(encoding="utf-8"))
    assert [run["analysis"] for run in data["analyses_run"]] == ["differential", "timecourse"]
    assert data["probe_counts"] == {"probes_analyzed": 10, "control_probes_removed": 2}
    assert data["parameters"]["MIN_ABS_LOG2FC"] == 1.0
    assert isinstance(data["parameters"]["SENSITIVITY_THRESHOLDS"], list)
    assert "python" in data["software_versions"] and "nicht testbar" in data["statistical_evidence"]


def test_build_dataset_info(tmp_path):
    samples = pd.DataFrame({"organism_ch1": ["Saccharomyces cerevisiae"] * 2, "platform_id": ["GPL1"] * 2,
                            "data_processing": ["RMA"] * 2}, index=["s0", "s30"])
    design = pd.DataFrame({"treatment": ["H2O2", "H2O2"], "time_min": [0.0, 30.0],
                           "concentration": ["0.3 mM", "0.3 mM"], "strain": ["Y262", "Y262"]},
                          index=["s0", "s30"])
    annotation_file = tmp_path / "GPL1.soft"
    annotation_file.write_text("x")
    annotation = pd.DataFrame({"Annotation Date": ["Mar 13, 2009", ""]})
    info = pv.build_dataset_info("GSE1", {"Series_title": "T"}, samples, design,
                                 tmp_path / "GSE1_series_matrix.txt.gz", annotation_file, annotation)
    assert info["treatments"]["H2O2"]["replicates_per_timepoint"] == 1
    assert info["pubmed_id_from_geo"] == "nicht in GEO-Metadaten angegeben"
    assert info["annotation"]["annotation_date_in_file"] == ["Mar 13, 2009"]
    assert "processed" in info["data_level"] and info["platform"] == ["GPL1"]
