"""Tests für die ROS-Genliste, die Zuordnung zu Sonden und die Auswertung."""
import re
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from src import ros_analysis as ra
from src import ros_genes as rg

TEST_GENES = (
    rg.RosGene("CTT1", "YGR088W", "antioxidant_defense", "catalase"),
    rg.RosGene("SOD1", "YJR104C", "antioxidant_defense", "superoxide dismutase"),
    rg.RosGene("RAD51", "YER095W", "dna_damage_response", "recombination"),
    rg.RosGene("FOO1", "YAL001C", "stress_response", "test"),      # nicht auf dem Chip
    rg.RosGene("XYZ1", "YBR001C", "stress_response", "test"),      # Symbol in Annotation anders
)


def make_result() -> pd.DataFrame:
    """Kleines Differential-Expression-Ergebnis (nur H2O2)."""
    rows = {
        # probe:      name,   ORF,       peak, up,   down,  responsive
        "p_ctt1":    ("CTT1", "YGR088W", 4.0, True, False, True),
        "p_sod1_a":  ("SOD1", "YJR104C", 0.4, False, False, False),   # schwache Sonde
        "p_sod1_b":  ("SOD1", "YJR104C", 1.5, True, False, True),     # starke Sonde = beste
        "p_rad51":   (None,   "YER095W", -0.2, False, False, False),  # Annotation ohne Symbol
        "p_xyz":     ("OTHER", "YBR001C", 2.0, True, False, True),
        "p_unknown": ("NEWG", "YZZ999W", 6.0, True, False, True),     # starker Antworter, nicht in Liste
        "p_quiet":   ("QUIET", "YQQ001W", 0.1, False, False, False),
    }
    result = pd.DataFrame(
        [{"display_name": name or orf, "gene_symbol": name, "systematic_name": orf,
          "H2O2_peak_log2fc": peak, "H2O2_peak_time_min": 60.0, "H2O2_n_timepoints_above": 3,
          "H2O2_up": up, "H2O2_down": down, "H2O2_responsive": resp,
          "H2O2_direction": "up" if up else "none"}
         for name, orf, peak, up, down, resp in rows.values()], index=list(rows))
    return result


def test_curated_list_is_consistent():
    orfs = [g.orf for g in rg.ROS_GENES]
    symbols = [g.symbol for g in rg.ROS_GENES]
    assert len(orfs) == len(set(orfs)) and len(symbols) == len(set(symbols))
    assert all(re.match(r"^Y[A-P][LR]\d{3}[CW](-[A-Z])?$", orf) for orf in orfs)   # ORF-Format
    assert set(Counter(g.category for g in rg.ROS_GENES)) == set(rg.CATEGORY_LABELS)


def test_known_examples_in_list():
    by_symbol = {g.symbol: g for g in rg.ROS_GENES}
    assert by_symbol["CTT1"].orf == "YGR088W"
    assert by_symbol["TSA2"].orf == "YDR453C"
    assert by_symbol["SRX1"].category == "antioxidant_defense"
    assert by_symbol["SOD1"].group == "superoxide dismutase"


def test_build_ros_table_matching_and_checks():
    table, missing = rg.build_ros_table(make_result(), genes=TEST_GENES)
    assert [g.symbol for g in missing] == ["FOO1"]
    assert table.loc["p_ctt1", "symbol_check"] == "ok"
    assert table.loc["p_rad51", "symbol_check"] == "annotation_has_no_symbol"
    assert table.loc["p_xyz", "symbol_check"] == "symbol_mismatch"
    assert table.loc["p_sod1_a", "n_probes"] == 2


def test_best_probe_is_the_strongest():
    table, _ = rg.build_ros_table(make_result(), genes=TEST_GENES)
    assert table.loc["p_sod1_b", "is_best_probe"] and not table.loc["p_sod1_a", "is_best_probe"]


def test_multi_gene_annotation_still_matches():
    result = make_result()
    result.loc["p_ctt1", "gene_symbol"] = "CTT1 /// OTHER"
    table, _ = rg.build_ros_table(result, genes=TEST_GENES)
    assert table.loc["p_ctt1", "symbol_check"] == "ok"


def test_missing_annotation_raises():
    with pytest.raises(rg.RosGeneError):
        rg.build_ros_table(make_result().drop(columns=["systematic_name"]), genes=TEST_GENES)


def test_nothing_found_raises():
    with pytest.raises(rg.RosGeneError):
        rg.build_ros_table(make_result(), genes=(rg.RosGene("A", "YAL001C", "stress_response", "x"),))


def test_category_summary_counts_best_probes_only():
    table, _ = rg.build_ros_table(make_result(), genes=TEST_GENES)
    summary = ra.category_summary(table, ["H2O2"])
    antioxidant = summary.loc[rg.CATEGORY_LABELS["antioxidant_defense"]]
    assert antioxidant["n_genes"] == 2            # SOD1 zählt nur einmal, obwohl 2 Sonden
    assert antioxidant["H2O2_up"] == 2            # CTT1 und beste SOD1-Sonde
    assert summary.loc[rg.CATEGORY_LABELS["dna_damage_response"], "H2O2_up"] == 0


def test_responsive_share():
    table, _ = rg.build_ros_table(make_result(), genes=TEST_GENES)
    in_list, rest = ra.responsive_share(make_result(), table, "H2O2")
    # In der Liste (beste Sonden): CTT1 (ja), SOD1_b (ja), RAD51 (nein), XYZ (ja) -> 3 von 4
    assert in_list == pytest.approx(0.75)
    assert rest == pytest.approx(1 / 3)          # Rest: p_sod1_a (nein), p_unknown (ja), p_quiet (nein)


def test_top_responsive_and_discovery():
    result = make_result()
    orfs = rg.known_orfs(TEST_GENES)
    top = ra.top_responsive(result, "H2O2", 3, orfs=orfs)
    assert list(top["systematic_name"]) == ["YZZ999W", "YGR088W", "YBR001C"]   # nach |peak|
    assert list(top["in_ros_list"]) == [False, True, True]
    assert list(top["rank"]) == [1, 2, 3]
    info = ra.discovery_summary(top, result, orfs)
    assert info["n_in_list"] == 2 and info["share_in_top"] == pytest.approx(2 / 3)
    assert info["base_rate"] == pytest.approx(5 / 7)   # 5 von 7 Sonden haben eine ORF aus TEST_GENES


def test_top_responsive_excludes_nonresponsive():
    top = ra.top_responsive(make_result(), "H2O2", 50, orfs=set())
    assert "YQQ001W" not in set(top["systematic_name"])


def test_ros_barplot_creates_file(tmp_path):
    pytest.importorskip("matplotlib")
    from src import visualization
    table, _ = rg.build_ros_table(make_result(), genes=TEST_GENES)
    path = visualization.ros_gene_barplot(table, ["H2O2"], rg.CATEGORY_LABELS, 1.0, tmp_path / "ros.png")
    assert path.exists() and path.stat().st_size > 0


def test_alias_in_annotation_counts_as_match():
    result = make_result()
    result.loc["p_ctt1", "gene_symbol"] = "OLDNAME"
    genes = (rg.RosGene("CTT1", "YGR088W", "antioxidant_defense", "catalase", aliases=("OLDNAME",)),)
    table, _ = rg.build_ros_table(result, genes=genes)
    assert table.loc["p_ctt1", "symbol_check"] == "ok"
    no_alias = (rg.RosGene("CTT1", "YGR088W", "antioxidant_defense", "catalase"),)
    table, _ = rg.build_ros_table(result, genes=no_alias)
    assert table.loc["p_ctt1", "symbol_check"] == "symbol_mismatch"


def test_gpx3_has_hyr1_alias():
    gpx3 = next(g for g in rg.ROS_GENES if g.symbol == "GPX3")
    assert "HYR1" in gpx3.aliases and gpx3.orf == "YIR037W"
