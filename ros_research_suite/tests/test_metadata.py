"""Tests für Metadaten-Parser und Datenimport mit einer kleinen Mini-Datei."""
import gzip

import pandas as pd
import pytest

from src import metadata

MINI_MATRIX = (
    '!Series_title\t"Mini test series"\n'
    '!Sample_title\t"ctrl 0 min"\t"H2O2 30 min"\t"H2O2 2 h"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n'
    '!Sample_characteristics_ch1\t"strain: BY4741"\t"strain: BY4741"\t"strain: BY4741"\n'
    '!series_matrix_table_begin\n'
    '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\n'
    '"YGR088W"\t1.0\t3.5\tnull\n'
    '"YLR109W"\t2.0\t2.5\t2.0\n'
    '!series_matrix_table_end\n'
)


@pytest.fixture
def matrix_file(tmp_path):
    path = tmp_path / "GSE0_series_matrix.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(MINI_MATRIX)
    return path


def test_header_and_samples(matrix_file):
    info, rows = metadata.read_header(matrix_file)
    assert info["Series_title"] == "Mini test series"
    table = metadata.build_sample_table(rows)
    assert list(table.index) == ["GSM1", "GSM2", "GSM3"]
    assert table.loc["GSM2", "strain"] == "BY4741"


def test_expression_table(matrix_file):
    expr = metadata.read_expression_table(matrix_file)
    assert expr.shape == (2, 3)
    assert expr.loc["YGR088W", "GSM2"] == 3.5
    assert expr.isna().loc["YGR088W", "GSM3"]  # 'null' -> NaN


def test_align_expression_orders_columns_and_checks_samples(matrix_file):
    expr = metadata.read_expression_table(matrix_file)
    design = pd.DataFrame(index=["GSM3", "GSM1"])
    assert list(metadata.align_expression(expr, design).columns) == ["GSM3", "GSM1"]
    with pytest.raises(metadata.MetadataError):
        metadata.align_expression(expr, pd.DataFrame(index=["GSM1", "GSM99"]))


def test_describe_expression(matrix_file):
    info = metadata.describe_expression(metadata.read_expression_table(matrix_file))
    assert info["n_probes"] == 2 and info["n_samples"] == 3
    assert info["n_missing"] == 1
    assert info["max"] == 3.5 and info["likely_log_scale"]


def test_wrong_format_gives_clear_error(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("irgendein Text ohne Tabelle\n")
    with pytest.raises(metadata.MetadataError):
        metadata.read_header(bad)
