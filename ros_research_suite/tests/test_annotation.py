"""Tests für die Plattform-Annotation (ohne Internet)."""
import pandas as pd
import pytest

import config
from src import annotation as ann
from src.geo_client import GeoClient, GeoInvalidIdError, normalize_platform_id, platform_query_url

PLATFORM_TEXT = (
    "^PLATFORM = GPL1\n!Platform_title = Test\n!platform_table_begin\n"
    "ID\tGB_ACC\tGene Symbol\tORF\n"
    "1_at\tX1\tCAT\tYGR088W\n"
    "2_at\tX2\t---\tYLR109W\n"
    "!platform_table_end\n^SAMPLE = GSM1\n"
)


@pytest.fixture
def platform_file(tmp_path):
    path = tmp_path / "GPL1.soft"
    path.write_text(PLATFORM_TEXT)
    return path


def test_platform_ids():
    assert normalize_platform_id(" gpl2529 ") == "GPL2529"
    with pytest.raises(GeoInvalidIdError):
        normalize_platform_id("GSE1")
    url = platform_query_url("GPL2529")
    assert "acc=GPL2529" in url and "view=full" in url


def test_parse_platform_table_stops_at_end_marker(platform_file):
    table = ann.parse_platform_table(platform_file)
    assert list(table.index) == ["1_at", "2_at"]
    assert table.loc["1_at", "Gene Symbol"] == "CAT"


def test_parse_without_table_gives_clear_error(tmp_path):
    bad = tmp_path / "bad.soft"
    bad.write_text("<html>Fehlerseite</html>\n")
    with pytest.raises(ann.AnnotationError):
        ann.parse_platform_table(bad)


def test_detect_gene_columns(platform_file):
    table = ann.parse_platform_table(platform_file)
    assert ann.detect_gene_columns(table) == {
        "gene_symbol": "Gene Symbol", "systematic_name": "ORF", "organism": None}


def test_detect_returns_none_if_nothing_matches():
    table = pd.DataFrame({"GB_ACC": ["x"]}, index=["1_at"])
    assert ann.detect_gene_columns(table) == {
        "gene_symbol": None, "systematic_name": None, "organism": None}


def test_config_override_must_exist(platform_file, monkeypatch):
    monkeypatch.setattr(config, "ANNOTATION_SYMBOL_COLUMN", "gibt_es_nicht")
    with pytest.raises(ann.AnnotationError):
        ann.detect_gene_columns(ann.parse_platform_table(platform_file))


def test_match_rate_and_attach(platform_file):
    table = ann.parse_platform_table(platform_file)
    probes = pd.Index(["1_at", "2_at", "3_at", "4_at"])
    assert ann.probe_match_rate(probes, table) == 0.5

    results = pd.DataFrame({"log2fc": [1.0, 2.0, 3.0]}, index=["1_at", "2_at", "3_at"])
    merged = ann.attach_annotation(results, table, ann.detect_gene_columns(table))
    assert merged.loc["1_at", "gene_symbol"] == "CAT"
    assert pd.isna(merged.loc["2_at", "gene_symbol"])      # '---' wird zu leer
    assert pd.isna(merged.loc["3_at", "systematic_name"])  # nicht in Annotation -> leer, nicht geraten
    assert list(merged["log2fc"]) == [1.0, 2.0, 3.0]


class _Response:
    status_code = 200

    def iter_content(self, size):
        yield PLATFORM_TEXT.encode()


class _Session:
    def __init__(self):
        self.calls = []

    def get(self, url, stream=False, timeout=None):
        self.calls.append(url)
        return _Response()


def test_fetch_platform_annotation_uses_cache(tmp_path):
    session = _Session()
    client = GeoClient(tmp_path, session=session)
    path = client.fetch_platform_annotation("GPL1")
    assert path == tmp_path / "GPL1" / "GPL1.soft" and path.exists()
    client.fetch_platform_annotation("GPL1")
    assert len(session.calls) == 1


def test_display_name_falls_back_to_systematic_name(platform_file):
    table = ann.parse_platform_table(platform_file)
    results = pd.DataFrame({"x": [1.0, 2.0]}, index=["1_at", "2_at"])
    merged = ann.attach_annotation(results, table, ann.detect_gene_columns(table))
    assert merged.loc["1_at", "display_name"] == "CAT"
    assert merged.loc["2_at", "display_name"] == "YLR109W"   # Symbol war '---'


def test_filter_by_organism_keeps_unknown_removes_other():
    table = pd.DataFrame({"organism": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe", None],
                          "v": [1, 2, 3]}, index=["a", "b", "c"])
    filtered, removed = ann.filter_by_organism(table, "Saccharomyces cerevisiae")
    assert list(filtered.index) == ["a", "c"] and removed == 1
    untouched, removed = ann.filter_by_organism(table.drop(columns="organism"), "x")
    assert removed == 0 and len(untouched) == 3
