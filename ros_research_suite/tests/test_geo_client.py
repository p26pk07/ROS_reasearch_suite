"""Tests für den GEO-Client. Kein Internet nötig: das Netzwerk wird simuliert."""
import pytest
import requests

from src import geo_client
from src.geo_client import (GeoClient, GeoInvalidIdError, GeoNetworkError,
                            GeoNotFoundError, normalize_geo_id,
                            parse_matrix_filenames, series_directory_url)


class FakeResponse:
    def __init__(self, status=200, text="", content=b""):
        self.status_code, self.text, self._content = status, text, content

    def iter_content(self, chunk_size):
        yield self._content

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(str(self.status_code))


class FakeSession:
    def __init__(self, routes):
        self.routes, self.calls = routes, []

    def get(self, url, stream=False, timeout=None):
        self.calls.append(url)
        result = self.routes.get(url)
        if isinstance(result, Exception):
            raise result
        return result if result is not None else FakeResponse(404)


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(geo_client.time, "sleep", lambda s: None)


def test_normalize_geo_id():
    assert normalize_geo_id(" gse12220 ") == "GSE12220"


@pytest.mark.parametrize("bad", ["12220", "GSM123", "GSE", "GSE12a", ""])
def test_invalid_ids_rejected(bad):
    with pytest.raises(GeoInvalidIdError):
        normalize_geo_id(bad)


def test_series_directory_url():
    url = series_directory_url("GSE12220")
    assert url.endswith("/GSE12nnn/GSE12220")


def test_parse_matrix_filenames_multiple_platforms():
    html = ('<a href="GSE1-GPL1_series_matrix.txt.gz">x</a>'
            '<a href="GSE1-GPL2_series_matrix.txt.gz">y</a><a href="other.txt">z</a>')
    assert parse_matrix_filenames(html) == [
        "GSE1-GPL1_series_matrix.txt.gz", "GSE1-GPL2_series_matrix.txt.gz"]


def _client_with_dataset(tmp_path):
    base = series_directory_url("GSE1")
    session = FakeSession({
        base + "/matrix/": FakeResponse(text='<a href="GSE1_series_matrix.txt.gz">f</a>'),
        base + "/matrix/GSE1_series_matrix.txt.gz": FakeResponse(content=b"DATA"),
    })
    return GeoClient(tmp_path, session=session), session


def test_download_then_cache(tmp_path):
    client, session = _client_with_dataset(tmp_path)
    files = client.fetch_series_matrix("GSE1")
    assert files[0].read_bytes() == b"DATA"
    calls_after_first = len(session.calls)

    client.fetch_series_matrix("GSE1")  # zweiter Aufruf: darf nichts laden
    assert len(session.calls) == calls_after_first


def test_force_redownloads(tmp_path):
    client, session = _client_with_dataset(tmp_path)
    client.fetch_series_matrix("GSE1")
    n = len(session.calls)
    client.fetch_series_matrix("GSE1", force=True)
    assert len(session.calls) > n


def test_no_part_file_left_behind(tmp_path):
    client, _ = _client_with_dataset(tmp_path)
    client.fetch_series_matrix("GSE1")
    assert not list(tmp_path.rglob("*.part"))


def test_dataset_not_found(tmp_path):
    client = GeoClient(tmp_path, session=FakeSession({}))
    with pytest.raises(GeoNotFoundError):
        client.fetch_series_matrix("GSE999999")


def test_network_down(tmp_path):
    base = series_directory_url("GSE1") + "/matrix/"
    session = FakeSession({base: requests.ConnectionError("offline")})
    with pytest.raises(GeoNetworkError):
        GeoClient(tmp_path, session=session).fetch_series_matrix("GSE1")


def test_http_403_gives_geo_error(tmp_path):
    base = series_directory_url("GSE1") + "/matrix/"
    client = GeoClient(tmp_path, session=FakeSession({base: FakeResponse(status=403)}))
    with pytest.raises(GeoNetworkError):
        client.fetch_series_matrix("GSE1")
