"""Client für öffentliche NCBI-GEO-Daten.

Lädt die "Series Matrix"-Dateien eines GEO-Datensatzes (GSExxxx) vom
offiziellen NCBI-Server. Eine Series-Matrix-Datei enthält in EINER Datei
sowohl die Sample-Metadaten (Zeilen mit '!') als auch die Expressionstabelle.

Wichtige Eigenschaften:
- Caching: Ist der Datensatz schon in data/raw/<GSE>/ vorhanden, wird nichts geladen.
- Sicheres Schreiben: erst in .part-Datei, dann umbenennen. So bleibt nach
  einem Abbruch keine halbe Datei als "gültiger Cache" zurück.
- Verständliche Fehler statt rohen Tracebacks.
"""
import logging
import re
import time
from pathlib import Path

import requests

import config
from src.utils import ensure_dir

logger = logging.getLogger(__name__)

GEO_ID_PATTERN = re.compile(r"^GSE\d+$")
PLATFORM_ID_PATTERN = re.compile(r"^GPL\d+$")
MATRIX_FILE_PATTERN = re.compile(r'href="(GSE\d+(?:-GPL\d+)?_series_matrix\.txt\.gz)"')


class GeoError(Exception):
    """Basisklasse für alle GEO-bezogenen Fehler."""


class GeoInvalidIdError(GeoError):
    """Die übergebene ID sieht nicht wie eine GEO-Series-ID aus."""


class GeoNotFoundError(GeoError):
    """Datensatz oder Datei existiert bei GEO nicht."""


class GeoNetworkError(GeoError):
    """Server nicht erreichbar (Internet, Timeout, Serverfehler)."""


def normalize_geo_id(raw_id: str) -> str:
    """Prüft eine GEO-ID und gibt sie normalisiert zurück ('gse12220' -> 'GSE12220')."""
    geo_id = raw_id.strip().upper()
    if not GEO_ID_PATTERN.match(geo_id):
        raise GeoInvalidIdError(
            f"'{raw_id}' ist keine gültige GEO-Series-ID. Erwartet wird z. B. 'GSE12220'."
        )
    return geo_id


def normalize_platform_id(raw_id: str) -> str:
    """Prüft eine GEO-Plattform-ID ('gpl2529' -> 'GPL2529')."""
    platform_id = raw_id.strip().upper()
    if not PLATFORM_ID_PATTERN.match(platform_id):
        raise GeoInvalidIdError(
            f"'{raw_id}' ist keine gültige GEO-Plattform-ID. Erwartet wird z. B. 'GPL2529'.")
    return platform_id


def platform_query_url(platform_id: str) -> str:
    """URL, die die Plattform samt Annotationstabelle als Text liefert."""
    return (f"{config.GEO_QUERY_URL}?acc={normalize_platform_id(platform_id)}"
            "&targ=self&form=text&view=full")


def series_directory_url(geo_id: str) -> str:
    """Baut die Basis-URL eines Datensatzes.

    NCBI gruppiert Datensätze in Ordner, bei denen die letzten drei Ziffern
    durch 'nnn' ersetzt werden: GSE12220 -> GSE12nnn/GSE12220.
    """
    geo_id = normalize_geo_id(geo_id)
    stub = geo_id[:-3] + "nnn"
    return f"{config.GEO_FTP_BASE}/{stub}/{geo_id}"


def parse_matrix_filenames(directory_html: str) -> list[str]:
    """Extrahiert die Series-Matrix-Dateinamen aus einer Ordner-Auflistung (HTML)."""
    return sorted(set(MATRIX_FILE_PATTERN.findall(directory_html)))


class GeoClient:
    """Lädt GEO-Datensätze und verwaltet den lokalen Cache."""

    def __init__(self, raw_dir: Path = config.RAW_DIR, session=None) -> None:
        self.raw_dir = raw_dir
        # Die Session kann in Tests durch ein Fake-Objekt ersetzt werden.
        self.session = session or requests.Session()

    # --- Cache --------------------------------------------------------------
    def dataset_dir(self, geo_id: str) -> Path:
        return self.raw_dir / normalize_geo_id(geo_id)

    def cached_matrix_files(self, geo_id: str) -> list[Path]:
        """Bereits lokal gespeicherte Series-Matrix-Dateien (leer, wenn keine)."""
        folder = self.dataset_dir(geo_id)
        if not folder.exists():
            return []
        return sorted(folder.glob("*_series_matrix.txt.gz"))

    # --- Netzwerk -----------------------------------------------------------
    def _get(self, url: str, stream: bool = False) -> requests.Response:
        """GET mit Wiederholungen. Übersetzt Fehler in GeoError-Klassen."""
        last_error: Exception | None = None
        for attempt in range(1, config.DOWNLOAD_RETRIES + 1):
            try:
                response = self.session.get(
                    url, stream=stream, timeout=config.REQUEST_TIMEOUT_SECONDS
                )
            except requests.RequestException as error:
                last_error = error
                logger.warning("Versuch %d/%d fehlgeschlagen: %s",
                               attempt, config.DOWNLOAD_RETRIES, error)
                time.sleep(min(2 ** attempt, 10) if attempt < config.DOWNLOAD_RETRIES else 0)
                continue
            if response.status_code == 404:
                raise GeoNotFoundError(f"Nicht gefunden bei GEO: {url}")
            if response.status_code >= 500:
                last_error = GeoNetworkError(f"Serverfehler {response.status_code}")
                continue
            if response.status_code >= 400:
                raise GeoNetworkError(
                    f"GEO hat die Anfrage abgelehnt (HTTP {response.status_code}): {url}")
            return response
        raise GeoNetworkError(
            "GEO ist nicht erreichbar. Prüfe deine Internetverbindung und versuche es erneut. "
            f"(Letzter Fehler: {last_error})"
        )

    def list_matrix_files(self, geo_id: str) -> list[str]:
        """Fragt bei GEO ab, welche Series-Matrix-Dateien es gibt.

        Mehrere Dateien gibt es, wenn ein Datensatz mehrere Plattformen (GPL...) nutzt.
        """
        url = series_directory_url(geo_id) + "/matrix/"
        response = self._get(url)
        names = parse_matrix_filenames(response.text)
        if not names:
            raise GeoNotFoundError(
                f"{geo_id}: Keine Series-Matrix-Datei gefunden. "
                "Der Datensatz hat evtl. keine aufbereiteten Expressionsdaten."
            )
        return names

    def _download_file(self, url: str, destination: Path) -> None:
        """Lädt eine Datei streamend herunter und schreibt sie atomar."""
        partial = destination.with_suffix(destination.suffix + ".part")
        response = self._get(url, stream=True)
        try:
            with open(partial, "wb") as file:
                for chunk in response.iter_content(config.DOWNLOAD_CHUNK_BYTES):
                    file.write(chunk)
            partial.replace(destination)
        finally:
            if partial.exists():
                partial.unlink()

    # --- Öffentliche Hauptfunktion -----------------------------------------
    def fetch_series_matrix(self, geo_id: str, force: bool = False) -> list[Path]:
        """Stellt sicher, dass die Series-Matrix-Dateien lokal vorliegen.

        Args:
            geo_id: z. B. 'GSE12220'.
            force: True = vorhandenen Cache ignorieren und neu laden.

        Returns:
            Liste lokaler Dateipfade.
        """
        geo_id = normalize_geo_id(geo_id)
        cached = self.cached_matrix_files(geo_id)
        if cached and not force:
            logger.info("%s: Cache gefunden, kein Download nötig (%d Datei/en).",
                        geo_id, len(cached))
            return cached

        folder = ensure_dir(self.dataset_dir(geo_id))
        base_url = series_directory_url(geo_id)
        local_files: list[Path] = []
        for name in self.list_matrix_files(geo_id):
            target = folder / name
            logger.info("Lade %s ...", name)
            self._download_file(f"{base_url}/matrix/{name}", target)
            local_files.append(target)
        return local_files

    def fetch_platform_annotation(self, platform_id: str, force: bool = False) -> Path:
        """Lädt die Plattform-Annotation (Sonden-ID -> Gen) und cached sie.

        Gespeichert wird unter data/raw/<GPL-ID>/<GPL-ID>.soft.
        """
        platform_id = normalize_platform_id(platform_id)
        target = self.raw_dir / platform_id / f"{platform_id}.soft"
        if target.exists() and not force:
            logger.info("%s: Annotation aus dem Cache.", platform_id)
            return target
        ensure_dir(target.parent)
        logger.info("Lade Plattform-Annotation %s ...", platform_id)
        self._download_file(platform_query_url(platform_id), target)
        return target
