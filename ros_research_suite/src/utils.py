"""Kleine Hilfsfunktionen, die von mehreren Modulen genutzt werden."""
import logging
from pathlib import Path


def setup_logging(verbose: bool = False) -> None:
    """Richtet einheitliches Logging für das ganze Programm ein."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )


def ensure_dir(path: Path) -> Path:
    """Legt einen Ordner (inkl. Elternordner) an, falls er fehlt."""
    path.mkdir(parents=True, exist_ok=True)
    return path
