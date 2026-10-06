import sys
from pathlib import Path

# Projektwurzel auf den Importpfad, damit 'import config' und 'import src...' funktionieren
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
