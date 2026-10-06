"""Zentrale Konfiguration der ROS Research Suite.

Alle Pfade und (später) alle Schwellenwerte stehen hier, nicht verstreut im Code.
"""
from pathlib import Path

SOFTWARE_VERSION = "0.1.0"

# --- Pfade -----------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"
REPORTS_DIR = RESULTS_DIR / "reports"

# --- NCBI GEO --------------------------------------------------------------
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
REQUEST_TIMEOUT_SECONDS = 60
DOWNLOAD_RETRIES = 3
DOWNLOAD_CHUNK_BYTES = 1024 * 256

# --- Schwellenwerte für spätere Phasen (Differential Expression) -----------
MIN_ABS_LOG2FC = 1.0
MAX_ADJUSTED_P = 0.05

# --- Plattform-Annotation (Phase 3) ----------------------------------------
GEO_QUERY_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
# Spaltennamen der Annotationstabelle. None = automatisch erkennen.
# Falls die Erkennung falsch liegt, trägst du hier den exakten Spaltennamen ein.
ANNOTATION_SYMBOL_COLUMN = None
ANNOTATION_SYSTEMATIC_COLUMN = None
# Mindestanteil der Sonden, der in der Annotation gefunden werden muss
MIN_ANNOTATION_MATCH_RATE = 0.9

# Die Plattform GPL2529 enthält Sonden für S. cerevisiae UND S. pombe.
# Nur Sonden dieses Organismus werden in den Ergebnistabellen behalten.
# None = kein Filter.
ORGANISM_FILTER = "Saccharomyces cerevisiae"

# --- Phase 4: Differential Expression ohne Replikate -------------------------
# "Responsiv" = |log2FC| >= MIN_ABS_LOG2FC in mindestens so vielen
# aufeinanderfolgenden Zeitpunkten (filtert Einzelwert-Ausreißer bei n=1).
MIN_CONSECUTIVE_TIMEPOINTS = 2
TOP_N_GENES = 20
CONTROL_PROBE_PREFIX = "AFFX-"   # Affymetrix-Kontrollsonden, keine Gene

# --- Phase 5: Zeitverlauf-Muster -----------------------------------------------
# Früh = Beginn der Antwort (erster Zeitpunkt der ersten Folge über der Schwelle) bis einschließlich
# dieser Minute, sonst spät.
EARLY_ONSET_MAX_MIN = 60
# Transient = am letzten Zeitpunkt unter der Schwelle UND höchstens dieser Bruchteil des Maximums.
TRANSIENT_RETURN_FRACTION = 0.5
# Sensitivitätsanalyse: 0.58 = 1.5x, 1.0 = 2x, 1.58 = 3x, 2.0 = 4x Veränderung
SENSITIVITY_THRESHOLDS = (0.58, 1.0, 1.58, 2.0)
SENSITIVITY_CONSECUTIVE = (1, 2, 3)

# --- Phase 6: ROS-Gene und Unbiased Discovery -----------------------------------
TOP_N_DISCOVERY = 50   # so viele der stärksten Antworter werden mit der ROS-Liste verglichen

# --- Phase 6b: Standards (Robustheit, Spezifität, Validierung) -------------------
# Spezifisch = in der einen Behandlung responsiv UND in der anderen max. |log2FC| < LOW_FRACTION * Schwelle.
# Dazwischen liegt eine Grauzone ("ambiguous"). Wert ist eine Konvention (Sensitivität: Grid).
SPECIFICITY_LOW_FRACTION = 0.5
SPECIFICITY_LOW_FRACTION_GRID = (0.3, 0.5, 0.7)
SPECIFICITY_THRESHOLD_GRID = (0.58, 1.0, 1.58)
# Robustheits-Score = Anteil der Parameterkombinationen (Schwelle x Zeitpunkte-in-Folge), in denen
# ein Gen responsiv bleibt. Klassengrenzen sind Konventionen.
ROBUSTNESS_ROBUST_MIN = 0.75
ROBUSTNESS_FRAGILE_MAX = 0.40
# Stabilität der Zeitverlaufs-Muster gegenüber den zwei Konventions-Parametern
PATTERN_EARLY_GRID = (30, 60, 100)
PATTERN_RETURN_GRID = (0.3, 0.5, 0.7)
# Positivkontrollen: erwartete Richtung unter H2O2 (Quelle muss vom Benutzer belegt werden)
