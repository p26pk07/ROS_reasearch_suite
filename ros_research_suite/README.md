# ROS Research Suite

Bioinformatik-Projekt zur computergestützten Untersuchung von **oxidativem Stress (H2O2) in
*Saccharomyces cerevisiae*** anhand öffentlicher Genexpressionsdaten (NCBI GEO, Start: GSE12220).

> **Stand: Phase 6b von 12** – wie Phase 6, plus Anpassung an die wissenschaftlichen Standards (siehe `docs/WISSENSCHAFTLICHE_STANDARDS.md`): Metadaten-JSON, Beispielrechnung, Spezifität mit Grauzone, Robustheits-Scores, GO-Abgleich der Genliste, Positivkontrollen.
> Analysen, Plots und Report folgen in den nächsten Phasen.

## Grundregeln
- Es werden nur tatsächlich heruntergeladene Daten verwendet. Nichts wird erfunden.
- Heruntergeladene Rohdaten liegen unter `data/raw/<GSE-ID>/` und werden wiederverwendet.

## Installation (VS Code, Terminal)
```bash
cd ros_research_suite
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # macOS / Linux
pip install -r requirements.txt
```

## Starten
```bash
python main.py --dataset GSE12220                    # laden (oder Cache) + Übersicht
python main.py --dataset GSE12220 --analysis log2fc  # log2FC je Zeitpunkt gegen t=0
python main.py --dataset GSE12220 --analysis differential  # responsive Gene, Top-Listen, H2O2 vs MMS, Grafiken
python main.py --dataset GSE12220 --analysis timecourse    # Zeitverlauf-Muster, Sensitivität, Verlaufsplots
python main.py --dataset GSE12220 --analysis ros           # bekannte ROS-Gene + Top-Antworter außerhalb der Liste
python main.py --dataset GSE12220 --analysis robustness    # Robustheit, Muster-Stabilität, Spezifität mit Grauzone
python main.py --dataset GSE12220 --analysis validation    # Beispielrechnung, GO-Abgleich, Positivkontrollen
python main.py --dataset GSE12220 --inspect-annotation  # Spalten der Gen-Annotation anzeigen
python main.py --dataset GSE12220 --force-download   # Cache ignorieren
python main.py --dataset GSE12220 --reconfirm        # Zuordnung neu bestätigen
python main.py --dataset GSE12220 --yes              # Zuordnung ohne Nachfrage (nur wenn eindeutig)
python main.py --dataset GSE12220 --verbose          # mehr Logging
python -m pytest                                     # Tests
```

## Struktur
| Pfad | Zweck |
|---|---|
| `config.py` | Pfade, URLs, Schwellenwerte |
| `main.py` | Kommandozeile |
| `src/geo_client.py` | Download + Cache + Fehlerklassen |
| `src/metadata.py` | Series-Matrix einlesen, Sample-Tabelle, Expressionstabelle |
| `src/statistics.py` | mean/median/sd/se, (log2-)Fold-Change, Welch-Test, Benjamini-Hochberg, log2FC gegen Baseline |
| `src/provenance.py` | analysis_metadata.json (Samples, Datenebene, Filter, Parameter, Versionen), Beispielrechnung |
| `src/robustness.py` | Robustheit responsiver Gene, Muster-Stabilität, Spezifitäts-Sensitivität |
| `src/ros_evidence.py` | GO-Abgleich der Genliste, Positivkontrollen |
| `src/ros_genes.py` | kuratierte ROS-/Stress-Genliste (Vorwissen, kein Beweis) + Zuordnung zu Sonden |
| `src/ros_analysis.py` | Kategorien-Auswertung, Top-Antworter, Vergleich Liste vs. Daten |
| `src/timecourse.py` | Muster: früh/spät, dauerhaft/transient/biphasisch, hoch/runter |
| `src/differential_expression.py` | responsive Gene (Stärke + Konsistenz), Top-Listen, Vergleich H2O2/MMS |
| `src/visualization.py` | MA-Plots, Antwortstärke-Plot (results/figures/) |
| `src/annotation.py` | Sonden-ID -> Gen (GEO-Plattform-Annotation) |
| `src/design.py` | Behandlung/Zeit erkennen, Bestätigung, manuelle Korrektur per CSV |
| `tests/` | pytest-Tests (ohne Internet) |

## Limitationen
- GSE12220 hat keine Replikate und keine unbehandelte Kontrolle: Referenz ist t=0, daher keine klassischen p-Werte.
- Die Zuordnung wird per Textmuster erkannt und muss von dir bestätigt werden.
- Genexpression ≠ ROS-Konzentration ≠ Proteinmenge ≠ Enzymaktivität. Spätere Ergebnisse sind entsprechend zu lesen.

## Datenquelle
NCBI Gene Expression Omnibus (GEO), https://www.ncbi.nlm.nih.gov/geo/
