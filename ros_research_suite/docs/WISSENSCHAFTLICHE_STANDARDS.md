# Wissenschaftliche Standards der ROS Research Suite

Diese Datei ist verbindlich für jede weitere Phase. Wo eine Anforderung für GSE12220
nicht valide umsetzbar ist, steht hier, was stattdessen gilt und warum.

## 0. Oberste Regel
Wissenschaftliche Korrektheit hat Vorrang vor visueller Komplexität. Kein Test, keine
Transformation und kein Filter wird eingebaut, nur weil er beeindruckend aussieht.
Ist eine Methode wegen fehlender Replikate oder Metadaten nicht valide, wird sie nicht
erzwungen: Die Einschränkung wird gekennzeichnet und eine valide Alternative genutzt.

## 1. Forschungsfrage
Wie reagiert *Saccharomyces cerevisiae* auf H2O2-induzierten oxidativen Stress auf Ebene der
Genexpression, wie verändert sich diese Reaktion über die Zeit, und welche Unterschiede
bestehen gegenüber der MMS-induzierten Stressantwort?

## 2. Was die Daten von GSE12220 erlauben (Datenlage)
| Eigenschaft | Befund (aus den Metadaten) | Folge |
|---|---|---|
| Design | H2O2 (0,3 mM) und MMS (0,1 %), je t = 0, 30, 60, 100, 140, 180 min | Zeitreihe je Behandlung |
| Replikate | keine (n = 1 je Behandlung und Zeitpunkt) | keine Varianzschätzung je Gen, **keine Gen-p-Werte, keine Gen-FDR** |
| Kontrolle | nur t = 0 vor der Behandlung, keine zeitgleiche unbehandelte Probe | Zeiteffekte (z. B. Wachstum) nicht abtrennbar |
| Datenebene | von den Autoren RMA-normalisierte, log2-skalierte Werte (**processed data**) | Die Pipeline normalisiert nicht selbst; CEL-Rohdateien werden nicht verwendet |
| Plattform | GPL2529, Affymetrix; Annotation vom 13.03.2009; enthält Kontroll- und Fremdorganismus-Sonden | Filter: AFFX-Sonden und andere Organismen entfernt |

## 3. Sprachregeln
- **Responsive** (nicht „signifikant“): |log2FC| >= Schwelle in >= k aufeinanderfolgenden Zeitpunkten.
  Das ist ein Kriterium für Effektstärke + Konsistenz, **kein Signifikanztest**. In Tabellen steht
  `statistical_evidence = not_testable_no_replicates`.
- **Signifikant** darf nur heißen: Test + FDR-Korrektur (Benjamini-Hochberg) wurden gerechnet.
  Das ist nur für Datensätze mit Replikaten möglich.
- **Spezifisch** (H2O2- oder MMS-spezifisch): starke Antwort in einer Behandlung UND geringe Antwort in der
  anderen (|log2FC| deutlich unter der Schwelle an allen Zeitpunkten). Dazwischen liegt eine Grauzone
  („ambig“), die nicht klassifiziert wird.
- **Neu** wird nie vor dem Literaturvergleich behauptet. Eigenleistung = reproduzierbare, zeitaufgelöste,
  vergleichende Reanalyse, keine Entdeckung.
- Jede Aussage trägt eine Stufe: MEASURED (Messwert), CALCULATED (Berechnung), INFERRED (Schluss aus Daten),
  HYPOTHESIS (Vermutung). Genexpression ist keine ROS-Konzentration, keine Proteinmenge, keine Enzymaktivität.

## 4. Mathematische Definitionen (werden im Code per Test geprüft)
- Eingabe: E = log2-skalierter Expressionswert (RMA). Daher ist die Differenz der Werte bereits der log2-Fold-Change.
- **M** = log2(E_behandelt) − log2(E_Baseline) = log2FC (Baseline = t = 0 derselben Behandlung)
- **A** = (log2(E_behandelt) + log2(E_Baseline)) / 2
- **Responsive**: Folge von >= k Zeitpunkten mit M >= +T (hoch) bzw. M <= −T (runter)
- **Muster**: Beginn *early* (erste Folge startet <= 60 min) oder *late*; Dauer *sustained* (am letzten Zeitpunkt >= T),
  *transient* (am Ende < T und <= 50 % des Maximums), *partial_decline*, *biphasic*
- Für ein Beispielgen wird automatisch eine vollständige Rechnung ausgegeben (Phase 6b).

## 5. Parameter und ihre Begründung
| Parameter | Wert | Status der Begründung |
|---|---|---|
| T (`MIN_ABS_LOG2FC`) | 1,0 | Konvention (2-fache Änderung); Sensitivität 0,58 / 1,0 / 1,58 / 2,0 |
| k (`MIN_CONSECUTIVE_TIMEPOINTS`) | 2 | Schutz vor Einzelwert-Ausreißern bei n = 1; Sensitivität 1 / 2 / 3 |
| `EARLY_ONSET_MAX_MIN` | 60 | **Konvention, nicht datengestützt** – Sensitivität: Raster 30 / 60 / 100 (`--analysis robustness`) |
| `TRANSIENT_RETURN_FRACTION` | 0,5 | **Konvention, nicht datengestützt** – Sensitivität: Raster 0,3 / 0,5 / 0,7 (`--analysis robustness`) |
| `SPECIFICITY_LOW_FRACTION` | 0,5 | Konvention („gering“ = < 0,5 × Schwelle) – Sensitivität: 0,3 / 0,5 / 0,7 |
| FDR-Schwelle | 0,05 | gilt nur bei Replikaten (Validierungsdatensatz) |

## 6. Hypothesenbasiert und unvoreingenommen getrennt
- **Hypothesenbasiert**: kuratierte ROS-/Stress-Genliste (`ros_genes.py`).
  *Offenlegung:* Die Liste wurde zusammengestellt, nachdem erste Top-Gene (u. a. SOL4, OYE3, DDR2, HSP30)
  schon in der Ausgabe zu sehen waren. Dass die Liste unter den stärksten Antwortern angereichert ist, ist
  deshalb **kein unabhängiger Beleg**. Abhilfe in Phase 6b: Begründung je Gen aus einer Datenbank
  (GO-Annotation der Plattform) statt aus dem Gedächtnis; Trennung in a-priori-Positivkontrollen und übrige Liste.
- **Unvoreingenommen**: alle Gene ohne Vorauswahl (Top-Antworter, Cluster, Anreicherung).

## 7. Quellenregel
Jede biologische Aussage braucht mindestens eine Quelle: Datensatz, reproduzierbare Berechnung, etablierte
Datenbank (mit Version/Datum) oder Primärliteratur. Literaturangaben werden nur eingetragen, wenn sie aus den
GEO-Metadaten stammen oder von dir geprüft wurden. Die KI erfindet keine Referenzen.

## 8. Umsetzungsstand der 22 Anforderungen
| # | Anforderung | Stand | Umsetzung / Anpassung |
|---|---|---|---|
| 1 | Forschungsfrage | erfüllt | steht oben, wird Titel des Reports |
| 2 | Reproduzierbarkeit | erfüllt (6b) | `results/<GSE>_analysis_metadata.json`: Samples, Plattform, Filter mit Sondenzahlen, Parameter, Methoden, Versionen, Annotationsdatum; wird bei jedem Lauf ergänzt |
| 3 | Raw vs. processed | erfüllt (6b) | Metadaten: `data_level = processed (RMA durch die Autoren)`; Report folgt in Phase 9 |
| 4 | Mathematische Transparenz | erfüllt (6b) | `statistics.ma_values` ist die einzige Stelle für M und A; Tests prüfen M = log2(Verhältnis) und Gleichheit mit der log2FC-Tabelle; `--analysis validation` gibt eine Beispielrechnung aus |
| 5 | Responsive definieren | angepasst, umgesetzt (6b) | Effektstärke + Konsistenz; Spalte `<Behandlung>_statistical_evidence` weist „nicht testbar“ aus |
| 6 | Multiple Testing | angepasst | Welch-Test + Benjamini-Hochberg gebaut und getestet, **nicht anwendbar auf GSE12220**, aktiv bei Replikaten |
| 7 | Zeitreihen ausbauen | teilweise | regelbasierte Muster vorhanden; Clustering mit begründetem k folgt (Phase 7b) |
| 8 | H2O2 vs. MMS | teilweise | Spezifität mit Grauzone umgesetzt (6b: `specificity_*`); UpSet/Heatmap folgen in Phase 8 |
| 9 | Hypothesen- vs. unvoreingenommen | erfüllt, mit Vorbehalt | beide Analysen vorhanden; Zirkularität der Liste offengelegt (Abschnitt 6) |
| 10 | GO / Pathways | geplant | Phase 7, lokale GO-Annotation der Plattform, Datum dokumentiert |
| 11 | GSEA / rangbasiert | geplant | Phase 7; wegen n = 1 Permutation über Gene, nicht über Samples (Limitation wird ausgewiesen) |
| 12 | Zeitliche Pathway-Analyse | geplant | Phase 7, Aussagen erst nach der Analyse |
| 13 | Netzwerk | geplant | später, mit Quelle und Evidenzart (experimentell vs. vorhergesagt) |
| 14/15 | Robustheit / Sensitivität | erfüllt (6b) | `--analysis robustness`: Score je Gen (12 Kombinationen), Muster-Stabilität gegenüber den 2 Konventions-Parametern, Spezifitäts-Stabilität |
| 16 | Nullergebnisse zulassen | teilweise | Ausgaben melden „keine responsiven Gene“; durchgängige Regel im Report (Phase 9) |
| 17 | Validierung / Positivkontrolle | teilweise (6b) | Mechanismus vorhanden (`--analysis validation`); die Quellen für die Erwartung musst du eintragen (`ros_evidence.POSITIVE_CONTROLS`) |
| 18 | Zweiter Datensatz | geplant | Datensatz mit Replikaten suchen; dort echte Tests + FDR |
| 19 | Facharbeit getrennt | geplant | OD600 als eigene Evidenzebene, nie mit Genexpression verrechnet |
| 20 | Keine KI-Halluzination | teilweise | Regeln hier festgelegt; Umsetzung im KI-Assistenten (Phase 10) |
| 21 | Report | geplant | Methods / Results / Interpretation / Limitations (Phase 9) |
| 22 | Neuheit | erfüllt | Eigenleistung wie in Abschnitt 3 formuliert |
