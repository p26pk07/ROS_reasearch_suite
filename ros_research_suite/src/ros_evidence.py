"""Belege für die kuratierte Genliste und Positivkontrollen (wissenschaftlicher Standard 17/20).

1. GO-Abgleich: Passt zu jedem Gen der Liste mindestens ein GO-Begriff der Plattform-Annotation
   zur Kategorie? 'supported' = ja, 'no_matching_term' = bitte manuell prüfen (kein Beweis für
   einen Fehler, die Annotation stammt von 2009 und der Abgleich ist eine Stichwortsuche).
2. Positivkontrollen: Stimmt die beobachtete Richtung unter H2O2 mit der erwarteten überein?
   Die Quellen für die Erwartung musst DU belegen; ohne Quelle steht 'Quelle fehlt'.
"""
import pandas as pd

from src.ros_genes import ROS_GENES

# Stichwörter (Kleinbuchstaben) in den GO-Begriffsnamen je Kategorie. Absichtlich sichtbar und einfach.
CATEGORY_KEYWORDS = {
    "antioxidant_defense": ("oxidative stress", "antioxidant", "superoxide", "peroxide", "glutathione",
                            "thioredoxin", "redox homeostasis", "peroxiredoxin", "catalase", "cell redox"),
    "ros_metabolism": ("nadph", "pentose", "oxidation reduction", "redox", "glucose-6-phosphate"),
    "dna_damage_response": ("dna repair", "dna damage", "excision repair", "mismatch repair",
                            "recombination", "checkpoint", "deoxyribonucleotide", "dna replication"),
    "protein_repair": ("protein folding", "heat", "chaperone", "refolding", "ubiquitin", "disulfide",
                       "protein repair", "response to stress"),
    "stress_response": ("stress", "osmotic", "heat", "trehalose"),
    "mitochondrial_response": ("mitochond", "respirat", "electron transport", "tricarboxylic",
                               "iron-sulfur", "aerobic", "citrate"),
}

# Erwartung: Richtung unter H2O2. Quelle = von dir zu belegen (SGD / Primärliteratur).
POSITIVE_CONTROLS = {
    "CTT1": ("up", ""), "SRX1": ("up", ""), "TRX2": ("up", ""), "TSA1": ("up", ""), "GPX2": ("up", ""),
}


def _go_terms(value) -> list[tuple[str, str]]:
    """'0006979 // response to oxidative stress // evidence /// ...' -> [(Name, Evidenz), ...]"""
    terms = []
    for entry in str(value).split("///"):
        parts = [p.strip() for p in entry.split("//")]
        if len(parts) >= 2 and parts[1]:
            terms.append((parts[1], parts[2] if len(parts) > 2 else ""))
    return terms


def go_support(annotation: pd.DataFrame, genes=ROS_GENES) -> pd.DataFrame:
    """Prüft jedes Gen der Liste gegen die GO-Begriffe der Annotation (siehe Modul-Doku)."""
    go_columns = [c for c in annotation.columns if "gene ontology" in str(c).lower()]
    if not go_columns:
        raise ValueError("Die Annotation enthält keine GO-Spalten.")
    by_orf = annotation.groupby("ORF").groups if "ORF" in annotation.columns else {}
    rows = []
    for gene in genes:
        probes = list(by_orf.get(gene.orf, []))
        terms = []
        for probe in probes:
            for column in go_columns:
                terms += _go_terms(annotation.loc[probe, column])
        keywords = CATEGORY_KEYWORDS[gene.category]
        matching = sorted({name for name, _ in terms if any(k in name.lower() for k in keywords)})
        rows.append({"symbol": gene.symbol, "orf": gene.orf, "category": gene.category,
                     "n_go_terms": len({name for name, _ in terms}),
                     "n_matching_terms": len(matching), "example_terms": "; ".join(matching[:3]),
                     "go_support": ("not_in_annotation" if not probes else
                                    "supported" if matching else "no_matching_term")})
    return pd.DataFrame(rows)


def validate_positive_controls(ros_table: pd.DataFrame, treatment: str = "H2O2",
                               controls: dict | None = None) -> pd.DataFrame:
    """Vergleicht erwartete und beobachtete Richtung der Positivkontrollen."""
    controls = POSITIVE_CONTROLS if controls is None else controls
    best = ros_table[ros_table["is_best_probe"]].set_index("symbol")
    rows = []
    for symbol, (expected, source) in controls.items():
        if symbol not in best.index:
            rows.append({"symbol": symbol, "expected": expected, "status": "not_on_chip",
                         "source": source or "Quelle fehlt"})
            continue
        row = best.loc[symbol]
        observed = row[f"{treatment}_direction"]
        rows.append({"symbol": symbol, "expected": expected, "observed": observed,
                     "peak_log2fc": round(float(row[f"{treatment}_peak_log2fc"]), 2),
                     "status": "agrees" if observed == expected and row[f"{treatment}_responsive"]
                     else "does_not_agree", "source": source or "Quelle fehlt"})
    return pd.DataFrame(rows)
