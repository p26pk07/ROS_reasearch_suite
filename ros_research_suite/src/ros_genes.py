"""Kuratierte Liste ROS-/Stress-relevanter Gene der Bäckerhefe (S. cerevisiae).

WICHTIG: Diese Liste ist Vorwissen, KEIN Beweis. Sie dient nur dazu, die
datengetriebenen Ergebnisse biologisch einzuordnen ("Wie verhalten sich bekannte
antioxidative Gene?"). Welche Gene tatsächlich reagieren, entscheiden allein die Daten.

Jeder Eintrag: Symbol (+ optional frühere Symbole), systematischer Name (ORF, z. B. YGR088W), Kategorie, Untergruppe.
Die Zuordnung zu den Chip-Sonden läuft über den systematischen Namen; das Programm
prüft dabei, ob das Symbol der Annotation zu deinem Eintrag passt, und meldet
Abweichungen und Gene, die auf dem Chip fehlen. Tippfehler fallen so auf.
Jedes Gen steht nur in EINER (der wichtigsten) Kategorie. Stichprobenartig gegen die
Saccharomyces Genome Database (yeastgenome.org) prüfen.
"""
import pandas as pd
from dataclasses import dataclass


class RosGeneError(Exception):
    """ROS-Analyse nicht möglich (z. B. fehlende Gennamen-Annotation)."""


@dataclass(frozen=True)
class RosGene:
    symbol: str
    orf: str
    category: str
    group: str
    aliases: tuple[str, ...] = ()   # frühere/alternative Symbole, die die Annotation benutzen könnte


CATEGORY_LABELS = {
    "antioxidant_defense": "Antioxidant defense",
    "ros_metabolism": "ROS production / metabolism (NADPH supply, redox enzymes)",
    "dna_damage_response": "DNA damage response",
    "protein_repair": "Protein repair / chaperones",
    "stress_response": "General stress response",
    "mitochondrial_response": "Mitochondrial response",
}


def _genes(category: str, group: str, *entries: tuple[str, ...]) -> list[RosGene]:
    """Einträge: (Symbol, ORF) oder (Symbol, ORF, 'Alias1,Alias2')."""
    genes = []
    for symbol, orf, *rest in entries:
        aliases = tuple(rest[0].split(",")) if rest else ()
        genes.append(RosGene(symbol, orf, category, group, aliases))
    return genes


ROS_GENES: tuple[RosGene, ...] = tuple(
    # --- Antioxidant defense ---------------------------------------------------
    _genes("antioxidant_defense", "superoxide dismutase",
           ("SOD1", "YJR104C"), ("SOD2", "YHR008C"), ("CCS1", "YMR038C"))
    + _genes("antioxidant_defense", "catalase", ("CTT1", "YGR088W"), ("CTA1", "YDR256C"))
    + _genes("antioxidant_defense", "peroxiredoxin",
             ("TSA1", "YML028W"), ("TSA2", "YDR453C"), ("AHP1", "YLR109W"),
             ("DOT5", "YIL010W"), ("PRX1", "YBL064C"))
    + _genes("antioxidant_defense", "glutathione peroxidase",
             ("GPX1", "YKL026C"), ("GPX2", "YBR244W"), ("GPX3", "YIR037W", "HYR1,ORP1"))
    + _genes("antioxidant_defense", "glutathione system",
             ("GSH1", "YJL101C"), ("GSH2", "YOL049W"), ("GLR1", "YPL091W"),
             ("GRX1", "YCL035C"), ("GRX2", "YDR513W"), ("GTT1", "YIR038C"), ("GTT2", "YLL060C"))
    + _genes("antioxidant_defense", "thioredoxin system",
             ("TRX1", "YLR043C"), ("TRX2", "YGR209C"), ("TRX3", "YCR083W"),
             ("TRR1", "YDR353W"), ("TRR2", "YHR106W"), ("SRX1", "YKL086W"))
    # --- ROS production / metabolism (Hefe hat keine klassische NADPH-Oxidase) ---
    + _genes("ros_metabolism", "NADPH supply / pentose phosphate pathway",
             ("ZWF1", "YNL241C"), ("GND1", "YHR183W"), ("GND2", "YGR256W"),
             ("SOL3", "YHR163W"), ("SOL4", "YGR248W"), ("TKL1", "YPR074C"),
             ("TAL1", "YLR354C"), ("ALD6", "YPL061W"), ("IDP2", "YLR174W"),
             ("POS5", "YPL188W"))
    + _genes("ros_metabolism", "NADPH oxidoreductase (old yellow enzymes)",
             ("OYE2", "YHR179W"), ("OYE3", "YPL171C"))
    # --- DNA damage response -------------------------------------------------------
    + _genes("dna_damage_response", "checkpoint / signaling",
             ("MEC1", "YBR136W"), ("RAD53", "YPL153C"), ("RAD9", "YDR217C"), ("DUN1", "YDL101C"))
    + _genes("dna_damage_response", "ribonucleotide reductase",
             ("RNR1", "YER070W"), ("RNR2", "YJL026W"), ("RNR3", "YIL066C"), ("RNR4", "YGR180C"))
    + _genes("dna_damage_response", "recombination / post-replication repair",
             ("RAD51", "YER095W"), ("RAD52", "YML032C"), ("RAD6", "YGL058W"), ("RAD18", "YCR066W"))
    + _genes("dna_damage_response", "base excision / alkylation repair",
             ("OGG1", "YML060W"), ("APN1", "YKL114C"), ("NTG1", "YAL015C"), ("NTG2", "YOL043C"),
             ("MAG1", "YER142C"), ("MGT1", "YDL200C"), ("RAD27", "YKL113C"), ("UNG1", "YML021C"))
    + _genes("dna_damage_response", "nucleotide excision / mismatch repair",
             ("RAD2", "YGR258C"), ("RAD14", "YMR201C"), ("RAD23", "YEL037C"),
             ("MSH2", "YOL090W"), ("MLH1", "YMR167W"))
    + _genes("dna_damage_response", "damage-responsive gene", ("DDR2", "YOL052C-A"))
    # --- Protein repair / chaperones --------------------------------------------------
    + _genes("protein_repair", "methionine sulfoxide reductase", ("MXR1", "YER042W"), ("MXR2", "YCL033C"))
    + _genes("protein_repair", "heat shock proteins / chaperones",
             ("HSP104", "YLL026W"), ("HSP26", "YBR072W"), ("HSP12", "YFL014W"),
             ("HSP42", "YDR171W"), ("HSP30", "YCR021C"), ("SSA1", "YAL005C"),
             ("SSA4", "YER103W"), ("HSP82", "YPL240C"), ("HSC82", "YMR186W"), ("STI1", "YOR027W"))
    + _genes("protein_repair", "ubiquitin / oxidative protein folding",
             ("UBI4", "YLL039C"), ("PDI1", "YCL043C"), ("ERO1", "YML130C"))
    # --- General stress response ----------------------------------------------------------
    + _genes("stress_response", "transcription factors",
             ("YAP1", "YML007W"), ("SKN7", "YHR206W"), ("MSN2", "YMR037C"),
             ("MSN4", "YKL062W"), ("HSF1", "YGL073W"), ("STB5", "YHR178W"))
    + _genes("stress_response", "signaling", ("HOG1", "YLR113W"), ("PBS2", "YJL128C"))
    + _genes("stress_response", "trehalose metabolism",
             ("TPS1", "YBR126C"), ("TPS2", "YDR074W"), ("TSL1", "YML100W"), ("NTH1", "YDR001C"))
    + _genes("stress_response", "stress-induced genes",
             ("GRE1", "YPL223C"), ("GRE2", "YOL151W"), ("SIP18", "YMR175W"),
             ("XBP1", "YIL101C"), ("YGP1", "YNL160W"))
    # --- Mitochondrial response -------------------------------------------------------------
    + _genes("mitochondrial_response", "oxidative-stress-related mitochondrial enzymes",
             ("CCP1", "YKR066C"), ("ACO1", "YLR304C"), ("ACO2", "YJL200C"), ("MDH1", "YKL085W"),
             ("SDH1", "YKL148C"), ("CIT1", "YNR001C"))
    + _genes("mitochondrial_response", "respiration / regulation",
             ("CYC1", "YJR048W"), ("CYC7", "YEL039C"), ("HAP4", "YKL109W"))
    + _genes("mitochondrial_response", "mitochondrial quality control / Fe-S",
             ("PIM1", "YBL022C"), ("YME1", "YPR024W"), ("HSP78", "YDR258C"),
             ("NFS1", "YCL017C"), ("ISU1", "YPL135W"))
)


def known_orfs(genes=ROS_GENES) -> set[str]:
    """Alle systematischen Namen der kuratierten Liste."""
    return {gene.orf for gene in genes}


def build_ros_table(result: pd.DataFrame, patterns: pd.DataFrame | None = None,
                    genes=ROS_GENES) -> tuple[pd.DataFrame, list[RosGene]]:
    """Ordnet die kuratierte Liste den Sonden zu und hängt die Analyse-Ergebnisse an.

    Args:
        result: Ergebnis von differential_expression.classify_all (mit Gennamen-Spalten).
        patterns: optional, Ergebnis von timecourse.classify_patterns_all.
    Returns:
        (Tabelle: eine Zeile pro gefundener Sonde, Liste der Gene ohne Sonde auf dem Chip)

    Spalten u. a.: symbol, orf, category, group, annotation_symbol, symbol_check
    ('ok' / 'annotation_has_no_symbol' / 'symbol_mismatch'), n_probes, is_best_probe.
    Hat ein Gen mehrere Sonden, gilt die mit dem stärksten |log2FC| als beste Sonde.
    """
    if not {"systematic_name", "gene_symbol"} <= set(result.columns):
        raise RosGeneError("Für die ROS-Analyse wird die Gennamen-Annotation benötigt "
                           "(Spalten systematic_name/gene_symbol fehlen). Internet/GEO prüfen.")
    probes_by_orf = result.groupby("systematic_name").groups
    treatments = [c[:-len("_peak_log2fc")] for c in result.columns if c.endswith("_peak_log2fc")]

    rows, missing = [], []
    for gene in genes:
        probes = list(probes_by_orf.get(gene.orf, []))
        if not probes:
            missing.append(gene)
            continue
        for probe in probes:
            annotation_symbol = result.loc[probe, "gene_symbol"]
            if pd.isna(annotation_symbol):
                check = "annotation_has_no_symbol"
            elif {gene.symbol, *gene.aliases} & {part.strip() for part in str(annotation_symbol).split("///")}:
                check = "ok"
            else:
                check = "symbol_mismatch"
            rows.append({"probe_id": probe, "symbol": gene.symbol, "orf": gene.orf,
                         "category": gene.category, "group": gene.group,
                         "annotation_symbol": annotation_symbol, "symbol_check": check})
    if not rows:
        raise RosGeneError("Keines der kuratierten Gene wurde auf dem Chip gefunden.")

    table = pd.DataFrame(rows).set_index("probe_id")
    stats = result.drop(columns=["display_name", "gene_symbol", "systematic_name"], errors="ignore")
    table = table.join(stats)
    if patterns is not None:
        table = table.join(patterns[[c for c in patterns.columns if c.endswith("_pattern")]])

    strength = result.loc[table.index, [f"{t}_peak_log2fc" for t in treatments]].abs().max(axis=1)
    table["strength"] = strength
    best = table.sort_values("strength", ascending=False).drop_duplicates("orf").index
    table["is_best_probe"] = table.index.isin(best)
    table["n_probes"] = table.groupby("orf")["orf"].transform("size")
    return table, missing
