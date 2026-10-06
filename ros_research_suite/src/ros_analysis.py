"""Auswertung der ROS-Gene und datengetriebene Entdeckung (Phase 6).

Zwei Blickrichtungen auf dieselben Daten:
  1. Vorwissen -> Daten:  Wie verhalten sich die BEKANNTEN ROS-Gene?
  2. Daten -> Vorwissen:  Welche Gene reagieren am stärksten, und wie viele davon stehen
                          überhaupt in der kuratierten Liste?
Die Liste ist nur Interpretationshilfe; sie fließt nicht in die Berechnung der Antwort ein.
"""
import pandas as pd

from src.ros_genes import CATEGORY_LABELS, known_orfs


def category_summary(ros_table: pd.DataFrame, treatments: list[str]) -> pd.DataFrame:
    """Pro Kategorie: Anzahl Gene auf dem Chip und responsive Gene je Behandlung.

    Pro Gen zählt nur die beste Sonde (sonst würden Gene mit 2 Sonden doppelt gezählt).
    """
    best = ros_table[ros_table["is_best_probe"]]
    rows = []
    for category, label in CATEGORY_LABELS.items():
        members = best[best["category"] == category]
        row = {"category": label, "n_genes": len(members)}
        for treatment in treatments:
            row[f"{treatment}_up"] = int(members[f"{treatment}_up"].sum())
            row[f"{treatment}_down"] = int(members[f"{treatment}_down"].sum())
            row[f"{treatment}_median_abs_peak"] = round(
                float(members[f"{treatment}_peak_log2fc"].abs().median()), 2) if len(members) else float("nan")
        rows.append(row)
    return pd.DataFrame(rows).set_index("category")


def responsive_share(result: pd.DataFrame, ros_table: pd.DataFrame, treatment: str) -> tuple[float, float]:
    """Anteil responsiver Gene: (in der ROS-Liste, im Rest des Chips)."""
    in_list = result.index.isin(ros_table.index[ros_table["is_best_probe"]])
    column = result[f"{treatment}_responsive"]
    return float(column[in_list].mean()), float(column[~in_list].mean())


def top_responsive(result: pd.DataFrame, treatment: str, n: int,
                   patterns: pd.DataFrame | None = None, orfs: set[str] | None = None) -> pd.DataFrame:
    """Die n responsiven Sonden mit dem größten |log2FC| am stärksten Zeitpunkt.

    Spalte 'in_ros_list' zeigt, ob das Gen in der kuratierten Liste steht.
    """
    orfs = known_orfs() if orfs is None else orfs
    responsive = result[result[f"{treatment}_responsive"]]
    order = responsive[f"{treatment}_peak_log2fc"].abs().sort_values(ascending=False).head(n).index
    top = responsive.loc[order, [c for c in ("display_name", "systematic_name") if c in result.columns]
                         + [f"{treatment}_direction", f"{treatment}_peak_log2fc",
                            f"{treatment}_peak_time_min", f"{treatment}_n_timepoints_above"]].copy()
    if patterns is not None:
        top[f"{treatment}_pattern"] = patterns.loc[order, f"{treatment}_pattern"]
    top["in_ros_list"] = top["systematic_name"].isin(orfs)
    top.insert(0, "rank", range(1, len(top) + 1))
    return top


def discovery_summary(top: pd.DataFrame, result: pd.DataFrame, orfs: set[str] | None = None) -> dict:
    """Vergleich "bekannte ROS-Gene" gegen "unerwartet stark reagierende Gene".

    base_rate = Anteil der Liste an ALLEN Sonden des Chips (so viele Treffer wären
    ohne Zusammenhang zu erwarten). Liegt der Anteil in den Top-Genen deutlich darüber,
    sind die bekannten Gene unter den stärksten Antwortern angereichert.
    """
    orfs = known_orfs() if orfs is None else orfs
    return {
        "n_top": len(top),
        "n_in_list": int(top["in_ros_list"].sum()),
        "share_in_top": float(top["in_ros_list"].mean()) if len(top) else float("nan"),
        "base_rate": float(result["systematic_name"].isin(orfs).mean()),
    }
