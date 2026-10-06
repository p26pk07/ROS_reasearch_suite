"""Wissenschaftliche Grafiken (matplotlib). Dateien landen in results/figures/.

Backend 'Agg' = nur Bilddateien schreiben, kein Fenster (funktioniert überall).
"""
from pathlib import Path

import matplotlib
import matplotlib.ticker  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src import statistics  # noqa: E402
from src.differential_expression import timepoint_columns  # noqa: E402


def ma_plots(expression: pd.DataFrame, design: pd.DataFrame, result: pd.DataFrame,
             treatment: str, threshold: float, path: Path) -> Path:
    """MA-Plots einer Behandlung: je Zeitpunkt ein Panel (gegen t=0).

    M = log2FC (Differenz der log2-Werte), A = mittlere log2-Expression beider Messungen.
    Zeigt, ob starke Änderungen nur bei schwach exprimierten (rauschigen) Genen
    auftreten. Rot = laut Kriterium responsiv.
    """
    genes = result.index
    rows = design[design["treatment"] == treatment]
    later = rows[~rows["is_baseline"].to_numpy(dtype=bool)]
    times = sorted(later["time_min"].unique())
    responsive = result[f"{treatment}_responsive"].to_numpy(dtype=bool)

    fig, axes = plt.subplots(1, len(times), figsize=(3.2 * len(times), 3.6),
                             sharey=True, squeeze=False)
    for axis, time_min in zip(axes[0], times):
        ma = statistics.ma_values(expression.loc[genes], design, treatment, time_min)
        m_values, a_values = ma["M"].to_numpy(), ma["A"].to_numpy()
        axis.scatter(a_values[~responsive], m_values[~responsive], s=3, c="lightgrey",
                     rasterized=True)
        axis.scatter(a_values[responsive], m_values[responsive], s=6, c="tab:red",
                     rasterized=True)
        for level in (-threshold, 0, threshold):
            axis.axhline(level, color="black", lw=0.7, ls="-" if level == 0 else "--")
        axis.set_title(f"{time_min:g} min")
        axis.set_xlabel("A (mittlere log2-Expression)")
    axes[0][0].set_ylabel("M (log2FC vs. t=0)")
    fig.suptitle(f"MA-Plots {treatment} (rot = responsiv)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def response_scatter(result: pd.DataFrame, treatment: str, threshold: float,
                     path: Path, label_top: int = 8) -> Path:
    """Ersatz für den Volcano-Plot (der braucht p-Werte): Stärke gegen Konsistenz.

    x = log2FC am stärksten Zeitpunkt, y = Anzahl Zeitpunkte mit |log2FC| >= Schwelle.
    Gene oben links/rechts reagieren stark UND über viele Zeitpunkte.
    """
    x = result[f"{treatment}_peak_log2fc"].to_numpy()
    y = result[f"{treatment}_n_timepoints_above"].to_numpy()
    responsive = result[f"{treatment}_responsive"].to_numpy(dtype=bool)
    jitter = np.random.default_rng(0).uniform(-0.15, 0.15, size=len(y))  # nur gegen Überlappung

    fig, axis = plt.subplots(figsize=(8, 5))
    axis.scatter(x[~responsive], y[~responsive] + jitter[~responsive], s=4, c="lightgrey",
                 rasterized=True)
    axis.scatter(x[responsive], y[responsive] + jitter[responsive], s=8, c="tab:red",
                 rasterized=True)
    for level in (-threshold, threshold):
        axis.axvline(level, color="black", lw=0.7, ls="--")

    # Beschriftung: viele Top-Gene haben fast gleiche Koordinaten (gleiche Zeitpunkt-Anzahl),
    # direkt am Punkt würden sich die Namen überlagern. Deshalb stehen die Namen in je einer
    # Spalte rechts (hochreguliert) bzw. links (herunterreguliert) und sind per Linie verbunden.
    names = result["display_name"] if "display_name" in result.columns else result.index.to_series()
    peak_column = f"{treatment}_peak_log2fc"
    count_column = f"{treatment}_n_timepoints_above"
    peaks = result.loc[responsive, peak_column]
    for sign, x_position, alignment in ((1, 1.03, "left"), (-1, -0.08, "right")):
        chosen = list(peaks[peaks * sign > 0].abs().nlargest(label_top).index)
        # nach Zeitpunkt-Anzahl und x sortiert -> Linien kreuzen sich möglichst wenig
        chosen.sort(key=lambda probe: (result.loc[probe, count_column],
                                       sign * result.loc[probe, peak_column]))
        for probe, slot in zip(chosen, np.linspace(0.08, 0.92, len(chosen))):
            axis.annotate(
                str(names[probe]),
                xy=(result.loc[probe, peak_column], result.loc[probe, count_column]),
                xytext=(x_position, slot), textcoords="axes fraction",
                ha=alignment, va="center", fontsize=7, annotation_clip=False,
                arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "grey"})
    axis.set_xlabel("log2FC am stärksten Zeitpunkt")
    # Kein y-Achsentitel: er würde mit den Namen links kollidieren. Die Erklärung steht im Titel.
    axis.set_title(f"Antwortstärke vs. Konsistenz – {treatment}\n"
                   "y-Achse: Anzahl Zeitpunkte mit |log2FC| ≥ Schwelle", fontsize=10)
    axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    fig.savefig(path, dpi=150, bbox_inches="tight")  # tight: Namen außerhalb der Achsen bleiben im Bild
    plt.close(fig)
    return path


def _trajectories(table: pd.DataFrame, treatment: str) -> tuple[np.ndarray, pd.DataFrame]:
    """Zeiten (inkl. t=0) und log2FC-Verläufe (Gene x Zeiten) einer Behandlung; t=0 ist 0."""
    columns = timepoint_columns(table, treatment)
    times = np.array([0.0] + [time for time, _ in columns])
    values = table[[name for _, name in columns]].copy()
    values.insert(0, "t0", 0.0)
    values.columns = times
    return times, values


def gene_timecourse_grid(table: pd.DataFrame, genes: list[str], treatments: list[str],
                         threshold: float, path: Path, columns: int = 4) -> Path:
    """Ein Panel pro Gen: log2FC über die Zeit, eine Linie je Behandlung.

    Hinweis: Es gibt keine unbehandelte Kontrolle; die Linien zeigen die zwei
    Behandlungen (jeweils gegen ihr eigenes t=0).
    """
    rows = int(np.ceil(len(genes) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.4 * columns, 2.6 * rows),
                             sharex=True, squeeze=False)
    names = table["display_name"] if "display_name" in table.columns else table.index.to_series()
    colors = plt.get_cmap("tab10")
    trajectories = {t: _trajectories(table, t) for t in treatments}
    for axis, gene in zip(axes.ravel(), genes):
        for index, treatment in enumerate(treatments):
            times, values = trajectories[treatment]
            axis.plot(times, values.loc[gene].to_numpy(), marker="o", ms=3, lw=1.4,
                      color=colors(index), label=treatment)
        for level in (-threshold, threshold):
            axis.axhline(level, color="black", lw=0.6, ls="--")
        axis.axhline(0, color="grey", lw=0.5)
        axis.set_title(str(names[gene]), fontsize=9)
    for axis in axes.ravel()[len(genes):]:
        axis.axis("off")
    for axis in axes[-1]:
        axis.set_xlabel("Zeit (min)")
    for axis in axes[:, 0]:
        axis.set_ylabel("log2FC vs. t=0")
    axes[0][0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def pattern_profiles(table: pd.DataFrame, patterns: pd.Series, treatment: str, threshold: float,
                     path: Path, min_genes: int = 5, max_panels: int = 8) -> Path:
    """Ein Panel pro Zeitverlaufs-Muster: alle Gene (grau) und der Median (rot).

    Zeigt, ob die Muster-Einteilung wirklich zusammengehörige Verläufe trennt.
    """
    counts = patterns[patterns != "none"].value_counts()
    chosen = [name for name, n in counts.items() if n >= min_genes][:max_panels]
    if not chosen:
        chosen = list(counts.index[:1])
    times, values = _trajectories(table, treatment)
    columns = min(4, len(chosen))
    rows = int(np.ceil(len(chosen) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.6 * columns, 3 * rows),
                             sharex=True, squeeze=False)
    rng = np.random.default_rng(0)
    for axis, name in zip(axes.ravel(), chosen):
        members = values.loc[patterns.index[patterns == name]]
        shown = members if len(members) <= 150 else members.iloc[rng.choice(len(members), 150, replace=False)]
        axis.plot(times, shown.to_numpy().T, color="lightgrey", lw=0.6)
        axis.plot(times, members.median().to_numpy(), color="tab:red", lw=2)
        for level in (-threshold, threshold):
            axis.axhline(level, color="black", lw=0.6, ls="--")
        axis.set_title(f"{name}\n(n = {len(members)})", fontsize=9)
    for axis in axes.ravel()[len(chosen):]:
        axis.axis("off")
    for axis in axes[-1]:
        axis.set_xlabel("Zeit (min)")
    for axis in axes[:, 0]:
        axis.set_ylabel("log2FC vs. t=0")
    fig.suptitle(f"Zeitverlaufs-Muster {treatment} (rot = Median)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def ros_gene_barplot(ros_table: pd.DataFrame, treatments: list[str], category_labels: dict,
                     threshold: float, path: Path) -> Path:
    """Balkendiagramm der kuratierten ROS-Gene: log2FC am stärksten Zeitpunkt je Behandlung.

    Gene sind nach Kategorie gruppiert und innerhalb der ersten Behandlung nach Stärke sortiert.
    Kräftige Balken = laut Kriterium responsiv, blasse Balken = nicht responsiv.
    """
    best = ros_table[ros_table["is_best_probe"]]
    lead = treatments[0]
    ordered = []
    for category in category_labels:
        members = best[best["category"] == category]
        ordered.append(members.sort_values(f"{lead}_peak_log2fc", ascending=False))
    data = pd.concat(ordered)
    n = len(data)
    height = 0.8 / len(treatments)
    colors = plt.get_cmap("tab10")

    fig, axis = plt.subplots(figsize=(8, max(5, 0.2 * n + 1.5)))
    positions = np.arange(n)
    for index, treatment in enumerate(treatments):
        values = data[f"{treatment}_peak_log2fc"].to_numpy()
        responsive = data[f"{treatment}_responsive"].to_numpy(dtype=bool)
        bars = axis.barh(positions + (index - (len(treatments) - 1) / 2) * height, values, height=height,
                         color=colors(index), label=treatment)
        for bar, is_responsive in zip(bars, responsive):
            bar.set_alpha(1.0 if is_responsive else 0.3)
    axis.set_yticks(positions)
    axis.set_yticklabels(data["symbol"], fontsize=6)
    axis.set_ylim(n - 0.5, -0.5)   # ohne Leerraum oben/unten, erstes Gen oben
    for level in (-threshold, 0, threshold):
        axis.axvline(level, color="black", lw=0.6, ls="-" if level == 0 else "--")

    start = 0
    for category, label in category_labels.items():
        count = int((data["category"] == category).sum())
        if count:
            axis.axhline(start - 0.5, color="grey", lw=0.5)
            axis.text(1.01, start + count / 2 - 0.5, label.split(" (")[0], transform=axis.get_yaxis_transform(),
                      fontsize=7, va="center")
            start += count
    axis.set_xlabel("log2FC am stärksten Zeitpunkt (vs. t=0)")
    axis.set_title("Kuratierte ROS-/Stress-Gene (kräftig = responsiv)")
    axis.legend(fontsize=7, loc="lower right")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path
