"""Standard plots from the tables made by collect_results.py.

Writes PNG files to analysis/plots/. The underlying numbers for every plot are in
results/recipes.csv (means and standard deviations) and results/curves.csv.

Usage (inside the thiolene conda environment, after collect_results.py):
    python analysis/make_plots.py
    python analysis/make_plots.py --results my_results --out my_plots
"""
import argparse
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent

# Reference palette (light mode): categorical slots 1 and 2, ink and chrome.
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

# The two series compared everywhere: homopolymerization off / on.
SERIES = [
    (0.0, BLUE, "thiol-ene only (homopolymerization 0)"),
    (0.15, ORANGE, "with homopolymerization (0.15)"),
]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "sans-serif", "font.size": 10,
    "text.color": INK, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 2, "lines.markersize": 6,
})


def main_sweep(recipes):
    """The group's main study: stiffness 100, reaction check every 100 steps."""
    return recipes[(recipes.angle_constant == 100) & (recipes.polymerize_period == 100)]


def plot_series(ax, data, column, label=True):
    """Mean +/- standard deviation of `column` against crosslinker %, one line per series."""
    for prob, color, name in SERIES:
        d = data[data.homopolymerization_prob == prob].dropna(subset=[f"{column}_mean"])
        ax.errorbar(d.crosslinker_percent, d[f"{column}_mean"], yerr=d[f"{column}_std"],
                    color=color, marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5,
                    elinewidth=1, capsize=3, label=name if label else None)
    ax.set_xlabel("Thiol on tetrathiol (%)")
    ax.set_xticks([0, 25, 50, 75, 100])


def finish(fig, path, title, note):
    """Title on top, wrapped footnote at the bottom, plots in between."""
    width_chars = int(fig.get_figwidth() * 15)
    note = textwrap.fill(note, width_chars)
    note_lines = note.count("\n") + 1
    bottom = (0.25 + 0.17 * note_lines) / fig.get_figheight()
    fig.tight_layout(rect=(0, bottom, 1, 1 - 0.45 / fig.get_figheight()))
    fig.suptitle(title, x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.01, 0.01, note, ha="left", va="bottom", fontsize=8.5, color=MUTED)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def gel_point(recipes, out):
    data = main_sweep(recipes)
    # 0 % tetrathiol without homopolymerization only makes linear chains: no gel, so no gel point
    data = data[~((data.homopolymerization_prob == 0) & (data.crosslinker_percent == 0))]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    x = np.linspace(5, 100, 100)
    ax.plot(x, 1 / np.sqrt(1 + 2 * x / 100), color=MUTED, linewidth=1.25, label="ideal theory, no loops (Flory–Stockmayer)")
    plot_series(ax, data, "gel_point_ene_conversion")
    ax.set_ylabel("C=C conversion at the gel point")
    ax.set_ylim(0.3, 1.0)
    ax.legend(loc="upper right")
    finish(fig, out / "gel_point.png", "Homopolymerization makes the network gel at much lower conversion",
           "Mean ± standard deviation of 5 replicas, 10,000 ene monomers, stiffness 100. Gel point = peak of the reduced weight-average "
           "molecule size, converted to C=C conversion. 0 % without homopolymerization forms no gel and is left out.")


def conversion(recipes, out):
    data = main_sweep(recipes)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    plot_series(axes[0], data, "final_ene_conversion")
    axes[0].set_title("C=C groups")
    plot_series(axes[1], data, "final_thiol_conversion", label=False)
    axes[1].set_title("SH groups")
    axes[0].set_ylabel("Conversion at the end")
    axes[0].set_ylim(0.4, 1.0)
    axes[0].legend(loc="lower left")
    finish(fig, out / "final_conversion.png", "With homopolymerization about half the thiol is left unreacted",
           "Mean ± standard deviation of 5 replicas. Every C=C that homopolymerizes uses no thiol, so SH conversion falls to about 0.5 "
           "while C=C conversion stays near 0.95.")


def network(recipes, out):
    data = main_sweep(recipes)
    panels = [
        ("gel_fraction", "Gel fraction", "Share of beads in the gel"),
        ("elastically_effective_fraction", "Load-carrying share of the gel", "Elastically effective fraction"),
        ("n_primary_loops", "Primary loops", "Number of loops"),
        ("effective_strand_density", "Load-carrying strand density", "Strands per unit volume"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.6))
    for i, (ax, (column, title, ylabel)) in enumerate(zip(axes.flat, panels)):
        plot_series(ax, data, column, label=(i == 0))
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
    axes[0, 0].legend(loc="lower right")
    finish(fig, out / "network_structure.png", "More tetrathiol gives a denser, more load-bearing network",
           "Mean ± standard deviation of 5 replicas, final frame. Primary loop = a strand that starts and ends on the same junction. "
           "Load-carrying (elastically effective) parts are found after shrinking all bonds (contract_bonds).")


def modulus(recipes, jobs, out):
    data = main_sweep(recipes)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    plot_series(axes[0], data, "youngs_modulus")
    axes[0].set_title("Against tetrathiol content")
    axes[0].set_ylabel("Young's modulus (simulation units)")
    axes[0].legend(loc="upper left")

    ax = axes[1]
    main = jobs[(jobs.angle_constant == 100) & (jobs.polymerize_period == 100)].dropna(subset=["youngs_modulus", "effective_strand_density"])
    for prob, color, name in SERIES:
        d = main[main.homopolymerization_prob == prob]
        ax.scatter(d.effective_strand_density, d.youngs_modulus, s=36, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    nu = np.linspace(0, main.effective_strand_density.max() * 1.05, 50)
    ax.plot(nu, 3 * nu, color=MUTED, linewidth=1.25)
    label_x = nu[-1] * 0.55
    ax.annotate("ideal rubber: E = 3νkT", xy=(label_x, 3 * label_x), xytext=(-10, 0), textcoords="offset points",
                ha="right", va="center", fontsize=8.5, color=INK2)
    ax.set_title("Against load-carrying strand density")
    ax.set_xlabel("Load-carrying strands per unit volume (ν)")
    ax.set_ylabel("Young's modulus (simulation units)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=min(0, main.youngs_modulus.min() * 1.1))
    finish(fig, out / "modulus.png", "Stiffness follows the load-carrying strand density",
           "Left: mean ± standard deviation of 5 replicas. Right: one dot per job. Modulus of the isolated gel stretched with bond and "
           "angle forces only, from the first 10 % strain.")


def curves(jobs, curves_table, out):
    """Conversion over time and molecule growth, 50 % tetrathiol, all replicas."""
    pick = jobs[(jobs.angle_constant == 100) & (jobs.polymerize_period == 100) & (jobs.crosslinker_percent == 50)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    for prob, color, name in SERIES:
        for i, job in enumerate(pick[pick.homopolymerization_prob == prob].job_id):
            c = curves_table[curves_table.job_id == job].dropna(subset=["ene_conversion"])
            label = name if i == 0 else None
            axes[0].plot(c.timestep / 1e6, c.ene_conversion, color=color, linewidth=1.25, alpha=0.8, label=label)
            axes[1].plot(c.ene_conversion, c.weight_average_molecule, color=color, linewidth=1.25, alpha=0.8)
    axes[0].set_title("C=C conversion over time")
    axes[0].set_xlabel("Timesteps (millions)")
    axes[0].set_ylabel("C=C conversion")
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc="lower right")
    axes[1].set_title("Molecule growth")
    axes[1].set_xlabel("C=C conversion")
    axes[1].set_ylabel("Weight-average molecule size (beads)")
    axes[1].set_yscale("log")
    axes[1].set_xlim(0, 1)
    finish(fig, out / "conversion_curves.png", "With homopolymerization the molecules grow to a gel much sooner",
           "50 % tetrathiol, stiffness 100, one line per replica (5 each). The weight-average size jumps when the gel spans the box.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=REPO / "analysis" / "results", help="folder with the csv files")
    parser.add_argument("--out", type=Path, default=REPO / "analysis" / "plots", help="output folder")
    args = parser.parse_args()

    recipes = pd.read_csv(args.results / "recipes.csv")
    jobs = pd.read_csv(args.results / "jobs.csv")
    curves_table = pd.read_csv(args.results / "curves.csv")
    args.out.mkdir(parents=True, exist_ok=True)

    gel_point(recipes, args.out)
    conversion(recipes, args.out)
    network(recipes, args.out)
    modulus(recipes, jobs, args.out)
    curves(jobs, curves_table, args.out)


if __name__ == "__main__":
    main()
