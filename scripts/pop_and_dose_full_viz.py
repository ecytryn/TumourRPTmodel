#!/usr/bin/env python3
"""
Plot tumour population and dose rate over time for a single simulation run.

Usage:
    python pop_and_dose_full_viz.py <run_dir>
    python pop_and_dose_full_viz.py <run_dir> --title "My plot"

Example:
    python pop_and_dose_full_viz.py results/TumourSizeThresholdSweep/TumourSizeThresholdSweep_2026-03-27_02-05-53/radius_50um_rep_1
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update(
    {
        "figure.figsize": (3.35, 2.4),
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.2,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# Cell type indices
NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4


def smooth_daily(arr, window=24):
    """Apply 24-point rolling mean to remove daily update artefacts."""
    if arr.ndim == 1:
        return np.convolve(arr, np.ones(window) / window, mode="same")
    return np.apply_along_axis(
        lambda x: np.convolve(x, np.ones(window) / window, mode="same"), 0, arr
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", help="Path to simulation run directory")
    parser.add_argument(
        "--title", default=None, help="Figure title (default: run directory name)"
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    title = args.title if args.title else run_dir.name

    # Load populations
    pop_file = run_dir / "populations.csv"
    if not pop_file.exists():
        print(f"ERROR: populations.csv not found in {run_dir}")
        sys.exit(1)
    populations = np.loadtxt(pop_file, delimiter=",", skiprows=1)
    print(f"Loaded populations: shape = {populations.shape}")

    # Load dose
    dose_file = run_dir / "dose.csv"
    if not dose_file.exists():
        print(f"ERROR: dose.csv not found in {run_dir}")
        sys.exit(1)
    doses = np.loadtxt(dose_file, delimiter=",", skiprows=1)
    if doses.ndim > 1:
        doses = doses[:, 0]
    print(f"Loaded dose: shape = {doses.shape}")

    # Align lengths
    n = min(len(populations), len(doses))
    populations = populations[:n]
    doses = doses[:n]
    # doses = smooth_daily(doses)

    # Time axis
    t_days = np.arange(n) / 24.0

    # Cell populations
    normoxic = populations[:, NORMAL]
    hypoxic = populations[:, HYPOXIC]
    necrotic = populations[:, NECROTIC]
    apoptotic = populations[:, APOPTOTIC]
    viable = normoxic + hypoxic

    # Plot
    fig, ax1 = plt.subplots()

    ax1.plot(
        t_days, normoxic, color="#C7AD84", linewidth=1.2, label="Normoxic", zorder=4
    )
    ax1.plot(t_days, hypoxic, color="#A76D34", linewidth=1.2, label="Hypoxic", zorder=3)
    ax1.plot(
        t_days,
        apoptotic,
        color="#82AC5D",
        linewidth=1.5,
        #        linestyle="--",
        label="Apoptotic",
        zorder=2,
    )

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Cell count")
    ax1.set_ylim(bottom=0)

    # Secondary axis for dose rate
    ax2 = ax1.twinx()
    ax2.fill_between(t_days, doses, alpha=0.2, color="#4477AA", step="pre")
    ax2.plot(
        t_days,
        doses,
        color="#4477AA",
        linewidth=1.2,
        alpha=0.6,
        drawstyle="steps-pre",
        label="Dose rate",
        zorder=1,
    )
    ax2.set_ylabel("Dose rate (Gy/h)")
    ax2.set_ylim(bottom=0)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", framealpha=0.9)

    #    ax1.set_title(title, fontsize=10, pad=3)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    out_pdf = run_dir / "pop_dose.pdf"
    out_png = run_dir / "pop_dose.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()
