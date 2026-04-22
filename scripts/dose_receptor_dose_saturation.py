#!/usr/bin/env python3
"""
Plot captive RL (N_b^hot + N_ic^hot) vs time from DR sweep simulations
at three receptor densities, across all injected amounts.
Three horizontal panels, shared y-axis scale.

Usage:
    python dose_receptor_dose_saturation.py
    python dose_receptor_dose_saturation.py 2026-03-19_10-00-00
"""

import glob
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update(
    {
        "figure.figsize": (6.7, 2.4),
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "lines.linewidth": 0.8,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# Receptor densities to plot (mol/cell)
# TARGET_RECEPS = [4.4e-19, 5.4e-19, 6.4e-19]
TARGET_RECEPS = [3.4e-19, 4.9e-19, 6.4e-19]
RECEP_TOL = 0.05e-19
RECEP_COLOURS = ["#4477AA", "#228833", "#EE6677"]
# RECEP_PEAK_REFS = [8.48e-8, 1.045e-7, 1.168e-7]
RECEP_PEAK_REFS = [6.8e-8, 9e-8, 1.168e-7]
RECEP_LABELS = [
    "$R_C=3.4\\times10^{-10}$ nmol/cell",
    "$R_C=4.9\\times10^{-10}$ nmol/cell",
    "$R_C=6.4\\times10^{-10}$ nmol/cell",
]

INJECTION_DAY = 5
MIN_PER_DAY = 1440.0


def find_sweep_dir(timestamp=None):
    base = "results/DoseReceptorSweep/DoseReceptorSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            print(f"ERROR: {d} not found")
            sys.exit(1)
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        print("ERROR: No DoseReceptorSweep directories found.")
        sys.exit(1)
    return dirs[-1]


def load_captive_rl(run_dir):
    """Load N_b_hot + N_ic_hot from a single run directory.
    Returns array in nmol, or None on failure."""
    pk_file = Path(run_dir) / "pkStateVariables.csv"
    if not pk_file.exists():
        return None
    try:
        df = pd.read_csv(pk_file)
        captive = (df["N_b_hot"] + df["N_ic_hot"]) * 1e9  # mol -> nmol
        return captive.values
    except Exception as e:
        print(f"  Warning: could not read {pk_file}: {e}")
        return None


def mean_with_nan_padding(arrays):
    """Average a list of arrays of unequal length,
    padding shorter ones with NaN so cured (early-stop) runs
    don't truncate the mean."""
    max_len = max(len(a) for a in arrays)
    padded = np.full((len(arrays), max_len), np.nan)
    for i, a in enumerate(arrays):
        padded[i, : len(a)] = a
    return np.nanmean(padded, axis=0)


def load_runs_for_recep(all_runs, target_recep):
    return [
        (dose, recep, rep, d)
        for dose, recep, rep, d in all_runs
        if abs(recep - target_recep) <= RECEP_TOL
    ]


def main():
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)
    print(f"Sweep dir: {sweep_dir}")

    # Parse all run directories once
    run_dirs = sorted(Path(sweep_dir).glob("dose_*_recep_*_rep_*"))
    if not run_dirs:
        print("ERROR: No run directories found.")
        sys.exit(1)

    pattern = re.compile(r"dose_(\d+\.?\d*)_recep_([\d.e+\-]+)_rep_(\d+)")
    all_runs = []
    for d in run_dirs:
        m = pattern.match(d.name)
        if not m:
            continue
        all_runs.append((float(m.group(1)), float(m.group(2)), int(m.group(3)), d))

    inj_hour = INJECTION_DAY * 24

    fig, axes = plt.subplots(1, 3, sharey=True, sharex=True)

    # Compute global y-max across all three groups for shared scale
    # (sharey handles this automatically, but we print it for reference)

    for ax, target_recep, colour, label in zip(
        axes, TARGET_RECEPS, RECEP_COLOURS, RECEP_LABELS
    ):
        runs = load_runs_for_recep(all_runs, target_recep)
        if not runs:
            print(f"WARNING: No runs found for recep ≈ {target_recep:.2e}")
            continue

        doses = sorted(set(r[0] for r in runs))
        print(
            f"recep={target_recep:.1e}: {len(doses)} dose levels, "
            f"{len(runs)} total runs"
        )

        for dose in doses:
            rep_data = []
            for dv, rv, rep, run_dir in runs:
                if dv == dose:
                    arr = load_captive_rl(run_dir)
                    if arr is not None:
                        rep_data.append(arr)

            if not rep_data:
                continue

            mean_captive = mean_with_nan_padding(rep_data)

            t_hours = np.arange(len(mean_captive)) - inj_hour
            t_days = t_hours / 24.0
            mask = t_days >= 0
            t_plot = t_days[mask]
            cap_plot = mean_captive[mask]

            ax.plot(t_plot, cap_plot, color=colour, linewidth=0.8, alpha=0.85, zorder=3)

        # Draw peak reference lines for all three receptor densities
        # (zorder=1 puts them behind the curves)
        for ref_peak, ref_colour in zip(RECEP_PEAK_REFS, RECEP_COLOURS):
            ax.axhline(
                ref_peak,
                color=ref_colour,
                linewidth=0.8,
                linestyle="--",
                alpha=0.7,
                zorder=1,
            )

        # Panel label as title
        ax.set_title(label, pad=3)
        ax.set_xlabel("Time since injection (days)")
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.set_xlim(left=0, right=40)
        ax.set_ylim(bottom=0, top=1.2e-7)

    # Only leftmost panel gets y-axis label
    axes[0].set_ylabel("$N_\\mathrm{captive}^H$ (nmol)")

    # Remove redundant y tick labels from middle and right panels
    # (sharey keeps the scale identical; just suppress the numbers)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    plt.tight_layout()

    out_pdf = Path(sweep_dir) / "pk_saturation_3recep.pdf"
    out_png = Path(sweep_dir) / "pk_saturation_3recep.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()
