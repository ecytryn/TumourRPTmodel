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
        "figure.figsize": (3.35, 2.4),
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
TARGET_RECEPS = [4.0e-19, 6.0e-19, 8.0e-19]
RECEP_TOL = 0.05e-19
# RECEP_PEAK_REFS = [7.35e-8, 1.12e-7, 1.375e-7]
RECEP_LABELS = [
    "$R_C=4.0\\times10^{-10}$ nmol/cell",
    "$R_C=6.0\\times10^{-10}$ nmol/cell",
    "$R_C=8.0\\times10^{-10}$ nmol/cell",
]
RECEP_COLOURS = ["#4477AA", "#228833", "#EE6677"]

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


def smooth_daily(arr, window=24):
    """Apply 24-point rolling mean to remove daily update artefacts."""
    if arr.ndim == 1:
        return np.convolve(arr, np.ones(window) / window, mode="same")
    return np.apply_along_axis(
        lambda x: np.convolve(x, np.ones(window) / window, mode="same"), 0, arr
    )


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

    fig, ax = plt.subplots()

    colour_max = {colour: 0.0 for colour in RECEP_COLOURS}

    for target_recep, colour, label in zip(
        reversed(TARGET_RECEPS), reversed(RECEP_COLOURS), reversed(RECEP_LABELS)
    ):
        runs = load_runs_for_recep(all_runs, target_recep)
        if not runs:
            print(f"WARNING: No runs found for recep ≈ {target_recep:.2e}")
            continue

        doses = sorted(set(r[0] for r in runs))

        first_dose = True
        for dose in doses:
            # if dose % 100 and dose % 200:
            if dose < 50:
                continue
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
            # cap_plot = smooth_daily(cap_plot)

            colour_max[colour] = max(colour_max[colour], np.nanmax(cap_plot))

            ax.plot(
                t_plot,
                cap_plot,
                color=colour,
                linewidth=0.8,
                alpha=0.85,
                zorder=3,
                label=label if first_dose else None,
            )
            first_dose = False

    # Auto-computed reference lines from data maxima
    print("\nPer-receptor captive RL peak values (nmol):")
    for ref_colour, label in zip(RECEP_COLOURS, RECEP_LABELS):
        peak = colour_max[ref_colour]
        print(f"  {label}: {peak:.4e} nmol")
        ax.axhline(
            peak,
            color=ref_colour,
            linewidth=0.8,
            linestyle="--",
            alpha=0.7,
            zorder=1,
        )

    ax.set_xlabel("Time since injection (days)")
    ax.set_ylabel("$N_\\mathrm{captive}^H$ ($\\times 10^{-7}$ nmol)")
    ax.set_xlim(left=0, right=40)
    #    ax.set_ylim(bottom=0, top=1.4e-7)
    global_max = max(colour_max.values())
    ax.set_ylim(bottom=0, top=1.05 * global_max)

    ax.yaxis.get_offset_text().set_visible(False)
    ax.yaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, _: f"{x * 1e7:.1f}")
    )
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    # Only leftmost panel gets y-axis label
    ax.set_ylabel("$N_\\mathrm{captive}^H$ ($\\times 10^{-7}$ nmol)")
    ax.yaxis.get_offset_text().set_visible(False)

    plt.tight_layout()

    out_pdf = Path(sweep_dir) / "pk_saturation_3recep.pdf"
    out_png = Path(sweep_dir) / "pk_saturation_3recep.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()
