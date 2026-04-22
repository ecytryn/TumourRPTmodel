#!/usr/bin/env python3
"""
Visualizes TumourSizeThresholdSweep results.
Plots cure rate vs tumour size at injection, with multiple sweep
directories overlaid on the same axes (one curve per capillary density).

Usage:
    python tumour_cellcount_vs_cure_plot.py                          # Auto-find latest
    python tumour_cellcount_vs_cure_plot.py 2026-03-12_12-37-05     # Single timestamp
    python tumour_cellcount_vs_cure_plot.py TS1 TS2 TS3 --labels "374" "470" "605"
"""

import argparse
import glob
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
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# COLOURS = ["#4477AA", "#228833", "#EE6677", "#CCBB44", "#AA3377", "#66CCEE"]
COLOURS = ["#CCBB44", "#4477AA"]

NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4
INJECTION_DAY = 5


def wilson_ci(n_cure, n_total, z=1.96):
    """Wilson score 95% confidence interval for a proportion."""
    if n_total == 0:
        return np.nan, np.nan
    p = n_cure / n_total
    denom = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denom
    half = z * np.sqrt(p * (1 - p) / n_total + z**2 / (4 * n_total**2)) / denom
    return centre - half, centre + half


def find_sweep_dirs(timestamps):
    base = "results/TumourSizeThresholdSweep/TumourSizeThresholdSweep"
    if not timestamps:
        dirs = sorted(glob.glob(f"{base}_*"))
        if not dirs:
            print("ERROR: No TumourSizeThresholdSweep directories found.")
            sys.exit(1)
        return [dirs[-1]]
    result = []
    for ts in timestamps:
        d = f"{base}_{ts}"
        if not Path(d).exists():
            print(f"ERROR: Directory not found: {d}")
            sys.exit(1)
        result.append(d)
    return result


def tumour_size_at_injection(sweep_dir, radius_um, replicate):
    """Read tumour size at injection day from populations.csv."""
    run_dir = Path(sweep_dir) / f"radius_{int(radius_um)}um_rep_{int(replicate)}"
    pop_file = run_dir / "populations.csv"
    if not pop_file.exists():
        return None
    try:
        pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        row_idx = INJECTION_DAY * 24
        if row_idx >= pops.shape[0]:
            row = pops[-1]
        else:
            row = pops[row_idx]
        return int(row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC])
    except Exception as e:
        print(f"  Warning: could not read {pop_file}: {e}")
        return None


def load_sweep(sweep_dir):
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    if not Path(csv_path).exists():
        print(f"ERROR: {csv_path} not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    df["outcome"] = df["outcome"].astype(str).str.strip().str.strip("'")
    df = df.dropna(subset=["radius_um", "replicate"])

    sizes = []
    for _, row in df.iterrows():
        s = tumour_size_at_injection(sweep_dir, row["radius_um"], row["replicate"])
        sizes.append(s)
    df["tumour_size_at_injection"] = sizes

    df["is_cure"] = df["outcome"].str.contains("CURE").astype(int)

    grouped = (
        df.groupby("radius_um")
        .agg(
            cure_rate=("is_cure", "mean"),
            n=("is_cure", "count"),
            n_cure=("is_cure", "sum"),
            median_size=("tumour_size_at_injection", "median"),
        )
        .reset_index()
        .sort_values("radius_um")
    )
    grouped["ci_lo"] = grouped.apply(
        lambda r: wilson_ci(r["n_cure"], r["n"])[0], axis=1
    )
    grouped["ci_hi"] = grouped.apply(
        lambda r: wilson_ci(r["n_cure"], r["n"])[1], axis=1
    )
    return grouped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "timestamps",
        nargs="*",
        default=None,
        help="Sweep timestamps (one per capillary density)",
    )
    parser.add_argument(
        "--labels", nargs="+", default=None, help="Legend labels e.g. '374' '470' '605'"
    )
    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.timestamps):
        print("ERROR: number of labels must match number of timestamps")
        sys.exit(1)

    sweep_dirs = find_sweep_dirs(args.timestamps)
    labels = (
        args.labels
        if args.labels
        else [Path(d).name.split("_", 2)[-1] for d in sweep_dirs]
    )

    fig, ax = plt.subplots()

    for sweep_dir, label, colour in zip(sweep_dirs, labels, COLOURS):
        print(f"Loading {sweep_dir}...")
        grouped = load_sweep(sweep_dir)

        print(f"  {len(grouped)} radius groups")

        ax.errorbar(
            grouped["median_size"],
            grouped["cure_rate"],
            yerr=[
                np.maximum(0, grouped["cure_rate"] - grouped["ci_lo"]),
                np.maximum(0, grouped["ci_hi"] - grouped["cure_rate"]),
            ],
            fmt="o-",
            color=colour,
            capsize=2,
            zorder=3,
            label=f"{label} cap/mm²",
        )

    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xlabel("Tumour size at injection (cells)")
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 50, 100, 200, 500, 1000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_ylabel("Cure Rate")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.legend(framealpha=0.9)

    plt.tight_layout()

    out_dir = sweep_dirs[0]
    out_pdf = Path(out_dir) / "cure_rate_vs_size_multidensity.pdf"
    out_png = Path(out_dir) / "cure_rate_vs_size_multidensity.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()
