#!/usr/bin/env python3
"""
Visualizes TumourSizeThresholdSweep results.

Plots cure rate vs initial tumour radius, with sample size annotated.

Usage:
    python tumour_size_threshold_plot.py                        # Auto-find latest
    python tumour_size_threshold_plot.py 2026-03-12_12-37-05   # Specific timestamp
"""

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


def find_sweep_dir(timestamp=None):
    base = "results/TumourSizeThresholdSweep/TumourSizeThresholdSweep"
    if timestamp:
        sweep_dir = f"{base}_{timestamp}"
        if not Path(sweep_dir).exists():
            print(f"ERROR: Directory not found: {sweep_dir}")
            sys.exit(1)
        return sweep_dir
    dirs = glob.glob(f"{base}_*")
    if not dirs:
        print("ERROR: No TumourSizeThresholdSweep directories found.")
        sys.exit(1)
    return sorted(dirs)[-1]


NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4
INJECTION_DAY = 5


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


def main():
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)
    csv_path = f"{sweep_dir}/sweep_summary.csv"

    if not Path(csv_path).exists():
        print(f"ERROR: {csv_path} not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    df = df.dropna(subset=["radius_um", "replicate"])

    print(f"Loaded {len(df)} rows from {csv_path}")

    # Extract tumour size at injection day for each run
    sizes = []
    for _, row in df.iterrows():
        s = tumour_size_at_injection(sweep_dir, row["radius_um"], row["replicate"])
        sizes.append(s)
    df["tumour_size_at_injection"] = sizes

    # Aggregate cure rate and median tumour size per radius
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)
    grouped = (
        df.groupby("radius_um")
        .agg(
            cure_rate=("is_cure", "mean"),
            n=("is_cure", "count"),
            cell_count=("cell_count", "first"),
            median_size_at_inj=("tumour_size_at_injection", "median"),
        )
        .reset_index()
    )
    grouped = grouped.sort_values("radius_um")

    print("\nSummary:")
    for _, row in grouped.iterrows():
        print(
            f"  r={row['radius_um']:.0f} µm: initial={row['cell_count']:.0f} cells, "
            f"at injection={row['median_size_at_inj']:.0f} cells, "
            f"cure rate={row['cure_rate']:.2f} ({row['cure_rate'] * row['n']:.0f}/{row['n']:.0f})"
        )

    fig, ax = plt.subplots()

    ax.plot(
        grouped["median_size_at_inj"],
        grouped["cure_rate"],
        "o-",
        color="steelblue",
        zorder=3,
    )

    # Annotate each point with initial radius
    #  for _, row in grouped.iterrows():
    #      ax.annotate(
    #          f"r={row['radius_um']:.0f} µm",
    #          xy=(row["median_size_at_inj"], row["cure_rate"]),
    #          xytext=(4, 4),
    #          textcoords="offset points",
    #          fontsize=6,
    #          color="dimgray",
    #      )

    ax.set_xlabel("Tumour size at injection (cells)")
    ax.set_ylabel("Cure Rate")
    ax.set_xlim(10, 125)
    ax.set_ylim(-0.05, 1.05)
    #    ax.set_xticks(grouped["median_size_at_inj"])
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", zorder=1)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    out_path = f"{sweep_dir}/cure_rate_vs_size_at_injection.pdf"
    plt.savefig(out_path, format="pdf", bbox_inches="tight")
    print(f"\nSaved to: {out_path}")


if __name__ == "__main__":
    main()
