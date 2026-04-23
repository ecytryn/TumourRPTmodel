#!/usr/bin/env python3
"""
Visualizes TumourSizeThresholdSweep results.

Plots cure rate vs initial tumour radius, with sample size annotated.

Usage:
    python tumour_size_threshold_plot.py                        # Auto-find latest
    python tumour_size_threshold_plot.py 2026-03-12_12-37-05   # Specific timestamp
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import sys
from pathlib import Path

mpl.rcParams.update({
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
    "ps.fonttype": 42
})


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


def main():
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)
    csv_path = f"{sweep_dir}/sweep_summary.csv"

    if not Path(csv_path).exists():
        print(f"ERROR: {csv_path} not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows from {csv_path}")

    # Aggregate cure rate per radius
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    grouped = df.groupby('radius_um').agg(
        cure_rate=('is_cure', 'mean'),
        n=('is_cure', 'count'),
        cell_count=('cell_count', 'first')
    ).reset_index()
    grouped = grouped.sort_values('radius_um')

    print("\nSummary:")
    for _, row in grouped.iterrows():
        print(f"  r={row['radius_um']:.0f} µm ({row['cell_count']:.0f} cells): "
              f"{row['cure_rate']*row['n']:.0f}/{row['n']:.0f} cures "
              f"= {row['cure_rate']:.2f}")

    fig, ax = plt.subplots()

    ax.plot(grouped['radius_um'], grouped['cure_rate'],
            'o-', color='steelblue', zorder=3)

    # Annotate each point with cell count
#    for _, row in grouped.iterrows():
#        ax.annotate(f"{row['cell_count']:.0f} cells",
#                    xy=(row['radius_um'], row['cure_rate']),
#                    xytext=(4, 4), textcoords='offset points',
#                    fontsize=6, color='dimgray')

    ax.set_xlabel('Initial Tumour Radius (µm)')
    ax.set_ylabel('Cure Rate')
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticks(grouped['radius_um'])
    ax.axhline(0.5, color='gray', linewidth=0.8, linestyle='--', zorder=1)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    out_path = f"{sweep_dir}/cure_rate_vs_radius.pdf"
    plt.savefig(out_path, format='pdf', bbox_inches='tight')
    print(f"\nSaved to: {out_path}")


if __name__ == "__main__":
    main()
