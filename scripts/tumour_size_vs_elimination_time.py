#!/usr/bin/env python3
"""
Plot mean time to tumour elimination (± SD) vs tumour size at injection,
for multiple capillary density sweeps overlaid on the same axes.
Only cured runs are included.

Usage:
    python tumour_size_vs_elimination_time.py TIMESTAMP1 TIMESTAMP2 TIMESTAMP3
    python tumour_size_vs_elimination_time.py --labels "605" "566" "374" TIMESTAMP1 TIMESTAMP2 TIMESTAMP3
"""

import argparse
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

# CB-safe colours for the three cap densities
COLOURS = ["#4477AA", "#228833", "#EE6677", "#CCBB44", "#AA3377", "#66CCEE"]

# populations.csv column indices
NORMAL = 1
HYPOXIC = 2
INJECTION_DAY = 5


def find_sweep_dir(timestamp):
    base = "results/TumourSizeThresholdSweep/TumourSizeThresholdSweep"
    d = f"{base}_{timestamp}"
    if not Path(d).exists():
        print(f"ERROR: {d} not found")
        sys.exit(1)
    return d


SIMULATION_DAYS = 45  # INJECTION_DAY + FOLLOW_UP_DAYS
FULL_RUN_ROWS = SIMULATION_DAYS * 24 + 1  # +1 for the initial zero row


def get_elimination_day(run_dir):
    pop_file = Path(run_dir) / "populations.csv"
    if not pop_file.exists():
        return None
    try:
        pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        if pops.shape[0] < FULL_RUN_ROWS:
            # Early stop — cured; subtract 1 for the zero row, convert to days
            return (pops.shape[0] - 1) / 24.0
        else:
            return None
    except Exception as e:
        print(f"  Warning: could not read {pop_file}: {e}")
        return None


def get_tumour_size_at_injection(run_dir):
    """Return tumour size (cells) at injection day."""
    pop_file = Path(run_dir) / "populations.csv"
    if not pop_file.exists():
        return None
    try:
        pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        row = pops[INJECTION_DAY * 24]
        return int(
            row[NORMAL] + row[HYPOXIC] + row[3] + row[4]
        )  # include necrotic/apoptotic for size
    except Exception:
        return None


def load_sweep(sweep_dir):
    """Load all runs from a sweep directory.
    Returns DataFrame with columns: radius_um, replicate, outcome,
    tumour_size_at_inj, elimination_day."""
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    summary = pd.read_csv(csv_path)
    summary["outcome"] = summary["outcome"].str.strip("'")

    print(f"  Sample row: {summary.iloc[0].to_dict()}")
    print(f"  CURE rows: {(summary['outcome'].str.contains('CURE')).sum()}")

    records = []
    for _, row in summary.iterrows():
        radius_um = row["radius_um"]
        replicate = int(row["replicate"])
        outcome = row["outcome"]

        run_dir = Path(sweep_dir) / f"radius_{int(radius_um)}um_rep_{replicate}"

        size = get_tumour_size_at_injection(run_dir)

        elim = get_elimination_day(run_dir) if "CURE" in outcome else None

        records.append(
            {
                "radius_um": radius_um,
                "replicate": replicate,
                "outcome": outcome,
                "tumour_size": size,
                "elimination_day": elim,
            }
        )

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "timestamps", nargs="+", help="Sweep timestamps (one per cap density)"
    )
    parser.add_argument(
        "--labels", nargs="+", default=None, help="Legend labels e.g. '605' '566' '374'"
    )
    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.timestamps):
        print("ERROR: number of labels must match number of timestamps")
        sys.exit(1)

    labels = args.labels if args.labels else args.timestamps

    fig, ax = plt.subplots()

    all_sizes = []

    for timestamp, label, colour in zip(args.timestamps, labels, COLOURS):
        sweep_dir = find_sweep_dir(timestamp)
        print(f"\nLoading {sweep_dir}...")
        df = load_sweep(sweep_dir)

        # Only cured runs with valid size and elimination day
        cured = df[
            (df["outcome"].str.contains("CURE"))
            & df["tumour_size"].notna()
            & df["elimination_day"].notna()
        ].copy()

        print(f"  {len(cured)} cured runs out of {len(df)} total")

        # Group by radius, compute mean/SD of elimination day
        # Use median tumour size per radius as x value
        grouped = (
            cured.groupby("radius_um")
            .agg(
                mean_elim=("elimination_day", "mean"),
                sd_elim=("elimination_day", "std"),
                n=("elimination_day", "count"),
                median_size=("tumour_size", "median"),
            )
            .reset_index()
            .sort_values("radius_um")
        )

        # Only plot radii with at least 3 cured runs for stable stats
        grouped = grouped[grouped["n"] >= 3]

        print(
            grouped[
                ["radius_um", "median_size", "mean_elim", "sd_elim", "n"]
            ].to_string(index=False)
        )

        all_sizes.extend(grouped["median_size"].tolist())

        ax.errorbar(
            grouped["median_size"],
            grouped["mean_elim"] - INJECTION_DAY,  # days after injection
            yerr=grouped["sd_elim"],
            fmt="o-",
            color=colour,
            capsize=2,
            linewidth=1.2,
            label=f"{label} cap/mm²",
            zorder=3,
        )

    ax.set_xlabel("Tumour size at injection (cells)")
    ax.set_ylabel("Days to elimination\n(after injection)")
    ax.set_ylim(bottom=0)
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    # Save to first sweep dir
    out_dir = find_sweep_dir(args.timestamps[0])
    out_pdf = Path(out_dir) / "elimination_time_vs_size.pdf"
    out_png = Path(out_dir) / "elimination_time_vs_size.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()
