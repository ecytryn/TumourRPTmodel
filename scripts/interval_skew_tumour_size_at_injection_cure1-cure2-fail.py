#!/usr/bin/env python3
"""
Extract tumour size at second injection from IS sweep and plot
overlapping histograms for cure vs failure outcomes.

Usage:
    python is_tumour_size_at_injection.py                          # auto-find latest, no filter
    python is_tumour_size_at_injection.py --min_interval 10        # exclude intervals < 10
    python is_tumour_size_at_injection.py 2026-03-19_10-00-00      # specific timestamp
    python is_tumour_size_at_injection.py 2026-03-19_10-00-00 --min_interval 14
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

# Cell type column indices in populations.csv
NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4

FIRST_INJECTION_DAY = 5  # must match IntervalSkewSweep


def find_sweep_dir(timestamp=None):
    base = "results/IntervalSkewSweep/IntervalSkewSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            print(f"ERROR: {d} not found")
            sys.exit(1)
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        print("ERROR: No IntervalSkewSweep directories found.")
        sys.exit(1)
    return dirs[-1]


def tumour_size_from_row(row):
    return int(row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC])


def extract_tumour_size_at_injection(sweep_dir):
    """Walk all run directories and extract tumour size at second injection."""
    records = []

    run_dirs = sorted(Path(sweep_dir).glob("interval_*_skew_*_rep_*"))
    if not run_dirs:
        print(f"ERROR: No run directories found in {sweep_dir}")
        sys.exit(1)

    print(f"Found {len(run_dirs)} run directories.")
    n_early_stop = 0

    for run_dir in run_dirs:
        parts = run_dir.name.split("_")
        try:
            interval = int(parts[1])
            skew = int(parts[3])
            rep = int(parts[5])
        except (IndexError, ValueError):
            print(f"  Skipping unrecognised directory: {run_dir.name}")
            continue

        pop_file = run_dir / "populations.csv"
        if not pop_file.exists():
            print(f"  Warning: no populations.csv in {run_dir.name}")
            continue

        try:
            pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        except Exception as e:
            print(f"  Warning: could not load {pop_file}: {e}")
            continue

        second_injection_day = FIRST_INJECTION_DAY + interval
        row_idx2 = second_injection_day * 24
        row_idx1 = FIRST_INJECTION_DAY * 24

        if row_idx2 >= pops.shape[0]:
            # Early stop - tumour eliminated before second injection
            # Record tumour size at first injection day instead
            early_stop = True
            n_early_stop += 1
            tumour_size = tumour_size_from_row(pops[row_idx1])
        else:
            early_stop = False
            tumour_size = tumour_size_from_row(pops[row_idx2])

        records.append(
            {
                "interval": interval,
                "skew": skew,
                "replicate": rep,
                "tumour_size_at_injection2": tumour_size,
                "early_stop": early_stop,
            }
        )

    df = pd.DataFrame(records)
    print(
        f"Extracted {len(df)} records ({n_early_stop} early stops - "
        f"size recorded at 1st injection instead)."
    )
    return df


def merge_with_outcomes(df, sweep_dir):
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found")
        sys.exit(1)

    summary = pd.read_csv(csv_path)
    summary["skew_int"] = summary["skew"].astype(int)

    merged = df.merge(
        summary[["interval", "skew_int", "replicate", "outcome"]],
        left_on=["interval", "skew", "replicate"],
        right_on=["interval", "skew_int", "replicate"],
        how="left",
    )

    missing = merged["outcome"].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} records could not be matched to outcomes.")

    return merged


def save_spreadsheet(df, sweep_dir, min_interval):
    cols = [
        "interval",
        "skew",
        "replicate",
        "tumour_size_at_injection2",
        "early_stop",
        "outcome",
    ]
    suffix = f"_minint{min_interval}" if min_interval > 0 else ""
    out_path = Path(sweep_dir) / f"tumour_size_at_injection2{suffix}.csv"
    df[cols].to_csv(out_path, index=False)
    print(f"Spreadsheet saved to: {out_path}")


def plot_histograms(df, sweep_dir, min_interval):
    early_cure_sizes = df.loc[
        (df["outcome"] == "CURE") & (df["early_stop"] == True),
        "tumour_size_at_injection2",
    ]
    late_cure_sizes = df.loc[
        (df["outcome"] == "CURE") & (df["early_stop"] == False),
        "tumour_size_at_injection2",
    ]
    #    cure_sizes    = df.loc[df['outcome'] == 'CURE',    'tumour_size_at_injection2']
    failure_sizes = df.loc[df["outcome"] == "FAILURE", "tumour_size_at_injection2"]

    n_early_cure = len(early_cure_sizes)
    n_late_cure = len(late_cure_sizes)
    n_failure = len(failure_sizes)
    n_total = n_early_cure + n_late_cure + n_failure

    print(
        f"\nAfter filtering (min_interval={min_interval}): {n_total} runs "
        f"({n_early_cure} early cure, {n_late_cure} late cure, {n_failure} failure)"
    )
    if n_early_cure > 0:
        print(
            f"  Cure:    median={early_cure_sizes.median():.0f}, "
            f"mean={early_cure_sizes.mean():.0f}"
        )
    if n_late_cure > 0:
        print(
            f"  Cure:    median={late_cure_sizes.median():.0f}, "
            f"mean={late_cure_sizes.mean():.0f}"
        )
    if n_failure > 0:
        print(
            f"  Failure: median={failure_sizes.median():.0f}, "
            f"mean={failure_sizes.mean():.0f}"
        )

    if n_early_cure == 0 or n_late_cure == 0 or n_failure == 0:
        print("  Warning: one outcome has 0 runs - histogram will be incomplete.")

    all_sizes = df["tumour_size_at_injection2"]
    bins = np.linspace(0, all_sizes.max() * 1.05, 40)
    bin_width = bins[1] - bins[0]

    fig, ax = plt.subplots()

    if n_early_cure > 0:
        weights = np.ones(n_early_cure) / (n_total * bin_width)
        ax.hist(
            early_cure_sizes,
            bins=bins,
            weights=weights,
            alpha=0.6,
            color="#EE6677",
            label=f"Cure at first inj. (n={n_early_cure})",
            zorder=4,
        )
    #        ax.hist(cure_sizes, bins=bins, density=True, alpha=0.6,
    #                color='#2E86AB', label=f'Cure (n={n_cure})', zorder=3)

    if n_late_cure > 0:
        weights = np.ones(n_late_cure) / (n_total * bin_width)
        ax.hist(
            late_cure_sizes,
            bins=bins,
            weights=weights,
            alpha=0.6,
            color="#AA3377",
            label=f"Cure at second inj. (n={n_late_cure})",
            zorder=3,
        )
    #        ax.hist(cure_sizes, bins=bins, density=True, alpha=0.6,
    #                color='#2E86AB', label=f'Cure (n={n_cure})', zorder=3)

    if n_failure > 0:
        weights = np.ones(n_failure) / (n_total * bin_width)
        ax.hist(
            failure_sizes,
            bins=bins,
            weights=weights,
            alpha=0.6,
            color="#4477AA",
            label=f"Failure (n={n_failure})",
            zorder=2,
        )
    #        ax.hist(failure_sizes, bins=bins, density=True, alpha=0.6,
    #                color='#E84855', label=f'Failure (n={n_failure})', zorder=2)

    filter_str = f", interval ≥ {min_interval} d" if min_interval > 0 else ""
    ax.set_xlabel(f"Tumour size at successful injection (cells){filter_str}")
    ax.set_ylabel("Probability density")
    ax.legend(framealpha=0.95)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    suffix = f"_minint{min_interval}" if min_interval > 0 else ""
    out_pdf = Path(sweep_dir) / f"tumour_size_at_injection2_hist{suffix}.pdf"
    out_png = Path(sweep_dir) / f"tumour_size_at_injection2_hist{suffix}.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Histogram saved to: {out_pdf}")

    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "timestamp", nargs="?", default=None, help="Sweep timestamp (optional)"
    )
    parser.add_argument(
        "--min_interval",
        type=int,
        default=0,
        help="Exclude runs with interval < this value (default: 0 = no filter)",
    )
    args = parser.parse_args()

    sweep_dir = find_sweep_dir(args.timestamp)
    print(f"Using sweep directory: {sweep_dir}")
    print(f"Minimum interval filter: {args.min_interval} days\n")

    df = extract_tumour_size_at_injection(sweep_dir)
    df = merge_with_outcomes(df, sweep_dir)
    # Apply skew filter
    df = df[df["skew"] == 0].copy()

    # Apply interval filter
    if args.min_interval > 0:
        n_before = len(df)
        df = df[df["interval"] >= args.min_interval].copy()
        print(
            f"Filtered from {n_before} to {len(df)} runs "
            f"(excluded interval < {args.min_interval})."
        )

    save_spreadsheet(df, sweep_dir, args.min_interval)
    plot_histograms(df, sweep_dir, args.min_interval)


if __name__ == "__main__":
    main()
