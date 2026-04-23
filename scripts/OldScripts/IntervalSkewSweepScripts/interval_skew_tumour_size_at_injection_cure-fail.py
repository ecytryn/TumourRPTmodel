#!/usr/bin/env python3
"""
Extract tumour size at second injection from IS sweep and plot
overlapping histograms for cure vs failure outcomes.

Usage:
    python is_tumour_size_at_injection.py                        # auto-find latest
    python is_tumour_size_at_injection.py 2026-03-19_10-00-00   # specific timestamp
"""

import numpy as np
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

# Cell type column indices in populations.csv
NORMAL   = 1
HYPOXIC  = 2
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


def extract_tumour_size_at_injection(sweep_dir):
    """Walk all run directories and extract tumour size at second injection."""
    records = []

    run_dirs = sorted(Path(sweep_dir).glob("interval_*_skew_*_rep_*"))
    if not run_dirs:
        print(f"ERROR: No run directories found in {sweep_dir}")
        sys.exit(1)

    print(f"Found {len(run_dirs)} run directories.")

    for run_dir in run_dirs:
        # Parse interval, skew, replicate from directory name
        parts = run_dir.name.split('_')
        try:
            interval = int(parts[1])
            skew     = int(parts[3])
            rep      = int(parts[5])
        except (IndexError, ValueError):
            print(f"  Skipping unrecognised directory: {run_dir.name}")
            continue

        pop_file     = run_dir / "populations.csv"
        summary_file = run_dir / "parameters.csv"  # not used here

        if not pop_file.exists():
            print(f"  Warning: no populations.csv in {run_dir.name}")
            continue

        # Load populations (skip header row)
        try:
            pops = np.loadtxt(pop_file, delimiter=',', skiprows=1)
        except Exception as e:
            print(f"  Warning: could not load {pop_file}: {e}")
            continue

        # Row index for second injection day
        second_injection_day = FIRST_INJECTION_DAY + interval
        row_idx = second_injection_day * 24  # 24 rows per day

        if row_idx >= pops.shape[0]:
            print(f"  Warning: row {row_idx} out of bounds for {run_dir.name} "
                  f"(only {pops.shape[0]} rows)")
            continue

        row = pops[row_idx]
        tumour_size = row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC]

        records.append({
            'interval': interval,
            'skew':     skew,
            'replicate': rep,
            'tumour_size_at_injection2': int(tumour_size),
        })

    df = pd.DataFrame(records)
    print(f"Extracted {len(df)} records.")
    return df


def merge_with_outcomes(df, sweep_dir):
    """Merge tumour size data with cure/failure outcomes from sweep_summary.csv."""
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found")
        sys.exit(1)

    summary = pd.read_csv(csv_path)
    # Normalise skew column: sweep_summary stores skew in nmol as float
    # directory names store skew as integer nmol
    summary['skew_int'] = summary['skew'].astype(int)

    merged = df.merge(
        summary[['interval', 'skew_int', 'replicate', 'outcome']],
        left_on=['interval', 'skew', 'replicate'],
        right_on=['interval', 'skew_int', 'replicate'],
        how='left'
    )

    missing = merged['outcome'].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} records could not be matched to outcomes.")

    return merged


def save_spreadsheet(df, sweep_dir):
    out_path = Path(sweep_dir) / "tumour_size_at_injection2.csv"
    cols = ['interval', 'skew', 'replicate', 'tumour_size_at_injection2', 'outcome']
    df[cols].to_csv(out_path, index=False)
    print(f"Spreadsheet saved to: {out_path}")
    return out_path


def plot_histograms(df, sweep_dir):
    cure_sizes    = df.loc[df['outcome'] == 'CURE',    'tumour_size_at_injection2']
    failure_sizes = df.loc[df['outcome'] == 'FAILURE', 'tumour_size_at_injection2']

    print(f"\nCure:    n={len(cure_sizes)}, "
          f"median={cure_sizes.median():.0f}, mean={cure_sizes.mean():.0f}")
    print(f"Failure: n={len(failure_sizes)}, "
          f"median={failure_sizes.median():.0f}, mean={failure_sizes.mean():.0f}")

    all_sizes = df['tumour_size_at_injection2']
    max_size  = all_sizes.max()
    bins = np.linspace(0, max_size * 1.05, 40)

    fig, ax = plt.subplots()

    ax.hist(failure_sizes, bins=bins, density=True, alpha=0.6,
            color='#E84855', label='Failure', zorder=2)
    ax.hist(cure_sizes, bins=bins, density=True, alpha=0.6,
            color='#2E86AB', label='Cure', zorder=3)

    ax.set_xlabel('Tumour size at 2nd injection (cells)')
    ax.set_ylabel('Probability density')
    ax.legend(framealpha=0.95)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    out_pdf = Path(sweep_dir) / "tumour_size_at_injection2_hist.pdf"
    out_png = Path(sweep_dir) / "tumour_size_at_injection2_hist.png"
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"Histogram saved to: {out_pdf}")

    plt.show()


def main():
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)
    print(f"Using sweep directory: {sweep_dir}\n")

    df = extract_tumour_size_at_injection(sweep_dir)
    df = merge_with_outcomes(df, sweep_dir)
    save_spreadsheet(df, sweep_dir)
    plot_histograms(df, sweep_dir)


if __name__ == "__main__":
    main()
