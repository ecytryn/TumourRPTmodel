#!/usr/bin/env python3
"""
Side-by-side comparison of cure rate vs tumour size at injection:
  Left panel:  IS sweep (skew=0) — binned cure rate with Wilson CI
  Right panel: TumourSizeThresholdSweep — cure rate vs median size at injection

Usage:
    python cure_rate_vs_size_comparison.py
    python cure_rate_vs_size_comparison.py --is_timestamp 2026-03-19_10-00-00
    python cure_rate_vs_size_comparison.py --is_timestamp TS1 --tst_timestamp TS2
    python cure_rate_vs_size_comparison.py --min_interval 10
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

# populations.csv column indices
NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4

FIRST_INJECTION_DAY = 5  # must match sweep Java code


# =============================================================================
# Shared utilities
# =============================================================================


def wilson_ci(n_cure, n_total, z=1.96):
    """Wilson score confidence interval for a proportion."""
    if n_total == 0:
        return np.nan, np.nan
    p = n_cure / n_total
    denom = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denom
    half = z * np.sqrt(p * (1 - p) / n_total + z**2 / (4 * n_total**2)) / denom
    return centre - half, centre + half


def tumour_size_from_row(row):
    return int(row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC])


# =============================================================================
# IS sweep data loading (reused from cure1-cure2-fail script)
# =============================================================================


def find_is_dir(timestamp=None):
    base = "results/IntervalSkewSweep/IntervalSkewSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            print(f"ERROR: IS sweep directory not found: {d}")
            sys.exit(1)
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        print("ERROR: No IntervalSkewSweep directories found.")
        sys.exit(1)
    return dirs[-1]


def load_is_data(sweep_dir, min_interval):
    """Load IS sweep, filter to skew=0 and min_interval, return df with
    tumour_size_at_last_injection and outcome columns."""
    run_dirs = sorted(Path(sweep_dir).glob("interval_*_skew_*_rep_*"))
    if not run_dirs:
        print(f"ERROR: No run directories found in {sweep_dir}")
        sys.exit(1)

    records = []
    for run_dir in run_dirs:
        parts = run_dir.name.split("_")
        try:
            interval = int(parts[1])
            skew = int(parts[3])
            rep = int(parts[5])
        except (IndexError, ValueError):
            continue

        if skew not in (-15, -10, -5, 0, 5, 10, 15):
            continue
        if interval < min_interval:
            continue

        pop_file = run_dir / "populations.csv"
        if not pop_file.exists():
            continue

        try:
            pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        except Exception:
            continue

        row_idx1 = FIRST_INJECTION_DAY * 24
        row_idx2 = (FIRST_INJECTION_DAY + interval) * 24

        if row_idx2 >= pops.shape[0]:
            # Early stop — use size at first injection
            tumour_size = tumour_size_from_row(pops[row_idx1])
        else:
            tumour_size = tumour_size_from_row(pops[row_idx2])

        records.append(
            {
                "interval": interval,
                "replicate": rep,
                "tumour_size_at_last_inj": tumour_size,
            }
        )

    df = pd.DataFrame(records)

    # Merge outcomes from sweep_summary.csv
    summary = pd.read_csv(Path(sweep_dir) / "sweep_summary.csv")
    summary = summary[summary["skew"].astype(int) == 0]
    df = df.merge(
        summary[["interval", "replicate", "outcome"]],
        on=["interval", "replicate"],
        how="left",
    )

    missing = df["outcome"].isna().sum()
    if missing > 0:
        print(f"  Warning: {missing} IS records missing outcome.")

    df = df[df["outcome"].isin(["CURE", "FAILURE"])].copy()
    print(
        f"IS sweep: {len(df)} runs (skew=0, interval≥{min_interval}), "
        f"{(df['outcome'] == 'CURE').sum()} cures, "
        f"{(df['outcome'] == 'FAILURE').sum()} failures."
    )
    return df


def bin_cure_rate(df, size_col, n_bins=20):
    """Bin cure rate by tumour size, return df with bin centres, cure rate, CI."""
    sizes = df[size_col]
    bins = np.linspace(0, sizes.max() * 1.05, n_bins + 1)
    bin_centres = 0.5 * (bins[:-1] + bins[1:])

    cure_rates, lo, hi, ns = [], [], [], []
    for i in range(len(bins) - 1):
        mask = (sizes >= bins[i]) & (sizes < bins[i + 1])
        subset = df[mask]
        n = len(subset)
        n_cure = (subset["outcome"] == "CURE").sum()
        if n == 0:
            cure_rates.append(np.nan)
            lo.append(np.nan)
            hi.append(np.nan)
        else:
            p = n_cure / n
            l, h = wilson_ci(n_cure, n)
            cure_rates.append(p)
            lo.append(l)
            hi.append(h)
        ns.append(n)

    result = pd.DataFrame(
        {
            "bin_centre": bin_centres,
            "cure_rate": cure_rates,
            "ci_lo": lo,
            "ci_hi": hi,
            "n": ns,
        }
    )
    return result[result["n"] > 0].reset_index(drop=True)


# =============================================================================
# Threshold sweep data loading (reused from tumour_cellcount_vs_cure_plot.py)
# =============================================================================


def find_tst_dir(timestamp=None):
    base = "results/TumourSizeThresholdSweep/TumourSizeThresholdSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            print(f"ERROR: TST sweep directory not found: {d}")
            sys.exit(1)
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        print("ERROR: No TumourSizeThresholdSweep directories found.")
        sys.exit(1)
    return dirs[-1]


def load_tst_data(sweep_dir):
    """Load threshold sweep, return df grouped by radius with cure rate and
    median tumour size at injection."""
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    df = pd.read_csv(csv_path)
    df = df.dropna(subset=["radius_um", "replicate"])

    # Read tumour size at injection from populations.csv for each run
    sizes = []
    for _, row in df.iterrows():
        run_dir = (
            Path(sweep_dir)
            / f"radius_{int(row['radius_um'])}um_rep_{int(row['replicate'])}"
        )
        pop_file = run_dir / "populations.csv"
        size = None
        if pop_file.exists():
            try:
                pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
                row_idx = FIRST_INJECTION_DAY * 24
                if row_idx < pops.shape[0]:
                    size = tumour_size_from_row(pops[row_idx])
            except Exception:
                pass
        sizes.append(size)

    df["tumour_size_at_inj"] = sizes
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)

    grouped = (
        df.groupby("radius_um")
        .agg(
            cure_rate=("is_cure", "mean"),
            n=("is_cure", "count"),
            n_cure=("is_cure", "sum"),
            median_size=("tumour_size_at_inj", "median"),
        )
        .reset_index()
        .sort_values("radius_um")
    )

    # Wilson CI per radius group
    grouped["ci_lo"] = grouped.apply(
        lambda r: wilson_ci(r["n_cure"], r["n"])[0], axis=1
    )
    grouped["ci_hi"] = grouped.apply(
        lambda r: wilson_ci(r["n_cure"], r["n"])[1], axis=1
    )

    print(f"Threshold sweep: {len(df)} runs across {len(grouped)} radii.")
    return grouped


# =============================================================================
# Plotting
# =============================================================================


def make_comparison_plot(is_binned, tst_grouped, min_interval, out_dir):

    COLOUR_IS = "#EE6677"  # rose
    COLOUR_TST = "#4477AA"  # blue

    x_min = max(is_binned["bin_centre"].min(), tst_grouped["median_size"].min()) * 0.92
    x_max = min(is_binned["bin_centre"].max(), tst_grouped["median_size"].max()) * 1.08

    fig, ax = plt.subplots()

    # IS sweep
    valid = is_binned.dropna(subset=["cure_rate"])
    ax.errorbar(
        valid["bin_centre"],
        valid["cure_rate"],
        yerr=[
            np.maximum(0, valid["cure_rate"] - valid["ci_lo"]),
            np.maximum(0, valid["ci_hi"] - valid["cure_rate"]),
        ],
        fmt="o-",
        color=COLOUR_IS,
        capsize=2,
        zorder=3,
        label="Interval-skew sweep (|S|<15)",
    )

    # Threshold sweep
    ax.errorbar(
        tst_grouped["median_size"],
        tst_grouped["cure_rate"],
        yerr=[
            np.maximum(0, tst_grouped["cure_rate"] - tst_grouped["ci_lo"]),
            np.maximum(0, tst_grouped["ci_hi"] - tst_grouped["cure_rate"]),
        ],
        fmt="s-",
        color=COLOUR_TST,
        capsize=2,
        zorder=3,
        label="Tumour size sweep",
    )

    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xlim(20, 170)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Tumour size at last injection (cells)")
    ax.set_ylabel("Cure rate")
    filter_str = f", interval ≥ {min_interval} d" if min_interval > 0 else ""
    ax.set_title(f"Cure rate vs tumour size at injection{filter_str}", pad=3)
    ax.legend(framealpha=0.95)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    suffix = f"_minint{min_interval}" if min_interval > 0 else ""
    out_pdf = Path(out_dir) / f"cure_rate_vs_size_comparison{suffix}.pdf"
    out_png = Path(out_dir) / f"cure_rate_vs_size_comparison{suffix}.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--is_timestamp",
        default=None,
        help="IntervalSkewSweep timestamp (default: latest)",
    )
    parser.add_argument(
        "--tst_timestamp",
        default=None,
        help="TumourSizeThresholdSweep timestamp (default: latest)",
    )
    parser.add_argument(
        "--min_interval",
        type=int,
        default=0,
        help="Exclude IS runs with interval < this value (default: 0)",
    )
    args = parser.parse_args()

    is_dir = find_is_dir(args.is_timestamp)
    tst_dir = find_tst_dir(args.tst_timestamp)
    print(f"IS sweep dir:  {is_dir}")
    print(f"TST sweep dir: {tst_dir}\n")

    is_df = load_is_data(is_dir, args.min_interval)
    is_binned = bin_cure_rate(is_df, "tumour_size_at_last_inj", n_bins=20)
    tst_grouped = load_tst_data(tst_dir)

    make_comparison_plot(is_binned, tst_grouped, args.min_interval, is_dir)


if __name__ == "__main__":
    main()
