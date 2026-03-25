#!/usr/bin/env python3
"""
Visualizes parameter sweep results as:
1. Clean heatmap with cure rate (greenish-yellow to green)
2. Small multiples grid showing population trajectories for each parameter combination

Usage:
    python visualize_interval_skew.py              # Uses default sweep_summary.csv
    python visualize_interval_skew.py _v2          # Uses sweep_summary_v2.csv
"""

import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams.update(
    {
        "figure.figsize": (3.35, 2.4),  # two-column
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,  # editable text in Illustrator
        "ps.fonttype": 42,
    }
)


# Usage: python visualize_interval_skew.py 2026-01-27_01-53-14
TIMESTAMP = sys.argv[1] if len(sys.argv) > 1 else None

if TIMESTAMP:
    CSV_FILE = (
        f"results/IntervalSkewSweep/IntervalSkewSweep_{TIMESTAMP}/sweep_summary.csv"
    )
    SWEEP_DIR = f"results/IntervalSkewSweep/IntervalSkewSweep_{TIMESTAMP}"
else:
    # Try to find most recent
    import glob

    sweep_dirs = glob.glob("results/IntervalSkewSweep/IntervalSkewSweep_*")
    if sweep_dirs:
        SWEEP_DIR = sorted(sweep_dirs)[-1]
        CSV_FILE = f"{SWEEP_DIR}/sweep_summary.csv"
    else:
        print(
            "ERROR: No sweep found. Pass timestamp: python interval_skew_visualize.py 2026-01-27_01-53-14"
        )
        sys.exit(1)

OUTPUT_DIR = SWEEP_DIR  # Put visualizations in same directory as sweep


def load_data(csv_path):
    """Load parameter sweep results"""
    df = pd.read_csv(csv_path)

    # Filter out invalid runs
    df = df[df["outcome"] != "INVALID"]

    print(f"Loaded {len(df)} valid simulation results")
    if "replicate" in df.columns:
        num_replicates = df["replicate"].max()
        print(f"  Replicates per parameter set: {num_replicates}")

    print("Parameter ranges:")
    print(f"  Intervals: {sorted(df['interval'].unique())}")
    print(f"  Skews: {sorted(df['skew'].unique())}")
    print("\nOutcome summary:")
    print(df["outcome"].value_counts())

    return df


def create_custom_colormap():
    """Create greenish-yellow to green colormap"""
    # Yellow-green (#CCFF00ish) to Green (#00AA00ish)
    colors = [
        (0.9, 0.2, 0.2),  # Red for 0 (failure)
        (1.0, 0.9, 0.4),  # Yellow for ~0.3
        (0.7, 0.9, 0.3),  # Yellow-green for ~0.5
        (0.4, 0.8, 0.3),  # Light green for ~0.7
        (0.0, 0.6, 0.2),  # Green for 1.0 (success)
    ]
    positions = [0.0, 0.3, 0.5, 0.7, 1.0]

    cmap = LinearSegmentedColormap.from_list("cure_rate", list(zip(positions, colors)))
    return cmap


def create_cure_rate_heatmap(df, output_path):
    """Create heatmap with cure rate coloring"""

    # Calculate cure rate per parameter combination
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)
    cure_rates = df.groupby(["interval", "skew"])["is_cure"].mean().reset_index()
    cure_rates.columns = ["interval", "skew", "cure_rate"]

    # Pivot for heatmap
    pivot_cure = cure_rates.pivot(index="skew", columns="interval", values="cure_rate")

    # Sort indices
    pivot_cure = pivot_cure.sort_index(ascending=False)
    pivot_cure = pivot_cure[sorted(pivot_cure.columns)]

    # Create figure
    fig, ax = plt.subplots()

    # Custom colormap
    cmap = create_custom_colormap()

    # Plot heatmap WITHOUT annotations
    sns.heatmap(
        pivot_cure,
        annot=False,  # Changed from annot=annot_data
        #                cmap=cmap,
        cmap="cividis",
        cbar_kws={"label": "Cure Rate"},
        linewidths=0.0,
        linecolor="white",
        vmin=0,
        vmax=1,
        ax=ax,
    )

    # Labels and title
    ax.set_xlabel("Inter-injection Interval (days)")
    ax.set_ylabel("Injection Skew (nmol)")

    # Y-axis: show only -20, -10, 0, 10, 20
    skew_vals = pivot_cure.index.tolist()  # already sorted descending
    ax.set_yticks(range(len(skew_vals)))
    ytick_labels = [
        str(int(v)) if int(v) in (-20, -10, 0, 10, 20) else "" for v in skew_vals
    ]
    ax.set_yticklabels(ytick_labels, rotation=0)

    # X-axis: show only 5, 10, 15, 20 (or whatever subset you want)
    interval_vals = pivot_cure.columns.tolist()
    ax.set_xticks(range(len(interval_vals)))
    xtick_labels = [
        str(int(v)) if int(v) in (5, 10, 15, 20) else "" for v in interval_vals
    ]
    ax.set_xticklabels(xtick_labels, rotation=0)

    title = "RPT Treatment Outcome"
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    print(f"\nCure rate heatmap saved to: {output_path}")

    return fig


def load_population_data(sweep_dir, interval, skew, num_replicates=5):
    """Load population trajectories from individual run folders"""
    trajectories = []

    for rep in range(1, num_replicates + 1):
        folder = f"{sweep_dir}/interval_{interval}_skew_{int(skew)}_rep_{rep}"
        pop_file = f"{folder}/populations.csv"

        if os.path.exists(pop_file):
            try:
                # populations.csv has columns for different cell types
                # We want total tumor = normoxic + hypoxic (columns 1 and 2, 0-indexed)
                data = np.loadtxt(pop_file, delimiter=",", skiprows=1)
                if len(data.shape) == 1:
                    data = data.reshape(1, -1)

                # Sum normoxic (index 1) and hypoxic (index 2) tumor cells
                if data.shape[1] >= 3:
                    tumor_pop = data[:, 1] + data[:, 2]
                    trajectories.append(tumor_pop)
            except Exception as e:
                print(f"Warning: Could not load {pop_file}: {e}")

    return trajectories


def create_small_multiples(df, sweep_dir, output_path, max_time=80, max_pop=3000):
    """Create small multiples grid showing population trajectories"""

    intervals = sorted(df["interval"].unique())
    skews = sorted(df["skew"].unique(), reverse=True)  # High skew at top

    n_rows = len(skews)
    n_cols = len(intervals)

    # Calculate cure rates for background coloring
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)
    cure_rates = df.groupby(["interval", "skew"])["is_cure"].mean()

    # Create figure - match main heatmap figure size
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.7, 4.8))

    cmap = create_custom_colormap()

    num_reps = df["replicate"].max() if "replicate" in df.columns else 1

    for i, skew in enumerate(skews):
        for j, interval in enumerate(intervals):
            ax = axes[i, j]

            # Get cure rate for background color
            try:
                cure_rate = cure_rates[(interval, skew)]
            except KeyError:
                cure_rate = 0.5

            # Set background color based on cure rate
            bg_color = cmap(cure_rate)
            ax.set_facecolor(bg_color)  # Darker version

            # Load and plot trajectories
            trajectories = load_population_data(sweep_dir, interval, skew, num_reps)

            if trajectories:
                for traj in trajectories:
                    # Time axis: each entry is 1 hour, plot in days
                    time_days = np.arange(len(traj)) / 24.0
                    ax.plot(time_days, traj, color="darkblue", alpha=0.5, linewidth=0.8)

            # Add dashed vertical lines at injection times
            # First injection at day 0, subsequent at 'interval' days
            # Total of 2 injections (pattern [75+s, 75-s])
            injection_times = [5, 5 + interval]
            for inj_time in injection_times:
                if inj_time <= max_time:
                    ax.axvline(
                        x=inj_time,
                        color="white",
                        linestyle="--",
                        linewidth=0.8,
                        alpha=0.7,
                    )

            # Minimal formatting
            ax.set_xlim(0, max_time)
            ax.set_ylim(0, max_pop)
            ax.set_xticks([])
            ax.set_yticks([])

            # Add interval label on bottom row - match main heatmap font size
            if i == n_rows - 1:
                ax.set_xlabel(f"{interval}", fontsize=8)

            # Add skew label on left column
            if j == 0:
                ax.set_ylabel(f"{int(skew)}", fontsize=8)

            # Add thin border
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_color("gray")

    # Add overall labels
    fig.text(0.5, 0.02, "Inter-injection Interval (days)", ha="center")
    fig.text(0.02, 0.5, "Injection skew (nmol)", va="center", rotation="vertical")

    # Title - single line only
    fig.suptitle("Tumor Population Trajectories", fontsize=8, y=0.98)

    plt.tight_layout()
    plt.subplots_adjust(left=0.08, bottom=0.08, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    print(f"Small multiples saved to: {output_path}")

    return fig


def main():
    """Main execution"""

    # Check if results file exists
    if not Path(CSV_FILE).exists():
        print(f"ERROR: Results file not found: {CSV_FILE}")
        print("Please run IntervalSkewSweep.java first.")
        return

    # Load data
    df = load_data(CSV_FILE)

    # Create output directory
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    # Generate visualizations
    print("\nGenerating visualizations...")

    # Figure 1: Clean heatmap with cure rate (no annotations)
    create_cure_rate_heatmap(df, f"{OUTPUT_DIR}/cure_rate_interval_skew_grid.pdf")

    # Figure 2: Small multiples with population trajectories
    if os.path.exists(SWEEP_DIR):
        create_small_multiples(
            df,
            SWEEP_DIR,
            f"{OUTPUT_DIR}/cure_rate_interval_skew_grid_with_tiny_plots.pdf",
        )
    else:
        print(f"Warning: Sweep directory not found: {SWEEP_DIR}")
        print("Skipping small multiples figure.")

    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    overall_cure_rate = (df["outcome"] == "CURE").mean() * 100
    print(f"Overall cure rate: {overall_cure_rate:.1f}%")

    print("\nCure rate by interval:")
    interval_cure = df.groupby("interval")["outcome"].apply(
        lambda x: (x == "CURE").mean() * 100
    )
    print(interval_cure.to_string())

    print("\nCure rate by skew:")
    skew_cure = df.groupby("skew")["outcome"].apply(
        lambda x: (x == "CURE").mean() * 100
    )
    print(skew_cure.to_string())

    print("\n=== VISUALIZATION COMPLETE ===")


if __name__ == "__main__":
    main()
