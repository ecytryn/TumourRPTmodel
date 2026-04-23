#!/usr/bin/env python3
"""
Visualizes dose-receptor sweep results as heatmap

Usage:
    python visualize_dose_receptor.py                    # Auto-find latest
    python visualize_dose_receptor.py 2026-01-27_15-30-00  # Specific timestamp
"""

import glob
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Match interval_skew formatting settings
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


def find_sweep_dir(timestamp=None):
    """Find sweep directory"""
    if timestamp:
        sweep_dir = f"results/DoseReceptorSweep/DoseReceptorSweep_{timestamp}"
        if not Path(sweep_dir).exists():
            print(f"ERROR: Directory not found: {sweep_dir}")
            sys.exit(1)
        return sweep_dir

    # Find most recent
    sweep_dirs = glob.glob("results/DoseReceptorSweep/DoseReceptorSweep_*")
    if not sweep_dirs:
        print("ERROR: No DoseReceptorSweep directories found!")
        sys.exit(1)

    return sorted(sweep_dirs)[-1]


def load_data(sweep_dir):
    """Load sweep results"""
    csv_path = f"{sweep_dir}/sweep_summary.csv"

    if not Path(csv_path).exists():
        print(f"ERROR: Results file not found: {csv_path}")
        sys.exit(1)

    df = pd.read_csv(csv_path)

    print(f"Loaded {len(df)} simulation results from {sweep_dir}")
    print("\nParameter ranges:")
    print(f"  Injected amounts: {sorted(df['dose_nmol'].unique())} nmol")
    print(f"  Receptor densities: {sorted(df['receptors_per_cell_mol'].unique())}")
    print("\nOutcome summary:")
    print(df["outcome"].value_counts())

    return df, sweep_dir


def create_heatmap(df, output_path):
    """Create cure rate heatmap: injected amount vs receptor density"""

    # Calculate cure rate per parameter combination
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)
    cure_rates = (
        df.groupby(["dose_nmol", "receptors_per_cell_mol"])["is_cure"]
        .mean()
        .reset_index()
    )
    cure_rates.columns = ["dose_nmol", "receptors_per_cell_mol", "cure_rate"]

    # Pivot for heatmap
    pivot = cure_rates.pivot(
        index="receptors_per_cell_mol", columns="dose_nmol", values="cure_rate"
    )

    # Sort (high receptors at top)
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot[sorted(pivot.columns)]

    # Create figure
    fig, ax = plt.subplots()

    # Custom colormap: red (failure) → yellow → green (cure)
    #    cmap = plt.cm.RdYlGn
    cmap = "cividis"

    # Plot heatmap
    sns.heatmap(
        pivot,
        annot=False,
        cmap=cmap,
        cbar_kws={"label": "Cure Rate"},
        linewidths=0.0,
        linecolor="white",
        vmin=0,
        vmax=1,
        ax=ax,
    )

    # Labels - updated terminology
    ax.set_xlabel("Injected Amount (nmol)")
    ax.set_ylabel("Receptor Density, $R_C$ (nmol/cell)")

    # Format y-axis labels with mantissa only, exponent as offset

    # Get the actual y-tick values from the pivot index
    ytick_values = pivot.index.values

    # Determine common exponent (assuming all values have similar order of magnitude)
    if len(ytick_values) > 0:
        common_exponent = int(np.floor(np.log10(np.abs(ytick_values[0]))))
        display_exponent = -10  # display in units of 1e-10 nmol/cell

        # Create labels with just the mantissa (scaled by the common exponent)
        # Show only wanted ticks to reduce crowding
        wanted_positions_y = []
        wanted_labels_y = []
        for i, val in enumerate(ytick_values):
            if i % 5 == 0:
                wanted_positions_y.append(i)
                wanted_labels_y.append(f"{val / 10**common_exponent:.1f}")

        ax.set_yticks(wanted_positions_y)
        ax.set_yticklabels(wanted_labels_y, rotation=0)

        # Add exponent as offset text at top of y-axis
        ax.text(
            0.12,
            1.05,
            f"$\\times 10^{{{display_exponent}}}$",
            transform=ax.transAxes,
            fontsize=7,
            ha="right",
        )

    xtick_values = pivot.columns.values

    wanted_positions = []
    wanted_labels = []
    for i, val in enumerate(xtick_values):
        if i % 4 == 0:
            wanted_positions.append(i)
            wanted_labels.append(f"{val:.0f}")

    ax.set_xticks(wanted_positions)
    ax.set_xticklabels(wanted_labels, rotation=0)

    if "replicate" in df.columns:
        num_reps = df["replicate"].max()
        title = "RPT Treatment Outcome"
    else:
        title = "RPT Treatment Outcome"

    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    print(f"\nHeatmap saved to: {output_path}")

    return fig


def create_dose_response_curves(df, output_path):
    """Create dose-response curves for each receptor density"""

    fig, ax = plt.subplots(figsize=(10, 7))

    # Calculate cure rates
    df["is_cure"] = (df["outcome"] == "CURE").astype(int)

    receptor_levels = sorted(df["receptors_per_cell_mol"].unique())

    for recep in receptor_levels:
        subset = df[df["receptors_per_cell_mol"] == recep]
        cure_by_dose = subset.groupby("dose_nmol")["is_cure"].mean()

        label = f"{recep:.1e} mol/cell"
        ax.plot(
            cure_by_dose.index,
            cure_by_dose.values,
            marker="o",
            linewidth=2,
            markersize=8,
            label=label,
        )

    # Updated terminology
    ax.set_xlabel("Total Injected Amount (nmol)", fontsize=12)
    ax.set_ylabel("Cure Rate", fontsize=12)
    ax.set_title("Dose-Response Curves by Receptor Density", fontsize=14)
    ax.legend(title="Receptor Density", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    print(f"Dose-response curves saved to: {output_path}")

    return fig


def print_summary_stats(df):
    """Print useful summary statistics"""
    print("\n=== SUMMARY STATISTICS ===")

    # Overall cure rate
    overall_cure = (df["outcome"] == "CURE").mean() * 100
    print(f"Overall cure rate: {overall_cure:.1f}%")

    # Cure rate by injected amount
    print("\nCure rate by injected amount:")
    dose_cure = df.groupby("dose_nmol")["outcome"].apply(
        lambda x: (x == "CURE").mean() * 100
    )
    for dose, rate in dose_cure.items():
        print(f"  {dose:.0f} nmol: {rate:.1f}%")

    # Cure rate by receptor density
    print("\nCure rate by receptor density:")
    recep_cure = df.groupby("receptors_per_cell_mol")["outcome"].apply(
        lambda x: (x == "CURE").mean() * 100
    )
    for recep, rate in recep_cure.items():
        print(f"  {recep:.1e} mol/cell: {rate:.1f}%")

    # Find approximate thresholds
    print("\n=== THRESHOLD ANALYSIS ===")

    # For baseline receptor density (1.0e-15), what injected amount gives ~50% cure?
    baseline_recep = 1.0e-15
    baseline_data = df[df["receptors_per_cell_mol"] == baseline_recep]
    if len(baseline_data) > 0:
        baseline_cure = baseline_data.groupby("dose_nmol")["outcome"].apply(
            lambda x: (x == "CURE").mean()
        )
        print(f"\nAt baseline receptor density ({baseline_recep:.1e} mol/cell):")
        for dose, rate in baseline_cure.items():
            print(f"  {dose:.0f} nmol → {rate * 100:.0f}% cure rate")

    # For a fixed injected amount (say 100 nmol), what receptor density gives ~50% cure?
    fixed_dose = 100.0
    dose_data = df[df["dose_nmol"] == fixed_dose]
    if len(dose_data) > 0:
        dose_cure = dose_data.groupby("receptors_per_cell_mol")["outcome"].apply(
            lambda x: (x == "CURE").mean()
        )
        print(f"\nAt {fixed_dose:.0f} nmol injected amount:")
        for recep, rate in dose_cure.items():
            print(f"  {recep:.1e} mol/cell → {rate * 100:.0f}% cure rate")


def main():
    """Main execution"""

    # Get sweep directory
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)

    # Load data
    df, sweep_dir = load_data(sweep_dir)

    # Create visualizations
    print("\nGenerating visualizations...")

    # Main heatmap
    heatmap_path = f"{sweep_dir}/cure_rate_dose_receptor.pdf"
    create_heatmap(df, heatmap_path)

    # Dose-response curves
    curves_path = f"{sweep_dir}/dose_response_curves.pdf"
    create_dose_response_curves(df, curves_path)

    # Print statistics
    print_summary_stats(df)

    print("\n=== VISUALIZATION COMPLETE ===")
    print(f"All outputs in: {sweep_dir}")


if __name__ == "__main__":
    main()
