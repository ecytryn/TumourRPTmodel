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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

# =============================================================================
# PK PARAMETERS (must match SimParams.java)
# =============================================================================

# K_M in nmol (manuscript units, consistent with dose axis)
k_on_manuscript = 1.5e-3  # L/(nmol·min)
k_off_manuscript = 1.2e-2  # /min
k_int_manuscript = 1e-3  # /min
V_cen_manuscript = 0.5  # L

K_M_nmol = V_cen_manuscript * (k_off_manuscript + k_int_manuscript) / k_on_manuscript
print(f"K_M = {K_M_nmol:.4f} nmol")

K_REL = 2e-4 / 60.0  # 1/s
LAMBDA_DECAY = 7.14e-5 / 60.0  # 1/s
V_CENTRAL = 0.5e-3  # m^3


def fit_boundary(pivot):
    """
    For each dose column, interpolate to find the R_C where cure_rate = 0.5.
    Then fit the single-parameter mechanistic curve:
        R_C = beta * (1 + K_M / A)
    using least squares.

    Returns beta (mol/cell) and the arrays used for the fit.
    """
    dose_vals = np.array(sorted(pivot.columns))  # nmol
    recep_vals = np.array(pivot.index)  # mol/cell (sorted descending)

    # pivot rows are sorted descending; for interpolation we need ascending
    recep_asc = recep_vals[::-1]  # ascending mol/cell

    half_recep = []  # interpolated R_C at cure_rate = 0.5 for each column
    valid_doses = []

    def logistic(x, R_half, k):
        return 1.0 / (1.0 + np.exp(-k * (x - R_half)))

    for dose in dose_vals:
        col = pivot[dose].values[::-1]  # cure rates, ascending with recep_asc

        # Only attempt fit if there's actual variation in this column
        if col.max() <= 0.0 or col.min() >= 1.0:
            continue
        if col.max() - col.min() < 0.1:
            continue

        try:
            # Initial guess: R_half at midpoint of receptor range, k positive
            R_mid = recep_asc[len(recep_asc) // 2]
            R_range = recep_asc.max() - recep_asc.min()
            k_init = (
                4.0 / R_range
            )  # logistic goes 0→1 over roughly 4/k, so this spans the data
            p0 = [R_mid, k_init]
            bounds = ([recep_asc.min(), 0], [recep_asc.max(), np.inf])
            mask = np.isfinite(col)
            if mask.sum() < 3:  # not enough points to fit
                return None, None, None
            recep_asc_clean = recep_asc[mask]
            col_clean = col[mask]
            popt, _ = curve_fit(
                logistic, recep_asc_clean, col_clean, p0=p0, bounds=bounds, maxfev=5000
            )
            # popt, _ = curve_fit(
            #    logistic, recep_asc, col, p0=p0, bounds=bounds, maxfev=5000
            # )

            R_half_fit = popt[0]

            # Only accept if the fitted R_half is within the receptor range
            if recep_asc.min() < R_half_fit < recep_asc.max():
                half_recep.append(R_half_fit)
                valid_doses.append(dose)
        except RuntimeError:
            # curve_fit failed to converge
            pass

    if len(valid_doses) < 2:
        print(
            "Warning: fewer than 2 columns cross cure_rate=0.5; boundary fit unreliable."
        )
        return None, None, None

    A = np.array(valid_doses)  # nmol
    RC = np.array(half_recep)  # mol/cell

    # Closed-form least-squares for R_C = beta * (1 + K_M/A)
    x = 1.0 + K_M_nmol / A  # predictor
    beta = np.dot(RC, x) / np.dot(x, x)

    print("\nBoundary fit:")
    print(f"  Columns crossing 0.5: {len(valid_doses)} / {len(dose_vals)}")
    print(f"  beta = {beta:.4e} mol/cell")
    print(f"  Horizontal asymptote R_C = {beta:.4e} mol/cell")

    return beta, A, RC


def overlay_boundary(ax, pivot, beta, A_fit, RC_fit, cmap):
    """
    Overlay R_C = beta*(1 + K_M/A) on the seaborn heatmap axes,
    plus dots at the (A, R_C) half-cure points used for the fit.
    """
    dose_cols = np.array(sorted(pivot.columns))  # nmol, ascending
    recep_rows = np.array(pivot.index)  # mol/cell, descending

    # Dense dose grid for smooth curve
    A_dense = np.linspace(dose_cols.min(), dose_cols.max(), 500)
    RC_curve = beta * (1.0 + K_M_nmol / A_dense)

    # Map dose → x pixel
    x_pixel = np.interp(A_dense, dose_cols, np.arange(len(dose_cols)) + 0.5)

    # Map R_C → y pixel (recep_rows is descending)
    recep_asc = recep_rows[::-1]
    row_idx_asc = np.arange(len(recep_rows))[::-1]

    y_pixel = np.interp(
        RC_curve, recep_asc, row_idx_asc + 0.5, left=np.nan, right=np.nan
    )

    mask = (RC_curve >= recep_rows.min()) & (RC_curve <= recep_rows.max())

    # White outline then coloured dashed curve
    ax.plot(
        x_pixel[mask],
        y_pixel[mask],
        color=cmap(0.8),
        linewidth=2.5,
        linestyle="-",
        zorder=4,
    )
    ax.plot(
        x_pixel[mask],
        y_pixel[mask],
        color=cmap(0.2),
        linewidth=1,
        linestyle="-",
        label=r"$R_C = \sigma(1 + K_M/A)$",
        zorder=5,
    )

    # --- Dots at the half-cure (A, R_C) points used for the fit ---
    x_dots = np.interp(A_fit, dose_cols, np.arange(len(dose_cols)) + 0.5)
    y_dots = np.interp(RC_fit, recep_asc, row_idx_asc + 0.5, left=np.nan, right=np.nan)

    ax.scatter(
        x_dots, y_dots, color=cmap(0.8), s=5, zorder=6, linewidths=1.5
    )  # outline
    ax.scatter(
        x_dots,
        y_dots,
        color=cmap(0.2),
        s=5,
        zorder=7,
        linewidths=0,
        label="50% cure level",
    )  # fill


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

    # Warp cividis to emphasise the 50% transition
    x = np.linspace(0, 1, 256)
    steepness = 7  # increase for sharper transition, try 6-10
    warped = 1 / (1 + np.exp(-steepness * (x - 0.5)))
    warped = (warped - warped.min()) / (warped.max() - warped.min())
    new_colors = plt.cm.cividis(warped)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "cividis_warped", new_colors, N=256
    )

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

    # Fit and overlay mechanistic boundary
    beta, A_fit, RC_fit = fit_boundary(pivot)
    if beta is not None:
        overlay_boundary(ax, pivot, beta, A_fit, RC_fit, cmap)

    # Labels - updated terminology
    ax.set_xlabel("Injected Amount, A (nmol)")
    ax.set_ylabel("Receptor Density, $R_C$ (nmol/cell)")

    # Format y-axis labels with mantissa only, exponent as offset

    # Desired y-axis ticks in mol/cell (actual stored units ~e-19)
    # desired_yticks_mol = [4e-19, 5e-19, 6e-19, 7e-19, 8e-19]
    desired_yticks_mol = [3e-19, 4e-19, 5e-19, 6e-19, 7e-19]
    recep_rows = np.array(pivot.index)  # descending, in mol/cell

    wanted_positions_y = []
    wanted_labels_y = []
    for tick_val in desired_yticks_mol:
        idx = np.argmin(np.abs(recep_rows - tick_val))
        wanted_positions_y.append(idx + 0.5)  # centre of cell
        wanted_labels_y.append(f"{tick_val / 1e-19:.0f}")  # display as integer × 1e-19

    ax.set_yticks(wanted_positions_y)
    ax.set_yticklabels(wanted_labels_y, rotation=0)

    # Exponent label
    ax.text(
        0.12,
        1.05,
        r"$\times 10^{-10}$",
        transform=ax.transAxes,
        fontsize=7,
        ha="right",
    )
    # Desired tick values in nmol
    desired_ticks = [12.5, 50, 100, 150, 200]
    dose_cols = np.array(sorted(pivot.columns))

    wanted_positions = []
    wanted_labels = []
    for tick_val in desired_ticks:
        # Find the column index closest to the desired tick value
        idx = np.argmin(np.abs(dose_cols - tick_val))
        if np.abs(dose_cols[idx] - tick_val) < 1.0:  # tolerance of 1 nmol
            wanted_positions.append(idx + 0.5)  # seaborn cell centres at i+0.5
            wanted_labels.append(f"{tick_val:g}")  # :g drops trailing zeros

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
