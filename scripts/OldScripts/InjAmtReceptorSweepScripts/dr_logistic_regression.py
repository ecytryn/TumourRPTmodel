#!/usr/bin/env python3
"""
Logistic regression on DR sweep outcomes.

Fits a linear decision boundary in (dose_nmol, receptors_per_cell_mol) space,
overlays it on the cure rate heatmap, and reports coefficients.

Usage:
    python dr_logistic_regression.py
    python dr_logistic_regression.py --timestamp 2026-03-14_21-39-50
    python dr_logistic_regression.py --output results/
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import glob
import argparse
import os
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score

mpl.rcParams.update({
    "figure.figsize": (6.7, 4.8),
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

# =============================================================================
# DATA
# =============================================================================

def find_sweep_dir(timestamp=None):
    base = "results/DoseReceptorSweep/DoseReceptorSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            raise FileNotFoundError(f"Not found: {d}")
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        raise FileNotFoundError("No DoseReceptorSweep directories found.")
    return dirs[-1]


def load_data(sweep_dir):
    csv_path = Path(sweep_dir) / "sweep_summary.csv"
    df = pd.read_csv(csv_path)
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    print(f"Loaded {len(df)} rows. Overall cure rate: {df['is_cure'].mean():.2%}")
    return df


# =============================================================================
# LOGISTIC REGRESSION
# =============================================================================

def fit_logistic_regression(df):
    X = df[['dose_nmol', 'receptors_per_cell_mol']].values
    y = df['is_cure'].values

    # Scale features for numerical stability
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    model = LogisticRegression(max_iter=1000)
    model.fit(X_scaled, y)

    y_pred = model.predict(X_scaled)
    acc  = accuracy_score(y, y_pred)
    auroc = roc_auc_score(y, model.predict_proba(X_scaled)[:, 1])

    print(f"\nLogistic Regression Results:")
    print(f"  Accuracy: {acc:.3f}")
    print(f"  AUROC:    {auroc:.3f}")

    # Unscale coefficients for interpretability
    coef_scaled = model.coef_[0]
    coef_original = coef_scaled / scaler.scale_
    intercept_original = (model.intercept_[0]
                          - np.dot(coef_scaled, scaler.mean_ / scaler.scale_))

    print(f"\nCoefficients (in original units):")
    print(f"  dose_nmol:              {coef_original[0]:.4e}")
    print(f"  receptors_per_cell_mol: {coef_original[1]:.4e}")
    print(f"  intercept:              {intercept_original:.4e}")

    # Relative importance (using scaled coefficients)
    total = np.sum(np.abs(coef_scaled))
    print(f"\nRelative importance (scaled):")
    print(f"  dose:      {abs(coef_scaled[0])/total*100:.1f}%")
    print(f"  receptors: {abs(coef_scaled[1])/total*100:.1f}%")

    return model, scaler, coef_original, intercept_original


def boundary_receptor_at_dose(dose, coef, intercept):
    """
    Decision boundary: coef[0]*dose + coef[1]*receptors + intercept = 0
    Solve for receptors given dose.
    """
    if abs(coef[1]) < 1e-30:
        return None
    return -(coef[0] * dose + intercept) / coef[1]


# =============================================================================
# HEATMAP WITH BOUNDARY
# =============================================================================

def create_heatmap_with_boundary(df, coef, intercept, output_path):
    # --- Build pivot table (same as dose_receptor_visualize.py) ---
    cure_rates = df.groupby(['dose_nmol', 'receptors_per_cell_mol'])['is_cure'].mean().reset_index()
    cure_rates.columns = ['dose_nmol', 'receptors_per_cell_mol', 'cure_rate']

    pivot = cure_rates.pivot(index='receptors_per_cell_mol',
                             columns='dose_nmol',
                             values='cure_rate')
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot[sorted(pivot.columns)]

    # Axis values in parameter space
    dose_vals   = np.array(sorted(pivot.columns))       # x-axis
    recep_vals  = np.array(sorted(pivot.index, reverse=True))  # y-axis (high at top)

    fig, ax = plt.subplots()

    cmap = plt.cm.RdYlGn
    sns.heatmap(pivot, annot=False, cmap=cmap,
                cbar_kws={'label': 'Cure Rate'},
                linewidths=1, linecolor='white',
                vmin=0, vmax=1, ax=ax)

    # --- Y-axis formatting (match dose_receptor_visualize.py) ---
    ytick_values = pivot.index.values
    common_exponent = int(np.floor(np.log10(np.abs(ytick_values[0]))))
    yticklabels = []
    for i, val in enumerate(ytick_values):
        yticklabels.append(f'{val / 10**common_exponent:.1f}' if i % 2 == 0 else '')
    ax.set_yticklabels(yticklabels, rotation=0)
    ax.text(0.03, 1.01, f'$\\times 10^{{{common_exponent}}}$',
            transform=ax.transAxes, fontsize=7, ha='right')

    xticklabels = []
    for i, label in enumerate(ax.get_xticklabels()):
        xticklabels.append(label.get_text() if i % 2 == 0 else '')
    ax.set_xticklabels(xticklabels, rotation=0)

    ax.set_xlabel('Total Injected Amount (nmol)')
    ax.set_ylabel('Receptor Density (mol/cell)')
    ax.set_title('RPT Treatment Outcome')

    # --- Decision boundary in heatmap coordinates ---
    # Seaborn heatmap cell centres are at 0.5, 1.5, 2.5, ...
    # x cell index: dose_vals[i] maps to x = i + 0.5
    # y cell index: recep_vals[j] maps to y = j + 0.5 (recep_vals sorted high-to-low)

    dose_fine = np.linspace(dose_vals.min(), dose_vals.max(), 300)
    recep_boundary = np.array([boundary_receptor_at_dose(d, coef, intercept)
                                for d in dose_fine])

    # Convert to heatmap coordinates
    # x: interpolate dose -> cell index
    x_hm = np.interp(dose_fine, dose_vals, np.arange(len(dose_vals))) + 0.5

    # y: recep_vals is high-to-low, so higher receptor = lower y index
    # interpolate receptor -> cell index (note reversed axis)
    recep_vals_asc = recep_vals[::-1]  # low to high
    y_indices_asc  = np.arange(len(recep_vals))[::-1]  # corresponding cell indices
    y_hm = np.interp(recep_boundary, recep_vals_asc, y_indices_asc[::-1].astype(float)) + 0.5

    # Mask points outside the receptor range
    valid = ((recep_boundary >= recep_vals.min()) &
             (recep_boundary <= recep_vals.max()))

    ax.plot(x_hm[valid], y_hm[valid], 'k--', linewidth=1.5,
            label='Decision boundary (p=0.5)', zorder=5)
    ax.legend(loc='lower right', framealpha=0.9)

    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--timestamp', default=None)
    parser.add_argument('--output',    default='results')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    sweep_dir = find_sweep_dir(args.timestamp)
    print(f"Using: {sweep_dir}\n")

    df = load_data(sweep_dir)
    model, scaler, coef, intercept = fit_logistic_regression(df)

    out_path = Path(args.output) / "dr_heatmap_with_boundary.pdf"
    create_heatmap_with_boundary(df, coef, intercept, out_path)


if __name__ == "__main__":
    main()
