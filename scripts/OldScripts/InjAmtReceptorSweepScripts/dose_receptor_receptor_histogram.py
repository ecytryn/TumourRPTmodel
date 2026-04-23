#!/usr/bin/env python3
"""
DR sweep analysis: total receptor count R_total = tumour_size * receptors_per_cell_mol

Tumour size is taken at the injection day (day 5).
Receptor density is read from sweep_summary.csv.

Optionally plot receptor density alone (--metric receptors) to compare
with the R_total version and assess how much tumour size variance matters.

Usage:
    python dr_receptor_histogram.py
    python dr_receptor_histogram.py --timestamp 2026-03-14_21-39-50
    python dr_receptor_histogram.py --metric receptors   # receptor density only
    python dr_receptor_histogram.py --output results/my_output
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import argparse
import os
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

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

FIRST_INJECTION_DAY = 5

# Cell type column indices in populations.csv
NORMAL    = 1
HYPOXIC   = 2
NECROTIC  = 3
APOPTOTIC = 4

# =============================================================================
# DATA EXTRACTION
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


def tumour_size_from_row(row):
    return int(row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC])


def load_tumour_size(pop_file, row_idx):
    pops = np.loadtxt(pop_file, delimiter=',', skiprows=1)
    row = pops[-1] if row_idx >= pops.shape[0] else pops[row_idx]
    return tumour_size_from_row(row)


def extract_data(sweep_dir):
    """Iterate over sweep_summary.csv rows for unambiguous folder matching."""
    summary = pd.read_csv(Path(sweep_dir) / "sweep_summary.csv")
    records = []
    row_idx = FIRST_INJECTION_DAY * 24
    n_missing = 0

    print(f"Processing {len(summary)} rows from sweep_summary.csv.")

    for _, row in summary.iterrows():
        dose      = row['dose_nmol']
        receptors = row['receptors_per_cell_mol']
        rep       = int(row['replicate'])
        outcome   = row['outcome']

        run_dir  = Path(sweep_dir) / f"dose_{dose:.0f}_recep_{receptors:.2e}_rep_{rep}"
        pop_file = run_dir / "populations.csv"

        if not pop_file.exists():
            n_missing += 1
            continue

        try:
            tumour_size = load_tumour_size(pop_file, row_idx)
        except Exception as e:
            print(f"  Warning: {run_dir.name}: {e}")
            continue

        records.append({
            'dose_nmol':              dose,
            'receptors_per_cell_mol': receptors,
            'tumour_size':            tumour_size,
            'R_total':                tumour_size * receptors,
            'replicate':              rep,
            'outcome':                outcome,
        })

    df = pd.DataFrame(records)
    if n_missing > 0:
        print(f"{n_missing} folders not found (skipped).")
    print(f"Extracted {len(df)} records.")
    return df


# =============================================================================
# PLOT
# =============================================================================

def plot_histogram(df, output_dir, metric='R_total'):
    """
    metric: 'R_total' = tumour_size * receptors_per_cell_mol
            'receptors' = receptors_per_cell_mol only
    """
    if metric == 'R_total':
        x_col  = 'R_total'
        xlabel = 'Total receptor count $R_{\\mathrm{total}}$ (mol)'
        stem   = 'dr_R_total_hist'
    else:
        x_col  = 'receptors_per_cell_mol'
        xlabel = 'Receptor density (mol/cell)'
        stem   = 'dr_receptors_hist'

    cure    = df.loc[df['outcome'] == 'CURE',    x_col]
    failure = df.loc[df['outcome'] == 'FAILURE', x_col]
    n_cure    = len(cure)
    n_failure = len(failure)
    n_total   = n_cure + n_failure

    print(f"\n{n_total} runs ({n_cure} cure, {n_failure} failure)")
    print(f"  Cure:    median={cure.median():.3e}, mean={cure.mean():.3e}")
    print(f"  Failure: median={failure.median():.3e}, mean={failure.mean():.3e}")

    all_x = df[x_col]
    bins = np.linspace(0, np.percentile(all_x, 99) * 1.05, 40)
    bin_width = bins[1] - bins[0]

    fig, ax = plt.subplots()

    ax.hist(failure, bins=bins,
            weights=np.ones(n_failure) / (n_total * bin_width),
            alpha=0.6, color='#E84855', label=f'Failure (n={n_failure})', zorder=2)
    ax.hist(cure, bins=bins,
            weights=np.ones(n_cure) / (n_total * bin_width),
            alpha=0.6, color='#2E86AB', label=f'Cure (n={n_cure})', zorder=3)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Probability density')
    ax.legend(framealpha=0.95)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

    plt.tight_layout()

    out_dir = Path(output_dir)
    plt.savefig(out_dir / f"{stem}.pdf", bbox_inches='tight')
    plt.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches='tight')
    print(f"Saved to {out_dir}/{stem}.pdf")

    df.to_csv(out_dir / f"{stem}.csv", index=False)

    plt.show()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--timestamp', default=None)
    parser.add_argument('--metric', choices=['R_total', 'receptors'], default='R_total',
                        help='X-axis metric: R_total (default) or receptors')
    parser.add_argument('--output', default='results')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    sweep_dir = find_sweep_dir(args.timestamp)
    print(f"Using: {sweep_dir}\n")

    df = extract_data(sweep_dir)
    plot_histogram(df, args.output, metric=args.metric)


if __name__ == "__main__":
    main()
