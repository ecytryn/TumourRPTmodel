#!/usr/bin/env python3
"""
Unified analysis of IS and DR sweeps via effective dose capture capacity:

    C_eff = R_total * f(r_eff)

where:
    R_total = receptors_per_cell_mol * tumour_size   (mol)
    r_eff   = sqrt(tumour_size / pi) * CELL_LENGTH   (mm)
    f       = TJ beta retention fraction at r_eff

Plots overlapping histograms of C_eff for cure vs failure,
pooling runs from both sweeps.

Usage:
    python unified_capture_capacity.py
    python unified_capture_capacity.py --is_timestamp 2026-03-19_10-00-00
                                       --dr_timestamp 2026-03-20_09-00-00
    python unified_capture_capacity.py --min_interval 10   # filter IS sweep
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import sys
import argparse
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
# PHYSICAL CONSTANTS (must match SimParams.java)
# =============================================================================

CELL_LENGTH_MM = 0.01      # 10 um in mm
FIRST_INJECTION_DAY = 5    # same in both sweeps

# Cell type column indices in populations.csv
NORMAL    = 1
HYPOXIC   = 2
NECROTIC  = 3
APOPTOTIC = 4

# =============================================================================
# TJ BETA RETENTION TABLE (from BetaRetention.java, RETENTION_FRACTIONS_TJ)
# =============================================================================

R_VALUES_MM = np.array([
    0.1000, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000,
    0.7000, 0.8000, 0.9000, 1.0000, 1.2000, 1.4000,
    1.6000, 1.8000, 2.0000, 2.5000, 3.0000, 3.5000,
    4.0000, 4.5000, 5.0000
])

RETENTION_TJ = np.array([
    0.1982, 0.2899, 0.3641, 0.4275, 0.4861, 0.5408,
    0.5921, 0.6420, 0.6898, 0.7342, 0.8137, 0.8760,
    0.9187, 0.9449, 0.9593, 0.9701, 0.9700, 0.9695,
    0.9700, 0.9714, 0.9734
])

C_ASYMP = R_VALUES_MM[-1] * (1.0 - RETENTION_TJ[-1])  # for large-R extrapolation


def retention_fraction(R_mm):
    """Python equivalent of BetaRetention.getRetentionFraction()"""
    if R_mm <= 0:
        return 0.0
    if R_mm <= R_VALUES_MM[0]:
        return RETENTION_TJ[0] * (R_mm / R_VALUES_MM[0])**3
    if R_mm >= R_VALUES_MM[-1]:
        return 1.0 - C_ASYMP / R_mm
    # Linear interpolation
    idx = np.searchsorted(R_VALUES_MM, R_mm) - 1
    t = (R_mm - R_VALUES_MM[idx]) / (R_VALUES_MM[idx+1] - R_VALUES_MM[idx])
    return RETENTION_TJ[idx] + t * (RETENTION_TJ[idx+1] - RETENTION_TJ[idx])


def tumour_size_from_row(row):
    return int(row[NORMAL] + row[HYPOXIC] + row[NECROTIC] + row[APOPTOTIC])


def r_eff_mm(tumour_size):
    """Effective radius in mm from 2D cell count (circle approximation)."""
    return np.sqrt(tumour_size / np.pi) * CELL_LENGTH_MM


def capture_capacity(tumour_size, receptors_per_cell_mol):
    """C_eff = R_total * f(r_eff) in mol."""
    R_total = receptors_per_cell_mol * tumour_size
    f = retention_fraction(r_eff_mm(tumour_size))
    return R_total * f


# =============================================================================
# DIRECTORY FINDERS
# =============================================================================

def find_latest(pattern):
    dirs = sorted(glob.glob(pattern))
    if not dirs:
        raise FileNotFoundError(f"No directories found matching: {pattern}")
    return dirs[-1]


def find_sweep_dir(base, timestamp=None):
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            raise FileNotFoundError(f"Not found: {d}")
        return d
    return find_latest(f"{base}_*")


# =============================================================================
# DATA EXTRACTION
# =============================================================================

def load_populations(pop_file, row_idx):
    """Load populations.csv and return tumour size at given row index."""
    pops = np.loadtxt(pop_file, delimiter=',', skiprows=1)
    if row_idx >= pops.shape[0]:
        # Early stop - use last available row
        row = pops[-1]
    else:
        row = pops[row_idx]
    return tumour_size_from_row(row)


def read_receptors_from_params(params_file):
    """Read receptors_per_cell_mol from a parameters.csv file."""
    params = pd.read_csv(params_file, header=0,
                         names=['parameter', 'value', 'units', 'description'])
    row = params[params['parameter'] == 'receptors_per_cell_mol']
    if len(row) == 0:
        raise ValueError(f"receptors_per_cell_mol not found in {params_file}")
    return float(row['value'].values[0])


def extract_is_data(sweep_dir, min_interval=0):
    """Extract tumour size at second injection from IS sweep."""
    records = []
    run_dirs = sorted(Path(sweep_dir).glob("interval_*_skew_*_rep_*"))
    print(f"IS sweep: found {len(run_dirs)} run directories.")

    for run_dir in run_dirs:
        parts = run_dir.name.split('_')
        try:
            interval = int(parts[1])
            skew     = int(parts[3])
            rep      = int(parts[5])
        except (IndexError, ValueError):
            continue

        if interval < min_interval:
            continue

        pop_file = run_dir / "populations.csv"
        if not pop_file.exists():
            continue

        try:
            second_injection_day = FIRST_INJECTION_DAY + interval
            row_idx = second_injection_day * 24
            tumour_size = load_populations(pop_file, row_idx)
        except Exception as e:
            print(f"  Warning: {run_dir.name}: {e}")
            continue

        records.append({
            'interval': interval,
            'skew':     skew,
            'replicate': rep,
            'tumour_size': tumour_size,
            'sweep': 'IS',
        })

    df = pd.DataFrame(records)

    # Merge outcomes from sweep_summary.csv
    summary = pd.read_csv(Path(sweep_dir) / "sweep_summary.csv")
    summary['skew_int'] = summary['skew'].astype(int)
    df = df.merge(
        summary[['interval', 'skew_int', 'replicate', 'outcome']],
        left_on=['interval', 'skew', 'replicate'],
        right_on=['interval', 'skew_int', 'replicate'],
        how='left'
    )

    # Receptor density is fixed in IS sweep - read from first parameters.csv
    param_files = list(Path(sweep_dir).glob("*/parameters.csv"))
    receptors = read_receptors_from_params(param_files[0])
    df['receptors_per_cell_mol'] = receptors

    print(f"IS sweep: extracted {len(df)} records "
          f"(min_interval={min_interval}, receptors={receptors:.3e} mol/cell).")
    return df


def extract_dr_data(sweep_dir):
    """Extract tumour size at injection from DR sweep.
    
    Iterates over sweep_summary.csv to guarantee unambiguous folder matching.
    """
    summary = pd.read_csv(Path(sweep_dir) / "sweep_summary.csv")
    records = []
    row_idx = FIRST_INJECTION_DAY * 24
    n_missing = 0

    print(f"DR sweep: processing {len(summary)} rows from sweep_summary.csv.")

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
            tumour_size = load_populations(pop_file, row_idx)
        except Exception as e:
            print(f"  Warning: {run_dir.name}: {e}")
            continue

        records.append({
            'tumour_size':            tumour_size,
            'receptors_per_cell_mol': receptors,
            'outcome':                outcome,
            'sweep':                  'DR',
        })

    df = pd.DataFrame(records)
    if n_missing > 0:
        print(f"DR sweep: {n_missing} folders not found (skipped).")
    print(f"DR sweep: extracted {len(df)} records.")
    return df



# =============================================================================
# PLOT
# =============================================================================

def plot_unified(df_combined, output_dir, label='combined'):
    df_combined['C_eff'] = df_combined.apply(
        lambda r: capture_capacity(r['tumour_size'], r['receptors_per_cell_mol']),
        axis=1
    )

    cure    = df_combined.loc[df_combined['outcome'] == 'CURE',    'C_eff']
    failure = df_combined.loc[df_combined['outcome'] == 'FAILURE', 'C_eff']
    n_cure    = len(cure)
    n_failure = len(failure)
    n_total   = n_cure + n_failure

    print(f"\n{label}: {n_total} runs ({n_cure} cure, {n_failure} failure)")
    if n_cure > 0:
        print(f"  Cure:    median={cure.median():.3e}, mean={cure.mean():.3e}")
    if n_failure > 0:
        print(f"  Failure: median={failure.median():.3e}, mean={failure.mean():.3e}")

    all_c = df_combined['C_eff']
    bins = np.linspace(0, np.percentile(all_c, 99) * 1.05, 40)
    bin_width = bins[1] - bins[0]

    fig, ax = plt.subplots()

    if n_failure > 0:
        ax.hist(failure, bins=bins,
                weights=np.ones(n_failure) / (n_total * bin_width),
                alpha=0.6, color='#E84855', label=f'Failure (n={n_failure})', zorder=2)
    if n_cure > 0:
        ax.hist(cure, bins=bins,
                weights=np.ones(n_cure) / (n_total * bin_width),
                alpha=0.6, color='#2E86AB', label=f'Cure (n={n_cure})', zorder=3)

    ax.set_xlabel('Effective dose capture capacity $C_{\\mathrm{eff}}$ (mol)')
    ax.set_ylabel('Probability density')
    ax.set_title(label)
    ax.legend(framealpha=0.95)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

    plt.tight_layout()

    out_dir = Path(output_dir)
    stem = f"capture_capacity_hist_{label.replace(' ', '_')}"
    plt.savefig(out_dir / f"{stem}.pdf", bbox_inches='tight')
    plt.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches='tight')
    print(f"Saved to {out_dir}/{stem}.pdf")

    df_combined.to_csv(out_dir / f"{stem}.csv", index=False)

    plt.show()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--is_timestamp',  default=None)
    parser.add_argument('--dr_timestamp',  default=None)
    parser.add_argument('--min_interval',  type=int, default=0,
                        help='Exclude IS runs with interval < this (default: 0)')
    parser.add_argument('--sweep', choices=['IS', 'DR', 'both'], default='both',
                        help='Which sweep(s) to include (default: both)')
    parser.add_argument('--output', default='results',
                        help='Output directory (default: results/)')
    args = parser.parse_args()

    import os
    os.makedirs(args.output, exist_ok=True)

    cols = ['tumour_size', 'receptors_per_cell_mol', 'outcome', 'sweep']

    if args.sweep in ('IS', 'both'):
        is_dir = find_sweep_dir("results/IntervalSkewSweep/IntervalSkewSweep",
                                args.is_timestamp)
        print(f"IS sweep: {is_dir}")
        df_is = extract_is_data(is_dir, min_interval=args.min_interval)
        df_is = df_is.dropna(subset=['outcome'])

    if args.sweep in ('DR', 'both'):
        dr_dir = find_sweep_dir("results/DoseReceptorSweep/DoseReceptorSweep",
                                args.dr_timestamp)
        print(f"DR sweep: {dr_dir}")
        df_dr = extract_dr_data(dr_dir)
        df_dr = df_dr.dropna(subset=['outcome'])

    if args.sweep == 'IS':
        plot_unified(df_is[cols], args.output, label='IS sweep')
    elif args.sweep == 'DR':
        plot_unified(df_dr[cols], args.output, label='DR sweep')
    else:
        # Both individually and combined
        plot_unified(df_is[cols], args.output, label='IS sweep')
        plot_unified(df_dr[cols], args.output, label='DR sweep')
        df_combined = pd.concat([df_is[cols], df_dr[cols]], ignore_index=True)
        plot_unified(df_combined, args.output, label='combined')


if __name__ == "__main__":
    main()
