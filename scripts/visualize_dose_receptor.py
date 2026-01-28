#!/usr/bin/env python3
"""
Visualizes dose-receptor sweep results as heatmap

Usage:
    python visualize_dose_receptor.py                    # Auto-find latest
    python visualize_dose_receptor.py 2026-01-27_15-30-00  # Specific timestamp
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import glob

def find_sweep_dir(timestamp=None):
    """Find sweep directory"""
    if timestamp:
        sweep_dir = f"results/DoseReceptorSweep_{timestamp}"
        if not Path(sweep_dir).exists():
            print(f"ERROR: Directory not found: {sweep_dir}")
            sys.exit(1)
        return sweep_dir
    
    # Find most recent
    sweep_dirs = glob.glob("results/DoseReceptorSweep_*")
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
    print(f"\nParameter ranges:")
    print(f"  Doses: {sorted(df['dose_nmol'].unique())} nmol")
    print(f"  Receptor densities: {sorted(df['receptors_per_cell_mol'].unique())}")
    print(f"\nOutcome summary:")
    print(df['outcome'].value_counts())
    
    return df, sweep_dir

def create_heatmap(df, output_path):
    """Create cure rate heatmap: dose vs receptor density"""
    
    # Calculate cure rate per parameter combination
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    cure_rates = df.groupby(['dose_nmol', 'receptors_per_cell_mol'])['is_cure'].mean().reset_index()
    cure_rates.columns = ['dose_nmol', 'receptors_per_cell_mol', 'cure_rate']
    
    # Pivot for heatmap
    pivot = cure_rates.pivot(index='receptors_per_cell_mol', 
                            columns='dose_nmol', 
                            values='cure_rate')
    
    # Sort (high receptors at top)
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot[sorted(pivot.columns)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Custom colormap: red (failure) → yellow → green (cure)
    cmap = plt.cm.RdYlGn
    
    # Plot heatmap
    sns.heatmap(pivot, 
                annot=True,
                fmt='.2f',
                cmap=cmap,
                cbar_kws={'label': 'Cure Rate'},
                linewidths=1,
                linecolor='white',
                vmin=0, vmax=1,
                ax=ax)
    
    # Labels
    ax.set_xlabel('Total Dose (nmol)', fontsize=14)
    ax.set_ylabel('Receptor Density (mol/cell)', fontsize=14)
    
    # Format y-axis labels as scientific notation
    yticklabels = [f'{float(label.get_text()):.1e}' for label in ax.get_yticklabels()]
    ax.set_yticklabels(yticklabels, rotation=0)
    
    if 'replicate' in df.columns:
        num_reps = df['replicate'].max()
        title = f'RPT Treatment Outcome vs Dose and Receptor Expression\n(n={num_reps} replicates per point)'
    else:
        title = 'RPT Treatment Outcome vs Dose and Receptor Expression'
    
    ax.set_title(title, fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nHeatmap saved to: {output_path}")
    
    return fig

def create_dose_response_curves(df, output_path):
    """Create dose-response curves for each receptor density"""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Calculate cure rates
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    
    receptor_levels = sorted(df['receptors_per_cell_mol'].unique())
    
    for recep in receptor_levels:
        subset = df[df['receptors_per_cell_mol'] == recep]
        cure_by_dose = subset.groupby('dose_nmol')['is_cure'].mean()
        
        label = f'{recep:.1e} mol/cell'
        ax.plot(cure_by_dose.index, cure_by_dose.values, 
                marker='o', linewidth=2, markersize=8, label=label)
    
    ax.set_xlabel('Total Dose (nmol)', fontsize=12)
    ax.set_ylabel('Cure Rate', fontsize=12)
    ax.set_title('Dose-Response Curves by Receptor Density', fontsize=14)
    ax.legend(title='Receptor Density', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dose-response curves saved to: {output_path}")
    
    return fig

def print_summary_stats(df):
    """Print useful summary statistics"""
    print("\n=== SUMMARY STATISTICS ===")
    
    # Overall cure rate
    overall_cure = (df['outcome'] == 'CURE').mean() * 100
    print(f"Overall cure rate: {overall_cure:.1f}%")
    
    # Cure rate by dose
    print("\nCure rate by dose:")
    dose_cure = df.groupby('dose_nmol')['outcome'].apply(
        lambda x: (x == 'CURE').mean() * 100
    )
    for dose, rate in dose_cure.items():
        print(f"  {dose:.0f} nmol: {rate:.1f}%")
    
    # Cure rate by receptor density
    print("\nCure rate by receptor density:")
    recep_cure = df.groupby('receptors_per_cell_mol')['outcome'].apply(
        lambda x: (x == 'CURE').mean() * 100
    )
    for recep, rate in recep_cure.items():
        print(f"  {recep:.1e} mol/cell: {rate:.1f}%")
    
    # Find approximate thresholds
    print("\n=== THRESHOLD ANALYSIS ===")
    
    # For baseline receptor density (1.0e-15), what dose gives ~50% cure?
    baseline_recep = 1.0e-15
    baseline_data = df[df['receptors_per_cell_mol'] == baseline_recep]
    if len(baseline_data) > 0:
        baseline_cure = baseline_data.groupby('dose_nmol')['outcome'].apply(
            lambda x: (x == 'CURE').mean()
        )
        print(f"\nAt baseline receptor density ({baseline_recep:.1e} mol/cell):")
        for dose, rate in baseline_cure.items():
            print(f"  {dose:.0f} nmol → {rate*100:.0f}% cure rate")
    
    # For a fixed dose (say 100 nmol), what receptor density gives ~50% cure?
    fixed_dose = 100.0
    dose_data = df[df['dose_nmol'] == fixed_dose]
    if len(dose_data) > 0:
        dose_cure = dose_data.groupby('receptors_per_cell_mol')['outcome'].apply(
            lambda x: (x == 'CURE').mean()
        )
        print(f"\nAt {fixed_dose:.0f} nmol dose:")
        for recep, rate in dose_cure.items():
            print(f"  {recep:.1e} mol/cell → {rate*100:.0f}% cure rate")

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
    heatmap_path = f"{sweep_dir}/heatmap_dose_receptor.png"
    create_heatmap(df, heatmap_path)
    
    # Dose-response curves
    curves_path = f"{sweep_dir}/dose_response_curves.png"
    create_dose_response_curves(df, curves_path)
    
    # Print statistics
    print_summary_stats(df)
    
    print("\n=== VISUALIZATION COMPLETE ===")
    print(f"All outputs in: {sweep_dir}")

if __name__ == "__main__":
    main()
