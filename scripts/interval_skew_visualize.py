#!/usr/bin/env python3
"""
Visualizes parameter sweep results as:
1. Clean heatmap with cure rate (greenish-yellow to green) + mean injections annotation
2. Small multiples grid showing population trajectories for each parameter combination

Usage:
    python visualize_interval_skew.py              # Uses default sweep_summary.csv
    python visualize_interval_skew.py _v2          # Uses sweep_summary_v2.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path
import sys
import os

# Usage: python visualize_interval_skew.py 2026-01-27_01-53-14
TIMESTAMP = sys.argv[1] if len(sys.argv) > 1 else None

if TIMESTAMP:
    CSV_FILE = f"results/IntervalSkewSweep/IntervalSkewSweep_{TIMESTAMP}/sweep_summary.csv"
    SWEEP_DIR = f"results/IntervalSkewSweep/IntervalSkewSweep_{TIMESTAMP}"
else:
    # Try to find most recent
    import glob
    sweep_dirs = glob.glob("results/IntervalSkewSweep/IntervalSkewSweep_*")
    if sweep_dirs:
        SWEEP_DIR = sorted(sweep_dirs)[-1]
        CSV_FILE = f"{SWEEP_DIR}/sweep_summary.csv"
    else:
        print("ERROR: No sweep found. Pass timestamp: python interval_skew_visualize.py 2026-01-27_01-53-14")
        sys.exit(1)

OUTPUT_DIR = SWEEP_DIR  # Put visualizations in same directory as sweep

def load_data(csv_path):
    """Load parameter sweep results"""
    df = pd.read_csv(csv_path)
    
    # Filter out invalid runs
    df = df[df['outcome'] != 'INVALID']
    
    print(f"Loaded {len(df)} valid simulation results")
    if 'replicate' in df.columns:
        num_replicates = df['replicate'].max()
        print(f"  Replicates per parameter set: {num_replicates}")
    
    print(f"Parameter ranges:")
    print(f"  Intervals: {sorted(df['interval'].unique())}")
    print(f"  Skews: {sorted(df['skew'].unique())}")
    print(f"\nOutcome summary:")
    print(df['outcome'].value_counts())
    
    return df

def create_custom_colormap():
    """Create greenish-yellow to green colormap"""
    # Yellow-green (#CCFF00ish) to Green (#00AA00ish)
    colors = [
        (0.9, 0.2, 0.2),   # Red for 0 (failure)
        (1.0, 0.9, 0.4),   # Yellow for ~0.3
        (0.7, 0.9, 0.3),   # Yellow-green for ~0.5
        (0.4, 0.8, 0.3),   # Light green for ~0.7
        (0.0, 0.6, 0.2),   # Green for 1.0 (success)
    ]
    positions = [0.0, 0.3, 0.5, 0.7, 1.0]
    
    cmap = LinearSegmentedColormap.from_list('cure_rate', list(zip(positions, colors)))
    return cmap

def create_cure_rate_heatmap(df, output_path):
    """Create heatmap with cure rate coloring and mean injections annotation"""
    
    # Calculate cure rate per parameter combination
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    cure_rates = df.groupby(['interval', 'skew'])['is_cure'].mean().reset_index()
    cure_rates.columns = ['interval', 'skew', 'cure_rate']
    
    # Calculate mean injections used per parameter combination
    mean_injections = df.groupby(['interval', 'skew'])['injectionsUsed'].mean().reset_index()
    mean_injections.columns = ['interval', 'skew', 'mean_injections']
    
    # Merge
    merged = cure_rates.merge(mean_injections, on=['interval', 'skew'])
    
    # Pivot for heatmap
    pivot_cure = merged.pivot(index='skew', columns='interval', values='cure_rate')
    pivot_inj = merged.pivot(index='skew', columns='interval', values='mean_injections')
    
    # Sort indices
    pivot_cure = pivot_cure.sort_index(ascending=False)
    pivot_cure = pivot_cure[sorted(pivot_cure.columns)]
    pivot_inj = pivot_inj.sort_index(ascending=False)
    pivot_inj = pivot_inj[sorted(pivot_inj.columns)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Custom colormap
    cmap = create_custom_colormap()
    
    # Create annotation array (mean injections, rounded to 1 decimal)
    annot_data = pivot_inj.round(1)
    
    # Plot heatmap
    sns.heatmap(pivot_cure, 
                annot=annot_data,
                fmt='.1f',
                cmap=cmap,
                cbar_kws={'label': 'Cure Rate'},
                linewidths=1,
                linecolor='white',
                vmin=0, vmax=1,
                ax=ax)
    
    # Labels and title
    ax.set_xlabel('Inter-injection Interval (days)', fontsize=14)
    ax.set_ylabel('Injection Skew s (nmol)', fontsize=14)
    
    if 'replicate' in df.columns:
        num_reps = df['replicate'].max()
        title = f'RPT Treatment Outcome\n(Color = Cure Rate, Number = Mean Injections Used, n={num_reps} replicates)'
    else:
        title = 'RPT Treatment Outcome\n(Color = Cure Rate, Number = Mean Injections Used)'
    
    ax.set_title(title + '\nInjection Pattern: [90+s, 90-s] nmol',
                fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
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
                data = np.loadtxt(pop_file, delimiter=',')
                if len(data.shape) == 1:
                    data = data.reshape(1, -1)
                
                # Sum normoxic (index 1) and hypoxic (index 2) tumor cells
                if data.shape[1] >= 3:
                    tumor_pop = data[:, 1] + data[:, 2]
                    trajectories.append(tumor_pop)
            except Exception as e:
                print(f"Warning: Could not load {pop_file}: {e}")
    
    return trajectories

def create_small_multiples(df, sweep_dir, output_path, max_time=70, max_pop=10000):
    """Create small multiples grid showing population trajectories"""
    
    intervals = sorted(df['interval'].unique())
    skews = sorted(df['skew'].unique(), reverse=True)  # High skew at top
    
    n_rows = len(skews)
    n_cols = len(intervals)
    
    # Calculate cure rates for background coloring
    df['is_cure'] = (df['outcome'] == 'CURE').astype(int)
    cure_rates = df.groupby(['interval', 'skew'])['is_cure'].mean()
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.5 * n_cols, 2 * n_rows))
    
    cmap = create_custom_colormap()
    
    num_reps = df['replicate'].max() if 'replicate' in df.columns else 1
    
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
#            ax.set_facecolor((*bg_color[:3], 0.3))  # Lighter version
            ax.set_facecolor(bg_color)  # Darker version            
            # Load and plot trajectories
            trajectories = load_population_data(sweep_dir, interval, skew, num_reps)
            
            if trajectories:
                for traj in trajectories:
                    # Time axis: each entry is 1 hour, plot in days
                    time_days = np.arange(len(traj)) / 24.0
                    ax.plot(time_days, traj, color='darkblue', alpha=0.5, linewidth=0.8)
            
            # Minimal formatting
            ax.set_xlim(0, max_time)
            ax.set_ylim(0, max_pop)
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Add interval label on bottom row
            if i == n_rows - 1:
                ax.set_xlabel(f'{interval}d', fontsize=9)
            
            # Add skew label on left column
            if j == 0:
                ax.set_ylabel(f's={int(skew)}', fontsize=9)
            
            # Add thin border
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_color('gray')
    
    # Add overall labels
    fig.text(0.5, 0.02, 'Inter-injection Interval', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'Injection Skew (nmol)', va='center', rotation='vertical', fontsize=14)
    
    # Title
    fig.suptitle(f'Tumor Population Trajectories\n(x: 0-{max_time} days, y: 0-{max_pop} cells, n={num_reps} replicates overlaid)',
                fontsize=14, y=1.02)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.08, bottom=0.08, top=0.95)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
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
    
    # Figure 1: Clean heatmap with cure rate + injections annotation
    create_cure_rate_heatmap(df, f"{OUTPUT_DIR}/cure_rate_grid.png")
    
    # Figure 2: Small multiples with population trajectories
    if os.path.exists(SWEEP_DIR):
        create_small_multiples(df, SWEEP_DIR, f"{OUTPUT_DIR}/sweep_summary.png")
    else:
        print(f"Warning: Sweep directory not found: {SWEEP_DIR}")
        print("Skipping small multiples figure.")
    
    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    overall_cure_rate = (df['outcome'] == 'CURE').mean() * 100
    print(f"Overall cure rate: {overall_cure_rate:.1f}%")
    
    print("\nCure rate by interval:")
    interval_cure = df.groupby('interval')['outcome'].apply(
        lambda x: (x == 'CURE').mean() * 100
    )
    print(interval_cure.to_string())
    
    print("\nCure rate by skew:")
    skew_cure = df.groupby('skew')['outcome'].apply(
        lambda x: (x == 'CURE').mean() * 100
    )
    print(skew_cure.to_string())
    
    print("\n=== VISUALIZATION COMPLETE ===")

if __name__ == "__main__":
    main()
