#!/usr/bin/env python3
"""
Compare interval-skew sweep sections across different vessel densities

Usage:
    python compare_sweep_sections.py --mode row --value 0
    python compare_sweep_sections.py --mode col --value 24

    # Specify exact sweep timestamps
    python compare_sweep_sections.py --mode row --value 0 \
        --timestamps 2026-02-09_14-30-00 2026-02-09_15-45-00

    # Change vessel configs to search for
    python compare_sweep_sections.py --mode row --value 0 \
        --vessel-configs Rrepel20 Rrepel40 Rrepel60
    
    # Use different replicate
    python compare_sweep_sections.py --mode row --value 0 --replicate 2
    
    # Custom output filename
    python compare_sweep_sections.py --mode row --value 0 \
        --output my_comparison.pdf
        
This creates a grid of plots comparing population + dose rate dynamics across vessel densities
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import glob

# Cell type indices
HEALTHY = 0
NORMAL = 1
HYPOXIC = 2
NECROTIC = 3
APOPTOTIC = 4
VESSEL = 5

def find_sweep_dir(timestamp=None, vessel_config=None):
    """Find the sweep directory"""
    if timestamp:
        if vessel_config:
            pattern = f"results/IntervalSkewSweep/IntervalSkewSweep_{timestamp}_*{vessel_config}*"
        else:
            pattern = f"results/IntervalSkewSweep/IntervalSkewSweep_{timestamp}"
        
        dirs = glob.glob(pattern)
        if dirs:
            return dirs[0]
    
    # Auto-find most recent
    if vessel_config:
        pattern = f"results/IntervalSkewSweep/IntervalSkewSweep_*{vessel_config}*"
    else:
        pattern = "results/IntervalSkewSweep/IntervalSkewSweep_*"
    
    dirs = sorted(glob.glob(pattern))
    if not dirs:
        raise FileNotFoundError(f"No sweep directories found matching {pattern}")
    
    return dirs[-1]

def load_run_data(sweep_dir, interval, skew, replicate):
    """Load populations and dose rate data for a single run"""
    run_dir = f"{sweep_dir}/interval_{interval}_skew_{int(skew)}_rep_{replicate}"
    
    pop_file = f"{run_dir}/populations.csv"
    dose_file = f"{run_dir}/dose.csv"
    
    if not Path(pop_file).exists():
        print(f"Warning: {pop_file} not found")
        return None, None
    
    populations = np.loadtxt(pop_file, delimiter=',')
    
    # Dose file might not exist
    dose_rate = None
    if Path(dose_file).exists():
        dose_rate = np.loadtxt(dose_file, delimiter=',')
        # Handle both single column and multi-column dose files
        if len(dose_rate.shape) > 1:
            dose_rate = dose_rate[:, 0]  # Take first column if multi-column
    
    return populations, dose_rate

def get_sweep_parameters(sweep_dir):
    """Extract sweep parameters from directory structure"""
    # Look at subdirectories to find parameter ranges
    subdirs = list(Path(sweep_dir).glob("interval_*_skew_*_rep_*"))
    
    intervals = set()
    skews = set()
    
    for d in subdirs:
        parts = d.name.split('_')
        interval = int(parts[1])
        skew = int(parts[3])
        
        intervals.add(interval)
        skews.add(skew)
    
    return sorted(intervals), sorted(skews, reverse=True)

def plot_single_run(ax_pop, ax_dose, populations, dose_rate, interval, skew, vessel_label):
    """Plot population and dose rate data for a single run"""
    
    t_hours = np.arange(populations.shape[0])
    t_days = t_hours / 24.0
    
    # Plot populations
    tumor_pop = populations[:, NORMAL] + populations[:, HYPOXIC]
    
    ax_pop.plot(t_days, populations[:, NORMAL], 'b-', linewidth=1, alpha=0.7, label='Normoxic')
    ax_pop.plot(t_days, populations[:, HYPOXIC], 'orange', linewidth=1, alpha=0.7, label='Hypoxic')
    ax_pop.plot(t_days, tumor_pop, 'k-', linewidth=2, label='Total tumor')
    
    ax_pop.set_ylabel('Cell count', fontsize=8)
    ax_pop.tick_params(labelsize=7)
    ax_pop.set_xlim(0, max(t_days))
    
    # Add title with parameters
    title = f"I={interval}d, S={skew}nm\n{vessel_label}"
    ax_pop.set_title(title, fontsize=8)
    
    # Only show legend on first plot
    if interval == 20 and skew == 20:  # Adjust based on your data
        ax_pop.legend(fontsize=6, loc='upper right')
    
    # Plot dose rate if available
    if dose_rate is not None:
        # Dose rate in Gy/hour
        ax_dose.fill_between(t_days, dose_rate, step='pre', alpha=0.3, color='red', label='Dose rate')
        ax_dose.plot(t_days, dose_rate, 'r-', linewidth=0.5, alpha=0.5)
        
        ax_dose.set_ylabel('Dose rate (Gy/hr)', fontsize=8)
        ax_dose.tick_params(labelsize=7)
        ax_dose.set_xlim(0, max(t_days))
        
        # Only show legend on first plot
        if interval == 20 and skew == 20:
            ax_dose.legend(fontsize=6, loc='upper right')
    else:
        ax_dose.text(0.5, 0.5, 'No dose data', 
                   ha='center', va='center', fontsize=8, transform=ax_dose.transAxes)
    
    # Only show x-label on bottom row
    ax_dose.set_xlabel('Time (days)', fontsize=8)

def main():
    parser = argparse.ArgumentParser(description='Compare sweep sections across vessel densities')
    parser.add_argument('--mode', choices=['row', 'col'], required=True,
                        help='Compare across row (fixed skew, varying interval) or column (fixed interval, varying skew)')
    parser.add_argument('--value', type=float, required=True,
                        help='Skew value (nmol) for row mode, or interval value (days) for col mode')
    parser.add_argument('--timestamps', nargs='+', default=None,
                        help='Timestamps of sweep directories (auto-detect if not provided)')
    parser.add_argument('--vessel-configs', nargs='+', default=['Rrepel20', 'Rrepel50', 'Rrepel60'],
                        help='Vessel configuration identifiers (default: Rrepel20 Rrepel50 Rrepel60)')
    parser.add_argument('--replicate', type=int, default=1,
                        help='Which replicate to visualize (default: 1)')
    parser.add_argument('--output', default=None,
                        help='Output filename (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Find sweep directories for each vessel density
    sweep_dirs = []
    for vessel_config in args.vessel_configs:
        timestamp = args.timestamps[len(sweep_dirs)] if args.timestamps and len(sweep_dirs) < len(args.timestamps) else None
        sweep_dir = find_sweep_dir(timestamp, vessel_config)
        sweep_dirs.append((vessel_config, sweep_dir))
        print(f"Using {vessel_config}: {sweep_dir}")
    
    # Get parameters from first sweep
    _, first_sweep = sweep_dirs[0]
    intervals, skews = get_sweep_parameters(first_sweep)
    
    print(f"\nAvailable intervals: {intervals}")
    print(f"Available skews: {skews}")
    
    # Determine which parameters to plot
    if args.mode == 'row':
        # Fixed skew, varying interval
        fixed_skew = args.value
        if fixed_skew not in skews:
            print(f"\nWarning: Skew {fixed_skew} not found. Using closest: {min(skews, key=lambda x: abs(x-fixed_skew))}")
            fixed_skew = min(skews, key=lambda x: abs(x-fixed_skew))
        
        plot_params = [(interval, fixed_skew) for interval in intervals]
        section_label = f"skew_{int(fixed_skew)}nm"
        
    else:  # col mode
        # Fixed interval, varying skew
        fixed_interval = int(args.value)
        if fixed_interval not in intervals:
            print(f"\nWarning: Interval {fixed_interval} not found. Using closest: {min(intervals, key=lambda x: abs(x-fixed_interval))}")
            fixed_interval = min(intervals, key=lambda x: abs(x-fixed_interval))
        
        plot_params = [(fixed_interval, skew) for skew in skews]
        section_label = f"interval_{fixed_interval}d"
    
    # Create figure
    n_vessels = len(sweep_dirs)
    n_params = len(plot_params)
    
    # Two rows per vessel density (population + dose)
    fig, axes = plt.subplots(n_vessels * 2, n_params, 
                             figsize=(3 * n_params, 3 * n_vessels))
    
    # Ensure axes is 2D
    if n_vessels == 1:
        if n_params == 1:
            axes = np.array([[axes[0]], [axes[1]]])
        else:
            axes = axes.reshape(2, n_params)
    elif n_params == 1:
        axes = axes.reshape(n_vessels * 2, 1)
    
    # Plot each vessel density
    for vessel_idx, (vessel_config, sweep_dir) in enumerate(sweep_dirs):
        print(f"\nProcessing {vessel_config}...")
        
        for param_idx, (interval, skew) in enumerate(plot_params):
            row_pop = vessel_idx * 2
            row_dose = vessel_idx * 2 + 1
            
            ax_pop = axes[row_pop, param_idx]
            ax_dose = axes[row_dose, param_idx]
            
            # Load data
            populations, dose_rate = load_run_data(sweep_dir, interval, skew, args.replicate)
            
            if populations is not None:
                plot_single_run(ax_pop, ax_dose, populations, dose_rate, 
                               interval, skew, vessel_config)
            else:
                ax_pop.text(0.5, 0.5, 'No data', ha='center', va='center', 
                           transform=ax_pop.transAxes)
                ax_dose.text(0.5, 0.5, 'No data', ha='center', va='center', 
                          transform=ax_dose.transAxes)
    
    # Overall title
    mode_str = f"Fixed skew={args.value}nm" if args.mode == 'row' else f"Fixed interval={args.value}d"
    fig.suptitle(f'Vessel Density Comparison: {mode_str} (replicate {args.replicate})', 
                 fontsize=12, y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save figure
    if args.output:
        output_file = args.output
    else:
        output_file = f"vessel_comparison_{section_label}_rep{args.replicate}.pdf"
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved to: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
