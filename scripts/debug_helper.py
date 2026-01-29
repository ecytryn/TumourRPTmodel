#!/usr/bin/env python3
"""
Debug Helper - Extract parameters from sweep heatmaps for investigation

This script helps you identify interesting points in your sweep heatmaps
and generates the command to debug them with full visualization.

USAGE:
    python debug_helper.py results/ParameterSweep_TIMESTAMP/
    python debug_helper.py results/DoseReceptorSweep_TIMESTAMP/
    
    # Interactive mode (select from heatmap)
    python debug_helper.py results/ParameterSweep_TIMESTAMP/ --interactive
    
    # Quick command generation
    python debug_helper.py results/ParameterSweep_TIMESTAMP/ --x 4 --y 6
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

def load_sweep_data(sweep_dir):
    """Load sweep results and identify sweep type"""
    sweep_path = Path(sweep_dir)
    
    # Find the summary CSV
    csv_files = list(sweep_path.glob("sweep_summary*.csv"))
    if not csv_files:
        print(f"ERROR: No sweep_summary.csv found in {sweep_dir}")
        sys.exit(1)
    
    csv_file = csv_files[0]
    print(f"Loading: {csv_file}")
    
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} sweep points\n")
    
    # Identify sweep type from columns
    # ParameterSweep columns: interval, skew, dose1_nmol, dose2_nmol, replicate, finalTumorCount, outcome, injectionsUsed
    # DoseReceptorSweep columns: dose_nmol, receptors_per_cell_mol, dose1_nmol, replicate, finalTumorCount, outcome
    
    if 'interval' in df.columns and 'skew' in df.columns:
        sweep_type = 'interval_skew'
        param1_name = 'interval'
        param2_name = 'skew'
    elif 'dose_nmol' in df.columns and 'receptors_per_cell_mol' in df.columns:
        sweep_type = 'dose_receptor'
        param1_name = 'dose_nmol'
        param2_name = 'receptors_per_cell_mol'
    else:
        print("ERROR: Could not identify sweep type from columns")
        print(f"Available columns: {df.columns.tolist()}")
        print("\nExpected columns:")
        print("  ParameterSweep: interval, skew, dose1_nmol, dose2_nmol, replicate, finalTumorCount, outcome")
        print("  DoseReceptorSweep: dose_nmol, receptors_per_cell_mol, dose1_nmol, replicate, finalTumorCount, outcome")
        sys.exit(1)
    
    return df, sweep_type, param1_name, param2_name

def create_heatmap(df, sweep_type, param1_name, param2_name):
    """Create heatmap showing cure/failure outcomes"""
    
    # Get unique parameter values
    param1_vals = sorted(df[param1_name].unique())
    param2_vals = sorted(df[param2_name].unique())
    
    # Create outcome matrix (average over replicates)
    outcome_matrix = np.zeros((len(param2_vals), len(param1_vals)))
    
    for i, p1 in enumerate(param1_vals):
        for j, p2 in enumerate(param2_vals):
            subset = df[(df[param1_name] == p1) & (df[param2_name] == p2)]
            if len(subset) > 0:
                # Cure rate (fraction of replicates with finalTumorCount < 10)
                cure_rate = (subset['finalTumorCount'] < 10).mean()
                outcome_matrix[j, i] = cure_rate
            else:
                outcome_matrix[j, i] = np.nan
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(outcome_matrix, aspect='auto', origin='lower', 
                   cmap='RdYlGn', vmin=0, vmax=1)
    
    # Set ticks
    ax.set_xticks(range(len(param1_vals)))
    ax.set_yticks(range(len(param2_vals)))
    
    if sweep_type == 'interval_skew':
        ax.set_xticklabels([f"{int(v)}" for v in param1_vals])
        ax.set_yticklabels([f"{int(v*1e9)}" for v in param2_vals])
        ax.set_xlabel('Interval (days)', fontsize=12)
        ax.set_ylabel('Skew (nmol)', fontsize=12)
    else:  # dose_receptor
        ax.set_xticklabels([f"{int(v)}" for v in param1_vals])
        # Format receptor density as percentage of baseline
        baseline = 6.64e-19
        ax.set_yticklabels([f"{int(v/baseline*100)}" for v in param2_vals])
        ax.set_xlabel('Dose (nmol)', fontsize=12)
        ax.set_ylabel('Receptor Density (% baseline)', fontsize=12)
    
    ax.set_title('Treatment Outcome Heatmap\n(Green = Cure, Red = Failure)', 
                 fontsize=14, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Cure Rate', rotation=270, labelpad=20)
    
    # Add text annotations
    for i in range(len(param1_vals)):
        for j in range(len(param2_vals)):
            if not np.isnan(outcome_matrix[j, i]):
                text_color = 'white' if outcome_matrix[j, i] < 0.5 else 'black'
                ax.text(i, j, f'{outcome_matrix[j, i]:.1f}',
                       ha="center", va="center", color=text_color, fontsize=8)
    
    plt.tight_layout()
    return fig, param1_vals, param2_vals

def interactive_selection(df, sweep_type, param1_name, param2_name):
    """Interactive selection of points from heatmap"""
    
    fig, param1_vals, param2_vals = create_heatmap(df, sweep_type, param1_name, param2_name)
    
    print("\n" + "="*60)
    print("INTERACTIVE POINT SELECTION")
    print("="*60)
    print("\nClick on a point in the heatmap to generate debug command")
    print("Close the window when done\n")
    
    selected_points = []
    
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            x_idx = int(round(event.xdata))
            y_idx = int(round(event.ydata))
            
            if 0 <= x_idx < len(param1_vals) and 0 <= y_idx < len(param2_vals):
                param1_val = param1_vals[x_idx]
                param2_val = param2_vals[y_idx]
                
                # Get actual outcomes at this point
                subset = df[(df[param1_name] == param1_val) & (df[param2_name] == param2_val)]
                
                if len(subset) > 0:
                    cure_rate = (subset['finalTumorCount'] < 10).mean()
                    avg_final = subset['finalTumorCount'].mean()
                    
                    print(f"\n--- Point ({x_idx}, {y_idx}) ---")
                    print(f"{param1_name}: {param1_val}")
                    print(f"{param2_name}: {param2_val}")
                    print(f"Cure rate: {cure_rate:.1f} ({int(cure_rate*len(subset))}/{len(subset)} replicates)")
                    print(f"Avg final pop: {avg_final:.1f} cells")
                    
                    # Generate command
                    cmd = generate_debug_command(sweep_type, param1_val, param2_val, 
                                                 x_idx, y_idx)
                    print(f"\nDebug command:")
                    print(f"  {cmd}")
                    
                    selected_points.append((x_idx, y_idx, param1_val, param2_val, cmd))
    
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    return selected_points

def generate_debug_command(sweep_type, param1_val, param2_val, x_idx=None, y_idx=None):
    """Generate the gradle command to debug this parameter combination"""
    
    if sweep_type == 'interval_skew':
        interval = int(param1_val)
        skew = int(param2_val * 1e9)  # Convert to nmol
        
        if x_idx is not None:
            cmd = f'./gradlew runDebug --args="sweep=interval x={x_idx} y={y_idx}"'
        else:
            cmd = f'./gradlew runDebug --args="interval={interval} skew={skew}"'
    
    else:  # dose_receptor
        dose = int(param1_val)
        baseline = 6.64e-19
        receptor_mult = param2_val / baseline
        
        if x_idx is not None:
            cmd = f'./gradlew runDebug --args="sweep=dose x={x_idx} y={y_idx}"'
        else:
            cmd = f'./gradlew runDebug --args="dose={dose} receptors={receptor_mult:.2f}"'
    
    return cmd

def show_summary(df, sweep_type, param1_name, param2_name):
    """Show summary statistics and interesting points"""
    
    print("\n" + "="*60)
    print("SWEEP SUMMARY")
    print("="*60)
    
    # Overall statistics
    total_sims = len(df)
    cure_count = (df['finalTumorCount'] < 10).sum()
    cure_rate = cure_count / total_sims
    
    print(f"\nTotal simulations: {total_sims}")
    print(f"Overall cure rate: {cure_rate:.2%} ({cure_count}/{total_sims})")
    
    # Best and worst parameter combinations
    grouped = df.groupby([param1_name, param2_name])['finalTumorCount'].agg(['mean', 'std', 'count'])
    grouped['cure_rate'] = df.groupby([param1_name, param2_name]).apply(
        lambda x: (x['finalTumorCount'] < 10).mean()
    ).values
    
    print("\n--- Best Parameter Combinations (highest cure rate) ---")
    best = grouped.nlargest(5, 'cure_rate')
    for idx, row in best.iterrows():
        param1, param2 = idx
        print(f"  {param1_name}={param1:.2e}, {param2_name}={param2:.2e}: " + 
              f"cure_rate={row['cure_rate']:.2f}, avg_final={row['mean']:.1f}")
    
    print("\n--- Most Variable Points (high stochastic variation) ---")
    variable = grouped.nlargest(5, 'std')
    for idx, row in variable.iterrows():
        param1, param2 = idx
        print(f"  {param1_name}={param1:.2e}, {param2_name}={param2:.2e}: " + 
              f"std={row['std']:.1f}, cure_rate={row['cure_rate']:.2f}")
    
    print("\n--- Interesting Transition Points (cure_rate ≈ 0.5) ---")
    transition = grouped.loc[(grouped['cure_rate'] > 0.3) & (grouped['cure_rate'] < 0.7)]
    for idx, row in transition.iterrows():
        param1, param2 = idx
        print(f"  {param1_name}={param1:.2e}, {param2_name}={param2:.2e}: " + 
              f"cure_rate={row['cure_rate']:.2f}")

def main():
    parser = argparse.ArgumentParser(
        description='Debug helper for sweep analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show summary and heatmap
  python debug_helper.py results/ParameterSweep_20260128_120000/
  
  # Interactive point selection
  python debug_helper.py results/ParameterSweep_20260128_120000/ --interactive
  
  # Generate command for specific point
  python debug_helper.py results/ParameterSweep_20260128_120000/ --x 4 --y 6
        """
    )
    
    parser.add_argument('sweep_dir', help='Path to sweep results directory')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='Interactive point selection from heatmap')
    parser.add_argument('--x', type=int, help='X index (param1) for direct command generation')
    parser.add_argument('--y', type=int, help='Y index (param2) for direct command generation')
    parser.add_argument('--summary-only', action='store_true',
                       help='Only show summary statistics, no plots')
    
    args = parser.parse_args()
    
    # Load data
    df, sweep_type, param1_name, param2_name = load_sweep_data(args.sweep_dir)
    
    # Show summary
    show_summary(df, sweep_type, param1_name, param2_name)
    
    # Get parameter arrays
    param1_vals = sorted(df[param1_name].unique())
    param2_vals = sorted(df[param2_name].unique())
    
    # Direct command generation
    if args.x is not None and args.y is not None:
        if args.x < len(param1_vals) and args.y < len(param2_vals):
            param1_val = param1_vals[args.x]
            param2_val = param2_vals[args.y]
            
            print(f"\n--- Point ({args.x}, {args.y}) ---")
            print(f"{param1_name}: {param1_val}")
            print(f"{param2_name}: {param2_val}")
            
            cmd = generate_debug_command(sweep_type, param1_val, param2_val, 
                                        args.x, args.y)
            print(f"\nDebug command:")
            print(f"  {cmd}\n")
        else:
            print(f"\nERROR: Indices out of bounds!")
            print(f"Valid x range: 0-{len(param1_vals)-1}")
            print(f"Valid y range: 0-{len(param2_vals)-1}\n")
        return
    
    # Interactive or static heatmap
    if not args.summary_only:
        if args.interactive:
            selected = interactive_selection(df, sweep_type, param1_name, param2_name)
            
            if selected:
                print("\n" + "="*60)
                print("SELECTED POINTS SUMMARY")
                print("="*60)
                for x, y, p1, p2, cmd in selected:
                    print(f"\nPoint ({x}, {y}): {param1_name}={p1}, {param2_name}={p2}")
                    print(f"  Command: {cmd}")
        else:
            fig, _, _ = create_heatmap(df, sweep_type, param1_name, param2_name)
            print("\nShowing heatmap... (close window to continue)")
            plt.show()
            
            print("\nTo debug a specific point, run:")
            print("  python debug_helper.py " + args.sweep_dir + " --x X_INDEX --y Y_INDEX")
            print("\nOr use interactive mode:")
            print("  python debug_helper.py " + args.sweep_dir + " --interactive\n")

if __name__ == '__main__':
    main()
