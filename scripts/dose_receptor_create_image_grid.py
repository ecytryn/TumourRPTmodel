#!/usr/bin/env python3
"""
Create comparison grid of tumor images across dose-receptor parameter sweep

This creates a PDF with images arranged in a grid:
- Rows: Different parameter values (doses or receptor densities)
- Columns: Time progression (days)

Usage:
    python create_image_grid_dose_receptor.py results/debug_sweeps/SWEEP_DIR/
    python create_image_grid_dose_receptor.py results/debug_sweeps/SWEEP_DIR/ --days 5 20 40 60
    python create_image_grid_dose_receptor.py results/debug_sweeps/SWEEP_DIR/ --replicate 1
"""

import argparse
import os
import sys
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

def create_summary_plot(sim_dir):
    """Create a summary plot of population and bound+internalized RL over time"""
    pop_file = sim_dir / 'populations.csv'
    pk_file = sim_dir / 'pkStateVariables.csv'
    
    if not pop_file.exists() or not pk_file.exists():
        return None
    
    try:
        # Load data (no headers)
        pop_data = np.loadtxt(pop_file, delimiter=',')
        pk_data = np.loadtxt(pk_file, delimiter=',')
        
        # Handle 1D case (single row)
        if len(pop_data.shape) == 1:
            pop_data = pop_data.reshape(1, -1)
        if len(pk_data.shape) == 1:
            pk_data = pk_data.reshape(1, -1)
        
        # Truncate to shorter length (they're logged slightly differently)
        min_len = min(pop_data.shape[0], pk_data.shape[0])
        pop_data = pop_data[:min_len, :]
        pk_data = pk_data[:min_len, :]
        
        # Time axis (hourly timesteps -> days)
        t_days = np.arange(min_len) / 24.0
        
        # Tumor population: NORMAL (col 1) + HYPOXIC (col 2)
        tumor_pop = pop_data[:, 1] + pop_data[:, 2]
        
        # Captive RL: NbHot (col 6) + NicHot (col 8)
        captive_rl = (pk_data[:, 6] + pk_data[:, 8]) * 1e9  # Convert to nmol
        
        # Create plot (bitmap, small size)
        fig, ax1 = plt.subplots(figsize=(3, 2))
        
        # Tumor population (blue, left axis)
        color1 = 'tab:blue'
        ax1.plot(t_days, tumor_pop, color=color1, linewidth=1.5, label='Tumor cells')
        ax1.set_xlabel('Days', fontsize=8)
        ax1.set_ylabel('Tumor cells', fontsize=8, color=color1)
        ax1.tick_params(axis='y', labelcolor=color1, labelsize=7)
        ax1.tick_params(axis='x', labelsize=7)
        ax1.set_ylim(bottom=0)
        
        # Captive RL (orange, right axis)
        ax2 = ax1.twinx()
        color2 = 'tab:orange'
        ax2.plot(t_days, captive_rl, color=color2, linewidth=1.5, label='Captive RL')
        ax2.set_ylabel('Bound+Int. RL (nmol)', fontsize=8, color=color2)
        ax2.tick_params(axis='y', labelcolor=color2, labelsize=7)
        ax2.set_ylim(bottom=0)
        
        plt.tight_layout()
        
        # Convert figure to image array (works on all platforms)
        fig.canvas.draw()
        
        # Use buffer_rgba() instead of tostring_rgb() for compatibility
        buf = fig.canvas.buffer_rgba()
        img = np.asarray(buf)
        
        # Convert RGBA to RGB (drop alpha channel)
        img = img[:, :, :3]
        
        plt.close(fig)
        
        return img
        
    except Exception as e:
        print(f"    Warning: Could not create plot for {sim_dir.name}: {e}")
        return None
        
        
def find_sweep_directories(sweep_dir):
    """Find all dose_X_rep_Y or receptor_X_rep_Y directories"""
    sweep_path = Path(sweep_dir)
    
    # Try both patterns
    dose_pattern = re.compile(r'dose_(\d+(?:\.\d+)?)_rep_(\d+)')
    receptor_pattern = re.compile(r'receptor_(\d+(?:\.\d+)?)_rep_(\d+)')
    
    dirs = {}
    sweep_type = None
    
    for item in sweep_path.iterdir():
        if item.is_dir():
            # Try dose pattern first
            match = dose_pattern.match(item.name)
            if match:
                param_value = float(match.group(1))
                rep = int(match.group(2))
                sweep_type = 'dose'
            else:
                # Try receptor pattern
                match = receptor_pattern.match(item.name)
                if match:
                    param_value = float(match.group(1))
                    rep = int(match.group(2))
                    sweep_type = 'receptor'
            
            if match:
                if param_value not in dirs:
                    dirs[param_value] = {}
                dirs[param_value][rep] = item
    
    return dirs, sweep_type

def find_image_for_day(sim_dir, day):
    """Find the tumor image for a specific day"""
    image_dir = sim_dir / 'tumour_images'
    
    if not image_dir.exists():
        return None
    
    # Look for day_XXX.png
    image_file = image_dir / f'day_{day:03d}.png'
    
    if image_file.exists():
        return image_file
    
    # If exact day not found, find closest
    all_images = sorted(image_dir.glob('day_*.png'))
    if not all_images:
        return None
    
    # Extract day numbers
    day_numbers = []
    for img in all_images:
        match = re.search(r'day_(\d+)', img.name)
        if match:
            day_numbers.append((int(match.group(1)), img))
    
    if not day_numbers:
        return None
    
    # Find closest
    closest = min(day_numbers, key=lambda x: abs(x[0] - day))
    return closest[1]

def get_outcome(sim_dir):
    """Determine if simulation resulted in cure or failure"""
    pop_file = sim_dir / 'populations.csv'
    
    if not pop_file.exists():
        return 'UNKNOWN'
    
    try:
        # Read last line
        with open(pop_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                return 'UNKNOWN'
            
            last_line = lines[-1].strip()
            parts = last_line.split(',')
            
            if len(parts) < 2:
                return 'UNKNOWN'
            
            final_pop = int(parts[1])  # Assuming second column is normoxic count
            
            return 'CURE' if final_pop < 10 else 'FAILURE'
    except:
        return 'UNKNOWN'

def create_grid_pdf(sweep_dir, output_file, days=None, replicate=None):
    """Create PDF with image comparison grid"""
    
    # Find all directories
    dirs, sweep_type = find_sweep_directories(sweep_dir)
    
    if not dirs:
        print(f"ERROR: No simulation directories found in {sweep_dir}")
        sys.exit(1)
    
    if sweep_type is None:
        print(f"ERROR: Could not determine sweep type (dose or receptor)")
        sys.exit(1)
    
    # Get parameter values and replicates
    param_values = sorted(dirs.keys())
    all_reps = set()
    for param_dict in dirs.values():
        all_reps.update(param_dict.keys())
    replicates = sorted(all_reps)
    
    param_label = 'Dose (nmol)' if sweep_type == 'dose' else 'Receptor density'
    
    print(f"Found {len(param_values)} {sweep_type} values: {param_values}")
    print(f"Found {len(replicates)} replicates: {replicates}")
    
    # If specific replicate requested, filter
    if replicate is not None:
        if replicate not in replicates:
            print(f"ERROR: Replicate {replicate} not found")
            sys.exit(1)
        replicates = [replicate]
        print(f"Using only replicate {replicate}")
    
    # Check for incomplete data and warn
    complete_reps = []
    for rep in replicates:
        complete = True
        for param_value in param_values:
            if param_value not in dirs or rep not in dirs[param_value]:
                complete = False
                break
        if complete:
            complete_reps.append(rep)
        else:
            print(f"WARNING: Replicate {rep} incomplete - will skip missing data")
    
    if complete_reps:
        print(f"Complete replicates: {complete_reps}")
    
    # Default days if not specified
    if days is None:
        # Try to infer from first directory
        first_param = param_values[0]
        first_rep = replicates[0]
        first_dir = dirs[first_param][first_rep]
        image_dir = first_dir / 'tumour_images'
        
        if image_dir.exists():
            all_images = sorted(image_dir.glob('day_*.png'))
            day_numbers = []
            for img in all_images:
                match = re.search(r'day_(\d+)', img.name)
                if match:
                    day_numbers.append(int(match.group(1)))
            
            if day_numbers:
                # Select ~6-8 evenly spaced timepoints
                n_timepoints = min(8, len(day_numbers))
                indices = np.linspace(0, len(day_numbers)-1, n_timepoints, dtype=int)
                days = [day_numbers[i] for i in indices]
            else:
                days = [0, 5, 10, 20, 40, 60, 80, 100]
        else:
            days = [0, 5, 10, 20, 40, 60, 80, 100]
    
    print(f"Using days: {days}")
    
    # Create PDF
    with PdfPages(output_file) as pdf:
        
        for rep in replicates:
            print(f"\nGenerating grid for replicate {rep}...")
            
            # Create figure
            n_rows = len(param_values)
            n_cols = len(days) + 1  # Extra column for summary plot
            
            fig, axes = plt.subplots(n_rows, n_cols, 
                                    figsize=(2*n_cols, 2*n_rows),
                                    squeeze=False)
            
            sweep_title = f'{param_label} Sweep, Replicate {rep}'
            fig.suptitle(sweep_title, fontsize=16, fontweight='bold')
            
            # Plot each cell in grid
            for i, param_value in enumerate(param_values):
                
                if param_value not in dirs or rep not in dirs[param_value]:
                    print(f"  WARNING: Missing data for {sweep_type}={param_value}, rep={rep}")
                    continue
                
                sim_dir = dirs[param_value][rep]
                outcome = get_outcome(sim_dir)
                
                for j, day in enumerate(days):
                    ax = axes[i, j]
                    
                    image_path = find_image_for_day(sim_dir, day)
                    
                    if image_path and image_path.exists():
                        try:
                            img = Image.open(image_path)
                            ax.imshow(img)
                            ax.axis('off')
                        except Exception as e:
                            ax.text(0.5, 0.5, 'Error\nloading\nimage',
                                  ha='center', va='center', fontsize=8)
                            ax.axis('off')
                    else:
                        ax.text(0.5, 0.5, 'No\nimage',
                              ha='center', va='center', fontsize=8)
                        ax.axis('off')
                    
                    # Column labels (days) - top row only
                    if i == 0:
                        ax.set_title(f'Day {day}', fontsize=10)
                    
                    # Row labels (parameter values) - first column only
                    if j == 0:
                        color = 'green' if outcome == 'CURE' else 'red'
                        if sweep_type == 'dose':
                            label_text = f'Dose={param_value:.0f}\n{outcome}'
                        else:
                            label_text = f'Receptor={param_value:.2f}\n{outcome}'
                        
                        ax.text(-0.12, 0.5, label_text,
                               transform=ax.transAxes,
                               ha='right', va='center',
                               fontsize=10, fontweight='bold',
                               color=color)

                ax = axes[i, n_cols-1]  # Last column
                
                plot_img = create_summary_plot(sim_dir)
                
                if plot_img is not None:
                    ax.imshow(plot_img)
                    ax.axis('off')
                else:
                    ax.text(0.5, 0.5, 'No\ndata',
                          ha='center', va='center', fontsize=8)
                    ax.axis('off')
                
                # Label this column (top row only)
                if i == 0:
                    ax.set_title('Population\n& Dose', fontsize=10)
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close()
            
            print(f"  Page added for replicate {rep}")
    
    print(f"\nPDF saved to: {output_file}")
    print(f"Pages: {len(replicates)} (one per replicate)")

def create_outcome_summary(sweep_dir, output_file):
    """Create a simple outcome summary table"""
    
    dirs, sweep_type = find_sweep_directories(sweep_dir)
    param_values = sorted(dirs.keys())
    
    all_reps = set()
    for param_dict in dirs.values():
        all_reps.update(param_dict.keys())
    replicates = sorted(all_reps)
    
    param_label = 'Dose' if sweep_type == 'dose' else 'Receptor'
    
    with open(output_file, 'w') as f:
        f.write("Outcome Summary\n")
        f.write("=" * 60 + "\n\n")
        
        # Header
        f.write(f"{param_label:<12}")
        for rep in replicates:
            f.write(f"Rep{rep:<8}")
        f.write(f"{'CureRate':<10}\n")
        f.write("-" * 60 + "\n")
        
        # Data
        for param_value in param_values:
            if sweep_type == 'dose':
                f.write(f"{param_value:<12.0f}")
            else:
                f.write(f"{param_value:<12.2f}")
            
            outcomes = []
            for rep in replicates:
                if param_value in dirs and rep in dirs[param_value]:
                    outcome = get_outcome(dirs[param_value][rep])
                    f.write(f"{outcome:<9}")
                    outcomes.append(1 if outcome == 'CURE' else 0)
                else:
                    f.write(f"{'MISSING':<9}")
            
            if outcomes:
                cure_rate = sum(outcomes) / len(outcomes)
                f.write(f"{cure_rate:.1%}\n")
            else:
                f.write("N/A\n")
    
    print(f"Outcome summary saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Create image comparison grid from dose-receptor debug sweep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create grid with all replicates, auto-detected days
  python create_image_grid_dose_receptor.py results/debug_sweeps/dose_row_receptor1.0_20260129/
  
  # Specify specific days
  python create_image_grid_dose_receptor.py results/debug_sweeps/dose_row_receptor1.0_20260129/ --days 5 20 40 60
  
  # Only show replicate 1
  python create_image_grid_dose_receptor.py results/debug_sweeps/dose_row_receptor1.0_20260129/ --replicate 1
        """
    )
    
    parser.add_argument('sweep_dir', help='Path to sweep directory')
    parser.add_argument('--days', type=int, nargs='+', 
                       help='Specific days to show (default: auto-detect)')
    parser.add_argument('--replicate', type=int,
                       help='Show only specific replicate (default: all)')
    parser.add_argument('--output', '-o',
                       help='Output PDF filename (default: comparison_grid.pdf in sweep dir)')
    
    args = parser.parse_args()
    
    # Validate sweep directory
    sweep_path = Path(args.sweep_dir)
    if not sweep_path.exists():
        print(f"ERROR: Directory not found: {args.sweep_dir}")
        sys.exit(1)
    
    # Set output filename
    if args.output:
        output_pdf = args.output
    else:
        output_pdf = sweep_path / 'comparison_grid.pdf'
    
    # Create grid PDF
    create_grid_pdf(args.sweep_dir, output_pdf, args.days, args.replicate)
    
    # Create outcome summary
    summary_file = sweep_path / 'outcome_summary.txt'
    create_outcome_summary(args.sweep_dir, summary_file)
    
    print("\n" + "="*60)
    print("COMPLETE")
    print("="*60)
    print(f"PDF: {output_pdf}")
    print(f"Summary: {summary_file}")

if __name__ == '__main__':
    main()
