#!/usr/bin/env python3
"""
Visualize RadioBio Survival Probability Lookup Tables
Generates heatmaps showing SF vs age for debugging
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def plot_lookup_table(data_file, output_file):
    """
    Plot survival probability lookup table as heatmap
    
    Args:
        data_file: CSV file with columns: age, SF_norm_orig, SF_hypo_orig, SF_norm_opt, SF_hypo_opt
        output_file: Path to save PNG image
    """
    
    # Read data
    try:
        data = np.loadtxt(data_file, delimiter=',', skiprows=1)
    except Exception as e:
        print(f"Error reading {data_file}: {e}")
        return
    
    if len(data) == 0:
        print(f"No data in {data_file}")
        return
    
    # Extract columns
    ages = data[:, 0]
    SF_norm_orig = data[:, 1]
    SF_hypo_orig = data[:, 2]
    SF_norm_opt = data[:, 3]
    SF_hypo_opt = data[:, 4]
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Survival Probability Lookup Table', fontsize=16, fontweight='bold')
    
    # Plot 1: Original Normoxic
    ax = axes[0, 0]
    ax.plot(ages, SF_norm_orig, 'b-', linewidth=2, label='Normoxic')
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Original Method - Normoxic', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 2: Original Hypoxic
    ax = axes[0, 1]
    ax.plot(ages, SF_hypo_orig, 'r-', linewidth=2, label='Hypoxic')
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Original Method - Hypoxic', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 3: Optimized Normoxic
    ax = axes[1, 0]
    ax.plot(ages, SF_norm_opt, 'b-', linewidth=2, label='Normoxic')
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Optimized Method - Normoxic', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 4: Optimized Hypoxic
    ax = axes[1, 1]
    ax.plot(ages, SF_hypo_opt, 'r-', linewidth=2, label='Hypoxic')
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Optimized Method - Hypoxic', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {output_file}")
    plt.close()


def plot_comparison(data_file, output_file):
    """
    Plot comparison of Original vs Optimized methods
    
    Args:
        data_file: CSV file with columns: age, SF_norm_orig, SF_hypo_orig, SF_norm_opt, SF_hypo_opt
        output_file: Path to save PNG image
    """
    
    # Read data
    try:
        data = np.loadtxt(data_file, delimiter=',', skiprows=1)
    except Exception as e:
        print(f"Error reading {data_file}: {e}")
        return
    
    if len(data) == 0:
        print(f"No data in {data_file}")
        return
    
    # Extract columns
    ages = data[:, 0]
    SF_norm_orig = data[:, 1]
    SF_hypo_orig = data[:, 2]
    SF_norm_opt = data[:, 3]
    SF_hypo_opt = data[:, 4]
    
    # Calculate differences
    diff_norm = np.abs(SF_norm_orig - SF_norm_opt)
    diff_hypo = np.abs(SF_hypo_orig - SF_hypo_opt)
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 12))
    fig.suptitle('Original vs Optimized Comparison', fontsize=16, fontweight='bold')
    
    # Plot 1: Normoxic comparison
    ax = axes[0]
    ax.plot(ages, SF_norm_orig, 'b-', linewidth=2, label='Original', alpha=0.7)
    ax.plot(ages, SF_norm_opt, 'r--', linewidth=2, label='Optimized', alpha=0.7)
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Normoxic Cells', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 2: Hypoxic comparison
    ax = axes[1]
    ax.plot(ages, SF_hypo_orig, 'b-', linewidth=2, label='Original', alpha=0.7)
    ax.plot(ages, SF_hypo_opt, 'r--', linewidth=2, label='Optimized', alpha=0.7)
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Survival Fraction', fontsize=12)
    ax.set_title('Hypoxic Cells', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 3: Absolute differences
    ax = axes[2]
    ax.semilogy(ages, diff_norm, 'b-', linewidth=2, label='Normoxic', alpha=0.7)
    ax.semilogy(ages, diff_hypo, 'r-', linewidth=2, label='Hypoxic', alpha=0.7)
    ax.axhline(y=1e-6, color='k', linestyle='--', alpha=0.5, label='1e-6 threshold')
    ax.set_xlabel('Age (days)', fontsize=12)
    ax.set_ylabel('Absolute Difference', fontsize=12)
    ax.set_title('Difference: |Original - Optimized|', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved comparison plot to {output_file}")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_lookup_table.py <data_file.csv> [output_prefix]")
        print("Generates two plots:")
        print("  1. <output_prefix>_lookup.png - 4-panel lookup table")
        print("  2. <output_prefix>_comparison.png - comparison and differences")
        sys.exit(1)
    
    data_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_prefix = sys.argv[2]
    else:
        # Use data file name without extension
        output_prefix = os.path.splitext(data_file)[0]
    
    # Generate plots
    plot_lookup_table(data_file, f"{output_prefix}_lookup.png")
    plot_comparison(data_file, f"{output_prefix}_comparison.png")
    
    print("\nDone! Generated:")
    print(f"  - {output_prefix}_lookup.png")
    print(f"  - {output_prefix}_comparison.png")
