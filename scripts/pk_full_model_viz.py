#!/usr/bin/env python3
"""
Visualize pharmacokinetic compartment data over time

Plots either concentrations or total amounts of hot radioligand in:
- Central compartment
- Tumor compartments (vasculature + interstitium + bound + internalized)

Reads state_variable.csv from a specified parameter sweep run.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Plot mode: 'concentration' or 'total_amount'
PLOT_MODE = 'concentration'  # Change to 'concentration' to plot concentrations

# =============================================================================
# COMPARTMENT VOLUMES (in liters)  - these need to be pulled from simulation output data
# =============================================================================

V_cen = 0.458                   # Central compartment
V_v_Domain = 3.885e-6          # Tumor vasculature (9713 vessels × 4e-10 L/vessel)
V_int_Domain = 2.56e-5         # Tumor interstitium (0.4 × domain volume)

# =============================================================================
# STATE VARIABLE INDICES (from state_variable.csv)
# =============================================================================

C_cen_hot_index = 0
C_cen_cold_index = 1
C_v_hot_index = 2
C_v_cold_index = 3
C_int_hot_index = 4
C_int_cold_index = 5
C_b_hot_index = 6
C_b_cold_index = 7
C_intern_hot_index = 8
C_intern_cold_index = 9

# =============================================================================
# MAIN CODE
# =============================================================================

def load_pk_data():
    """Load PK state variable data"""
    
    pk_file = "results/single_runs/pkStateVariables.csv"
    
    if not Path(pk_file).exists():
        raise FileNotFoundError(f"PK data file not found: {pk_file}")
    
    # Load data (each row is one time point, columns are state variables)
    data = np.loadtxt(pk_file, delimiter=',')
    
    print(f"Loaded PK data from: {pk_file}")
    print(f"  Shape: {data.shape} (time points × state variables)")
    print(f"  Duration: {len(data)} hours = {len(data)/24:.1f} days")
    
    return data

def calculate_amounts(data):
    """Calculate total drug amounts from concentrations"""
    
    # Extract hot concentrations (only plotting hot drug)
    C_cen = data[:, C_cen_hot_index]
    C_v = data[:, C_v_hot_index]
    C_int = data[:, C_int_hot_index]
    C_b = data[:, C_b_hot_index]
    C_intern = data[:, C_intern_hot_index]
    
    # Calculate amounts (nmol)
    A_cen = C_cen * V_cen
    A_v = C_v * V_v_Domain
    A_int = C_int * V_int_Domain
    A_b = C_b * V_int_Domain
    A_intern = C_intern * V_int_Domain
    
    # Total tumor amount (all tumor compartments)
    A_tumor = A_v + A_int + A_b + A_intern
    
    return {
        'central': A_cen,
        'vasculature': A_v,
        'interstitium': A_int,
        'bound': A_b,
        'internalized': A_intern,
        'tumor_total': A_tumor
    }

def calculate_concentrations(data):
    """Extract concentrations directly"""
    
    # Extract hot concentrations
    C_cen = data[:, C_cen_hot_index]
    C_v = data[:, C_v_hot_index]
    C_int = data[:, C_int_hot_index]
    C_b = data[:, C_b_hot_index]
    C_intern = data[:, C_intern_hot_index]
    
    # For plotting, sum tumor concentrations (note: not physically meaningful
    # to add concentrations with different volumes, but useful for comparison)
    # Better: show them separately or as fractions
    
    return {
        'central': C_cen,
        'vasculature': C_v,
        'interstitium': C_int,
        'bound': C_b,
        'internalized': C_intern,
    }

def plot_pk_data(data, mode='total_amount'):
    """Create PK compartment plots"""
    
    # Time axis (each row is 1 hour)
    time_hours = np.arange(len(data))
    time_days = time_hours / 24.0
    
    if mode == 'total_amount':
        amounts = calculate_amounts(data)
        
        # Create figure with two subplots
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # Plot 1: Central vs Total Tumor
        ax1 = axes[0]
        ax1.plot(time_days, amounts['central'], 'r-', linewidth=2, label='Central')
        ax1.plot(time_days, amounts['tumor_total'], 'b-', linewidth=2, label='Tumor (total)')
        ax1.set_xlabel('Time (days)', fontsize=12)
        ax1.set_ylabel('Amount (nmol)', fontsize=12)
        ax1.set_title(f'Hot Radioligand Distribution: Central vs Tumor\n', fontsize=14)
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, max(time_days))
        ax1.set_yscale('log')

        # Plot 2: Tumor compartments breakdown
        ax2 = axes[1]
        ax2.plot(time_days, amounts['vasculature'], linewidth=2, label='Vasculature')
        ax2.plot(time_days, amounts['interstitium'], linewidth=2, label='Extracellular')
        ax2.plot(time_days, amounts['bound'], linewidth=2, label='Bound')
        ax2.plot(time_days, amounts['internalized'], linewidth=2, label='Inrtracellular')
        ax2.set_xlabel('Time (days)', fontsize=12)
        ax2.set_ylabel('Amount (nmol)', fontsize=12)
        ax2.set_title('Tumor Compartments Breakdown', fontsize=14)
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, max(time_days))
        ax2.set_yscale('log')
        
        # Print summary statistics
        print("\n=== DRUG DISTRIBUTION SUMMARY ===")
        print(f"Peak central amount: {amounts['central'].max():.2f} nmol")
        print(f"Peak total tumor amount: {amounts['tumor_total'].max():.2f} nmol")
        print(f"\nTumor compartment peaks:")
        print(f"  Vasculature: {amounts['vasculature'].max():.2e} nmol")
        print(f"  Interstitium: {amounts['interstitium'].max():.2e} nmol")
        print(f"  Bound: {amounts['bound'].max():.2e} nmol")
        print(f"  Internalized: {amounts['internalized'].max():.2e} nmol")
        
        # Calculate what fraction of total tumor drug is bound+internalized
        tumor_peak_idx = np.argmax(amounts['tumor_total'])
        bound_intern_fraction = (amounts['bound'][tumor_peak_idx] + 
                                amounts['internalized'][tumor_peak_idx]) / amounts['tumor_total'][tumor_peak_idx]
        print(f"\nAt peak tumor loading:")
        print(f"  Bound + Internalized fraction: {bound_intern_fraction:.1%}")
        
    else:  # concentration mode
        concs = calculate_concentrations(data)
        
        # Create figure with two subplots
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # Plot 1: Central concentration
        ax1 = axes[0]
        ax1.plot(time_days, concs['central'], 'r-', linewidth=2)
        ax1.set_xlabel('Time (days)', fontsize=12)
        ax1.set_ylabel('Concentration (nmol/L)', fontsize=12)
        ax1.set_title(f'Central Concentration\n', fontsize=14)
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, max(time_days))
        
        # Plot 2: Tumor compartment concentrations
        ax2 = axes[1]
        ax2.plot(time_days, concs['vasculature'], linewidth=2, label='Vasculature')
        ax2.plot(time_days, concs['interstitium'], linewidth=2, label='Extracellular')
        ax2.plot(time_days, concs['bound'], linewidth=2, label='Bound')
        ax2.plot(time_days, concs['internalized'], linewidth=2, label='Intracellular')
        ax2.set_xlabel('Time (days)', fontsize=12)
        ax2.set_ylabel('Concentration (nmol/L)', fontsize=12)
        ax2.set_title('Tumor Compartment Concentrations', fontsize=14)
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, max(time_days))
    
    plt.tight_layout()
    
    # Save figure
    output_file = "results/single_runs/pkStateVariables.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_file}")
    
    return fig

def main():
    """Main execution"""
    
    print("=== PK COMPARTMENT VISUALIZATION ===")
    print(f"Plot mode: {PLOT_MODE}")
    print()
    
    # Load data
    data = load_pk_data()
    
    # Create plots
    plot_pk_data(data, mode=PLOT_MODE)
    
    print("\n=== VISUALIZATION COMPLETE ===")

if __name__ == "__main__":
    main()
