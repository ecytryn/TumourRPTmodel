#!/usr/bin/env python3
"""
Plot analytical solutions to the simplified PK model

Based on the supplemental material equations:
- C_cen(t) = C_cen0 * exp(-lambda_bio * t)
- C_cen^H(t) = C_cen0^H * exp(-(lambda_bio + lambda_decay) * t)
- C_ic(t) = k_int * R_T_tilde * exp(-k_rel * t) * integral[...]
- C_ic^H(t) = k_int * R_T_tilde * exp(-k_rel * t) * exp(-lambda_decay * t) * integral[...]
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# =============================================================================
# PARAMETERS from Table in supplemental material
# =============================================================================

# Rate constants (all but k_on in sec^-1)
lambda_bio = 1.6e-4 / 60      # Clearance from central compartment
lambda_decay = 7.14e-5 / 60   # Lu-177 decay
k_on = 0.046 * 1e6 / 60.0     # Binding rate (m^3/(mol*s))
k_off = 0.368 / 60            # Unbinding rate <--- Arman's value
k_int = 0.001 / 60            # Internalization rate
k_rel = 2e-4 / 60             # Release rate

# Derived parameter
beta = (k_off + k_int) / k_on  # Dissociation constant-like parameter (mol/m^3)

# =============================================================================
# SIMULATION GEOMETRY - MATCH THESE TO JAVA OUTPUT
# =============================================================================

# **COPY THESE VALUES FROM JAVA CONSOLE OUTPUT**
tumour_cells_2D = 3285      # Number of cells in 2D cross-section (from Grid.countTumorCells())
numberOfVessels = 365       # Number of vessels near tumor (from Grid.countVesselsNearTumor())

CELL_LENGTH = 1e-5          # m (10 μm)
INTERSTITIAL_FRACTION = 0.4 # Fraction of tumor that is extracellular

# =============================================================================
# DERIVE VOLUMES (matching Java SimParams.computeTumorVolume logic)
# =============================================================================

# Calculate tumor geometry (matching Java)
radius_cells = np.sqrt(tumour_cells_2D / np.pi)  # radius in cell lengths
radius_m = radius_cells * CELL_LENGTH             # radius in meters
height_m = 2.0 * radius_m                          # cylindrical extrusion

# Tumor volume (cylindrical)
tumourVolume_m3 = np.pi * radius_m**2 * height_m   # m³
tumourVolume_L = tumourVolume_m3 * 1000            # L (for display)

# Extracellular/interstitial volume
V_ec_m3 = INTERSTITIAL_FRACTION * tumourVolume_m3  # m³
V_ec_L = V_ec_m3 * 1000                            # L

# Vascular volume (each vessel is a cylinder: cross-section × height)
vesselVolumeInTumor_m3 = (CELL_LENGTH**2) * height_m  # m³ per vessel
V_vasc_m3 = numberOfVessels * vesselVolumeInTumor_m3  # m³
V_vasc_L = V_vasc_m3 * 1000                            # L

# Central compartment (body blood volume - fixed)
V_cen_m3 = 0.458e-3  # m³

# Flow/permeability rates (L/min)  - these drop out of the reduced model
#F_T = 2.9e-3*60
#PS = 3.5e-3*60

# Tumour sizes:
# N=number of cells. V = 1e-12 N L. r = 10 (N/4)^(1/3) um
# N=1e11. V = 1e-12 1e11 L = 100 mL. r = 10 (1e11/4)^(1/3) um = 30 mm = 3cm
# N=3.2e7. V = 0.033 mL. r = 2 mm.
# N=4e6. V = 4 uL. r = 1 mm.

# =============================================================================
# RECEPTOR CALCULATIONS
# =============================================================================

AVOGADRO = 6.022e23
RECEPTORS_PER_CELL = 5e5

# Total cells in 3D (with cylindrical extrusion)
# Matches: totalCells3D = tumorCells * height (in Java)
height_cells = height_m / CELL_LENGTH
totalCells3D = tumour_cells_2D * height_cells

# Total receptors
R_total_mol = totalCells3D * RECEPTORS_PER_CELL / AVOGADRO  # mol
R_total_nmol = R_total_mol * 1e9                             # nmol

# Receptor concentration scaled by V_ec (because it always appears in that form in the ODE)
R_T_tilde_mol_m3 = R_total_mol / V_ec_m3      # mol/m³
R_T_tilde_nmol_L = R_T_tilde_mol_m3 * 1e9 / 1000  # nmol/L for the figure

# Initial conditions
dose_mol = 100e-9       # Total injection (mol)
dose_nmol = dose_mol*1e-9       # Total injection (nmol) for the figure title
hot_fraction = 0.1       # Fraction that is radioactive
C_cen0 = dose_mol / V_cen_m3          # Initial central concentration (mol/m^3)
C_cen0_H = hot_fraction * C_cen0    # Initial hot concentration (mol/m^3)

gamma = beta / C_cen0    # Dimensionless parameter

# =============================================================================
# DIAGNOSTIC OUTPUT - Compare with Java
# =============================================================================

print("\n" + "="*70)
print("GEOMETRY AND PARAMETER COMPARISON")
print("="*70)

print("\nInput from Java:")
print(f"  Tumor cells (2D): {tumour_cells_2D}")
print(f"  Vessels near tumor: {numberOfVessels}")

print("\nCalculated volumes:")
print(f"  V_ec = {V_ec_m3:.6e} m³")
print(f"  V_v = {V_vasc_m3:.6e} m³")
print(f"  V_cen = {V_cen_m3:.6e} m³")
print(f"  Tumor volume = {tumourVolume_m3:.6e} m³")

print("\nReceptor calculations:")
print(f"  Total cells (3D): {totalCells3D:.0f}")
print(f"  R_total = {R_total_mol:.6e} mol")
print(f"  R_T_tilde = {R_T_tilde_mol_m3:.6e} mol/m³")

print("\nDerived parameters:")
print(f"  beta = {beta:.6e} mol/m³")
print(f"  gamma = {gamma:.6e}")

print("\nInitial conditions:")
print(f"  C_cen0 = {C_cen0:.6e} mol/m³")
print(f"  C_cen0_H = {C_cen0_H:.6e} mol/m³")

print("\n" + "="*70)
print("Copy the values above and compare with Java output")
print("="*70 + "\n")

# =============================================================================
# TIME ARRAY
# =============================================================================

t_max_days = 30                 # Simulate for 30 days
t_max_sec = t_max_days * 24 * 60 * 60   # Convert to seconds
dt = 3600                       # Time step (hours)
t_sec = np.arange(0, t_max_sec, dt)
t_days = t_sec / (24*3600)      # For plotting in days

# =============================================================================
# ANALYTICAL SOLUTIONS
# =============================================================================

# Central compartment concentrations
C_cen = C_cen0 * np.exp(-lambda_bio * t_sec)
C_cen_H = C_cen0_H * np.exp(-(lambda_bio + lambda_decay) * t_sec)

# Bound concentrations (from QSS)
C_b = R_T_tilde_mol_m3 * C_cen / (C_cen + beta)
C_b_H = R_T_tilde_mol_m3 * C_cen_H / (C_cen + beta)

# Intracellular concentrations - need to evaluate the integral numerically
# C_ic = k_int * R_T_tilde * exp(-k_rel * t) * integral[exp(k_rel*s) * exp(-lambda_bio*s) / (exp(-lambda_bio*s) + gamma) ds from 0 to t]

# Define the integrand
def integrand(s):
    return np.exp(k_rel * s) * np.exp(-lambda_bio * s) / (np.exp(-lambda_bio * s) + gamma)

# Compute the integral at each time point using cumulative integration
integral_values = cumulative_trapezoid(integrand(t_sec), t_sec, initial=0)

# Calculate C_ic and C_ic^H
C_ic = k_int * R_T_tilde_mol_m3 * np.exp(-k_rel * t_sec) * integral_values
C_ic_H = np.exp(-lambda_decay * t_sec) * C_ic

# =============================================================================
# CONVERT TO TOTAL AMOUNTS (nmol)
# =============================================================================

# Central compartment
A_cen = C_cen * V_cen_m3
A_cen_H = C_cen_H * V_cen_m3

# Vascular compartment
A_vasc = C_cen * V_vasc_m3
A_vasc_H = C_cen_H * V_vasc_m3

# Interstitial / extracellular compartment
A_ec = C_cen * V_ec_m3
A_ec_H = C_cen_H * V_ec_m3

# Bound (use V_ec as per the model)
A_b = C_b * V_ec_m3
A_b_H = C_b_H * V_ec_m3

# Intracellular
A_ic = C_ic * V_ec_m3
A_ic_H = C_ic_H * V_ec_m3

A_tumour_H = A_vasc_H + A_ec_H + A_b_H + A_ic_H
A_captive_H = A_b_H + A_ic_H
A_free = A_vasc_H + A_ec_H

# =============================================================================
# PLOTTING
# =============================================================================

fig, axes = plt.subplots(3, 1, figsize=(8, 12))

# Upper plot: Central vs Tumor
ax1 = axes[0]
ax1.plot(t_days, A_cen_H * 1e9, linewidth=2.5, label='Central')
ax1.plot(t_days, A_tumour_H * 1e9, linewidth=2.5, label='Tumor (total)')
ax1.set_xlabel('Time (days)', fontsize=14)
ax1.set_ylabel('Amount (nmol)', fontsize=14)
ax1.set_title('Radiation Focusing: Central vs Tumor', fontsize=15, fontweight='bold')
ax1.legend(fontsize=12, loc='best')
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# Middle plot: Captive (bound+internalized) vs Free (vascular+extracellular)
ax2 = axes[1]
ax2.plot(t_days, A_captive_H * 1e9, '#2ca02c', linewidth=2.5, label='Captive (bound + internalized)')
ax2.plot(t_days, A_free * 1e9, '#d62728', linewidth=2.5, label='Free (vascular + extracellular)')
ax2.set_xlabel('Time (days)', fontsize=14)
ax2.set_ylabel('Amount (nmol)', fontsize=14)
ax2.set_title('Ligand Targeting: Captive vs Free', fontsize=15, fontweight='bold')
ax2.legend(fontsize=12, loc='best')
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

# Lower plot: bound, internalized, vascular, extracellular
ax3 = axes[2]
#ax3.plot(t_days, A_cen_H * 1e9, linewidth=2.5, label='Central')
ax3.plot(t_days, A_vasc_H * 1e9, linewidth=2.5, label='Vascular')
ax3.plot(t_days, A_ec_H * 1e9, linewidth=2.5, label='Extracellular')
ax3.plot(t_days, A_b_H * 1e9, linewidth=2.5, label='Bound')
ax3.plot(t_days, A_ic_H * 1e9, linewidth=2.5, label='Intracellular')
ax3.set_xlabel('Time (days)', fontsize=14)
ax3.set_ylabel('Amount (nmol)', fontsize=14)
ax3.set_title('Tumour Compartments', fontsize=15, fontweight='bold')
ax3.legend(fontsize=12, loc='best')
ax3.grid(True, alpha=0.3)
ax3.set_yscale('log')

plt.suptitle(f'Reduced Model (PK)\n' + 
             f'{tumour_cells_2D} cells, R̃_T = {R_T_tilde_nmol_L:.2e} nmol/L, ' +
             f'V_ec = {V_ec_L:.2e} L',
             fontsize=13, fontweight='bold')

fig.tight_layout()

# Save figure
output_file = 'results/compare_models/pk_reduced_model.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_file}")

plt.show()
