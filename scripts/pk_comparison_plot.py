#!/usr/bin/env python3
"""
Compare three PK model implementations:
1. Java QSS (from simulation CSV)
2. Python reduced model (analytical QSS)
3. Python full ODE (no approximations)

Creates figure showing N_cen_hot and N_ic_hot from all three models.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import matplotlib as mpl

mpl.rcParams.update({
    "figure.figsize": (3.35, 2.4),   # single column
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,   # editable text in Illustrator
    "ps.fonttype": 42
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# **USER: SET THESE PATHS**
# Path to Java simulation output directory
java_output_dir = "results/single_runs"  # Will find most recent CustomRunToMakeAFig

# Cell count (from Java console output "Tumor cells (2D): ")
TUMOR_CELLS_2D = 3285  # **UPDATE THIS FROM JAVA OUTPUT**

# =============================================================================
# FIND JAVA OUTPUT
# =============================================================================

print("\n" + "="*70)
print("PK MODEL COMPARISON")
print("="*70)

# Find most recent CustomRunToMakeAFig directory
if os.path.exists(java_output_dir):
    subdirs = [d for d in os.listdir(java_output_dir) if "CustomRunToMakeAFig" in d]
    if subdirs:
        # Sort by modification time, get most recent
        subdirs.sort(key=lambda x: os.path.getmtime(os.path.join(java_output_dir, x)), reverse=True)
        latest_dir = os.path.join(java_output_dir, subdirs[0])
        csv_path = os.path.join(latest_dir, "pkStateVariables.csv")
        
        if os.path.exists(csv_path):
            print(f"\nFound Java output: {csv_path}")
        else:
            print(f"\nERROR: pkStateVariables.csv not found in {latest_dir}")
            sys.exit(1)
    else:
        print(f"\nERROR: No CustomRunToMakeAFig directories found in {java_output_dir}")
        sys.exit(1)
else:
    print(f"\nERROR: Directory not found: {java_output_dir}")
    sys.exit(1)

# =============================================================================
# LOAD JAVA DATA
# =============================================================================

print("\nLoading Java simulation data...")

# Read CSV with header
df_java = pd.read_csv(csv_path)

# Fixed column order in pkStateVariables.csv
# (must match Java output exactly)
column_mapping = [
    "N_cen_hot",
    "N_cen_cold",
    "N_v_hot",
    "N_v_cold",
    "N_ec_hot",
    "N_ec_cold",
    "N_b_hot",
    "N_b_cold",
    "N_ic_hot",
    "N_ic_cold",
    # "A_blob",  # include if present
]

# Safety check
if df_java.shape[1] < len(column_mapping):
    print("\nERROR: CSV has fewer columns than expected")
    print(f"Expected ≥ {len(column_mapping)}, found {df_java.shape[1]}")
    sys.exit(1)

# Assign column names by position
df_java = df_java.iloc[:, :len(column_mapping)]
df_java.columns = column_mapping

# Create time array (hourly data)
num_hours = len(df_java)
java_time_hours = np.arange(num_hours)
java_time_days = java_time_hours / 24.0

java_N_cen_hot = df_java["N_cen_hot"].values * 1e9  # mol --> nmol conversion
java_N_ic_hot  = df_java["N_ic_hot"].values * 1e9   # mol --> nmol conversion

print(f"  Loaded {len(df_java)} time points (hours)")
print(f"  Time range: {java_time_days[0]:.1f} to {java_time_days[-1]:.1f} days")
print(f"  Total simulation time: {num_hours/24:.1f} days")

# =============================================================================
# RUN REDUCED MODEL
# =============================================================================

print("\nRunning reduced QSS model...")

# Import required modules
from scipy.integrate import odeint

# Parameters (match SimParams.java exactly)
lambda_bio = 1.6e-4 / 60.0
lambda_decay = 7.14e-5 / 60.0
#k_on_mol = 0.046 * 1e6 / 60.0 # units m³/(mol·s) - in the Java simulation, it's m³/(mol·s)
k_on_mol = 0.0015 * 1e6 / 60.0 # units m³/(mol·s) - in the Java simulation, it's m³/(mol·s)
k_on_nmol = k_on_mol * 1e-9 # units m³/(nmol·s) , mol --> nmol conversion
#k_off = 0.368 / 60.0
k_off = 0.012 / 60.0
k_int = 0.001 / 60.0
k_rel = 2e-4 / 60.0
beta_nmol = (k_off + k_int) / k_on_nmol

# Geometry
CELL_LENGTH = 1e-5
INTERSTITIAL_FRACTION = 0.4
AVOGADRO = 6.022e23
RECEPTORS_PER_CELL = 3e5

radius_cells = np.sqrt(TUMOR_CELLS_2D / np.pi)
radius_m = radius_cells * CELL_LENGTH
height_m = 2.0 * radius_m
tumourVolume_m3 = np.pi * radius_m**2 * height_m
V_ec = INTERSTITIAL_FRACTION * tumourVolume_m3
V_cen = 0.5e-3
height_cells = height_m / CELL_LENGTH
totalCells3D = TUMOR_CELLS_2D * height_cells
R_total_mol = totalCells3D * RECEPTORS_PER_CELL / AVOGADRO
R_total_nmol = R_total_mol * 1e9 # mol --> nmol conversion

# Initial conditions
dose_mol = 100e-9
dose_nmol = dose_mol* 1e9  # mol --> mnmol conversion
hot_fraction = 0.1
N_cen0_nmol = dose_nmol
N_cen0_H_nmol = hot_fraction * dose_nmol
beta_scaled_nmol = V_cen * beta_nmol

# Time array
t_max_days = 30
t_max_sec = t_max_days * 24 * 3600
t_sec = np.arange(0, t_max_sec, 3600)
t_days = t_sec / 86400

from scipy.integrate import solve_ivp

# ODE system
#def odes_reduced(t, y):
def odes_reduced(y, t):
    N_cen_H, N_cen_C, N_ic_H, N_ic_C = y
    N_cen = N_cen_H + N_cen_C
    N_b_H = R_total_nmol * N_cen_H / (N_cen + beta_scaled_nmol)
    N_b_C = R_total_nmol * N_cen_C / (N_cen + beta_scaled_nmol)
    
    dN_cen_H = -lambda_bio * N_cen_H - k_int * N_b_H - lambda_decay * N_cen_H
    dN_cen_C = -lambda_bio * N_cen_C - k_int * N_b_C + lambda_decay * N_cen_H
    dN_ic_H = k_int * N_b_H - k_rel * N_ic_H - lambda_decay * N_ic_H
    dN_ic_C = k_int * N_b_C - k_rel * N_ic_C + lambda_decay * N_ic_H
    
    return [dN_cen_H, dN_cen_C, dN_ic_H, dN_ic_C]

# Integrate
y0_reduced_nmol = [N_cen0_H_nmol, N_cen0_nmol - N_cen0_H_nmol, 0.0, 0.0]
#sol_reduced = solve_ivp(odes_reduced, (0, t_max_sec), y0_reduced, t_eval=t_sec, method='Radau', rtol=1e-10, atol=1e-12)
#reduced_N_cen_hot = sol_reduced.y[0, :]
#reduced_N_ic_hot  = sol_reduced.y[2, :]
sol_reduced = odeint(odes_reduced, y0_reduced_nmol, t_sec)
reduced_N_cen_hot = sol_reduced[:, 0]
reduced_N_ic_hot = sol_reduced[:, 2]

print(f"  Integrated {len(t_sec)} time points")

# =============================================================================
# RUN FULL ODE MODEL
# =============================================================================

print("\nRunning full ODE model...")

from scipy.integrate import solve_ivp

# Flow/permeability parameters
numberOfVessels = 365
VESSEL_FLOW = 5e-9 / 1000.0
VESSEL_PS = (6.0/5.0) * VESSEL_FLOW
vesselVolume = (CELL_LENGTH**2) * height_m
V_v = numberOfVessels * vesselVolume
F_T = numberOfVessels * VESSEL_FLOW
PS = numberOfVessels * VESSEL_PS

# Full ODE system
def odes_full(t, y):
    N_cen_H, N_cen_C = y[0], y[1]
    N_v_H, N_v_C = y[2], y[3]
    N_ec_H, N_ec_C = y[4], y[5]
    N_b_H, N_b_C = y[6], y[7]
    N_ic_H, N_ic_C = y[8], y[9]
    
    C_cen_H = N_cen_H / V_cen
    C_cen_C = N_cen_C / V_cen
    C_v_H = N_v_H / V_v
    C_v_C = N_v_C / V_v
    C_ec_H = N_ec_H / V_ec
    C_ec_C = N_ec_C / V_ec
    
    N_R_free = R_total_nmol - (N_b_H + N_b_C)
    C_R_free = max(0, N_R_free / V_ec)
    
    dy = np.zeros(10)
    
    # Hot compartments
    dy[0] = -lambda_bio * N_cen_H - F_T * (C_cen_H - C_v_H) - lambda_decay * N_cen_H
    dy[2] = F_T * (C_cen_H - C_v_H) - PS * (C_v_H - C_ec_H) - lambda_decay * N_v_H
    dy[4] = PS * (C_v_H - C_ec_H) - k_on_nmol * C_ec_H * C_R_free * V_ec + k_off * N_b_H - lambda_decay * N_ec_H
    dy[6] = k_on_nmol * C_ec_H * C_R_free * V_ec - k_off * N_b_H - k_int * N_b_H - lambda_decay * N_b_H
    dy[8] = k_int * N_b_H - k_rel * N_ic_H - lambda_decay * N_ic_H
    
    # Cold compartments
    dy[1] = -lambda_bio * N_cen_C - F_T * (C_cen_C - C_v_C) + lambda_decay * N_cen_H
    dy[3] = F_T * (C_cen_C - C_v_C) - PS * (C_v_C - C_ec_C) + lambda_decay * N_v_H
    dy[5] = PS * (C_v_C - C_ec_C) - k_on_nmol * C_ec_C * C_R_free * V_ec + k_off * N_b_C + lambda_decay * N_ec_H
    dy[7] = k_on_nmol * C_ec_C * C_R_free * V_ec - k_off * N_b_C - k_int * N_b_C + lambda_decay * N_b_H
    dy[9] = k_int * N_b_C - k_rel * N_ic_C + lambda_decay * N_ic_H
    
    return dy

# Initial conditions and integration
y0_full = np.zeros(10)
y0_full[0] = hot_fraction * dose_nmol
y0_full[1] = (1 - hot_fraction) * dose_nmol

sol_full = solve_ivp(
    odes_full,
    (0, t_max_sec),
    y0_full,
    t_eval=t_sec,
    method='BDF',
    max_step=3600.0,
    rtol=1e-7,
    atol = 1e-8 #* np.maximum(np.abs(y0_full), 1.0)
    )

full_N_cen_hot = sol_full.y[0, :]
full_N_ic_hot = sol_full.y[8, :]

print(f"  Integrated {len(sol_full.t)} time points")

# =============================================================================
# CREATE COMPARISON PLOT
# =============================================================================

print("\nCreating comparison plot...")

fig, axes = plt.subplots(2, 1)

# N_cen_hot comparison
ax1 = axes[0]
ax1.plot(java_time_days, java_N_cen_hot, label='Java QSS', markersize=1, linewidth=2, alpha=0.7)
ax1.plot(t_days, reduced_N_cen_hot, '--', label='Reduced Model (QSS)', linewidth=1)
ax1.plot(t_days, full_N_cen_hot, ':', label='Full ODE', linewidth=1)

#ax1.set_xlabel('Time (days)', fontsize=12)
ax1.set_ylabel(r"$N_{\mathrm{cen}}^{\mathrm{hot}}$ (nmol)")
#ax1.set_ylabel('N_cen_hot (mol)', fontsize=12)
#ax1.set_title('Central Compartment (Hot)', fontsize=13, fontweight='bold')
ax1.set_title('PK Model Validation', fontweight='bold')

#ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# N_ic_hot comparison
ax2 = axes[1]
ax2.plot(java_time_days, java_N_ic_hot, label='Java QSS', markersize=1, linewidth=2, alpha=0.7)
ax2.plot(t_days, reduced_N_ic_hot, '--', label='Reduced Model (QSS)', linewidth=1)
ax2.plot(t_days, full_N_ic_hot, ':', label='Full ODE', linewidth=1)

ax2.set_xlabel('Time (days)')
ax2.set_ylabel(r"$N_{\mathrm{ic}}^{\mathrm{hot}}$ (nmol)")
#ax2.set_ylabel('N_ic_hot (mol)', fontsize=12)
#ax2.set_title('Intracellular Compartment (Hot)', fontsize=13, fontweight='bold')
#ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

#plt.suptitle(f'PK Model Validation\n{TUMOR_CELLS_2D} cells (2D), frozen tumor, 100 nmol injection',
#             fontsize=14, fontweight='bold')
plt.tight_layout()

# Save figure
output_dir = 'results/compare_models'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'pk_model_comparison.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_file}")

# =============================================================================
# PRINT COMPARISON STATISTICS
# =============================================================================

print("\n" + "="*70)
print("COMPARISON STATISTICS")
print("="*70)

# Find matching time points (hourly in Python, hourly in Java)
# Interpolate Java to Python time grid for comparison
from scipy.interpolate import interp1d
java_interp_cen = interp1d(java_time_days, java_N_cen_hot, kind='linear', fill_value='extrapolate')
java_interp_ic = interp1d(java_time_days, java_N_ic_hot, kind='linear', fill_value='extrapolate')

java_cen_matched = java_interp_cen(t_days)
java_ic_matched = java_interp_ic(t_days)

# Compute relative errors
rel_err_reduced_cen = np.abs(reduced_N_cen_hot - java_cen_matched) / java_cen_matched
rel_err_reduced_ic = np.abs(reduced_N_ic_hot - java_ic_matched) / java_ic_matched
rel_err_full_cen = np.abs(full_N_cen_hot - java_cen_matched) / java_cen_matched
rel_err_full_ic = np.abs(full_N_ic_hot - java_ic_matched) / java_ic_matched

print("\nN_cen_hot relative error (compared to Java QSS):")
print(f"  Reduced model: mean = {np.mean(rel_err_reduced_cen)*100:.2f}%, max = {np.max(rel_err_reduced_cen)*100:.2f}%")
print(f"  Full ODE:      mean = {np.mean(rel_err_full_cen)*100:.2f}%, max = {np.max(rel_err_full_cen)*100:.2f}%")

print("\nN_ic_hot relative error (compared to Java QSS):")
print(f"  Reduced model: mean = {np.mean(rel_err_reduced_ic)*100:.2f}%, max = {np.max(rel_err_reduced_ic)*100:.2f}%")
print(f"  Full ODE:      mean = {np.mean(rel_err_full_ic)*100:.2f}%, max = {np.max(rel_err_full_ic)*100:.2f}%")

# Sample values at key time points
print("\nSample values at key time points:")
for day in [0, 1, 7, 14, 30]:
    idx = np.argmin(np.abs(t_days - day))
    if idx < len(t_days):
        print(f"\nDay {day}:")
        print(f"  N_cen_hot: Java={java_cen_matched[idx]:.4e}, Reduced={reduced_N_cen_hot[idx]:.4e}, Full={full_N_cen_hot[idx]:.4e}")
        print(f"  N_ic_hot:  Java={java_ic_matched[idx]:.4e}, Reduced={reduced_N_ic_hot[idx]:.4e}, Full={full_N_ic_hot[idx]:.4e}")

print("\n" + "="*70)
print("COMPARISON COMPLETE")
print("="*70 + "\n")

plt.show()
