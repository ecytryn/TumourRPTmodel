#!/usr/bin/env python3
"""
Compare Full Java PK Model vs Reduced Analytical Model

Loads:
- Java simulation output (state_variable.csv)
- Reduced analytical solution (computed here)

Plots them overlaid with error metrics.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from pathlib import Path

# =============================================================================
# PARAMETERS - MUST MATCH BOTH MODELS
# =============================================================================

# Rate constants (SI units: 1/s, m³/(mol·s))
lambda_bio = 1.6e-4 / 60.0
lambda_decay = 7.14e-5 / 60.0
k_on = 0.046 * 1e6 / 60.0
k_off = 0.368 / 60.0
k_int = 0.001 / 60.0
k_rel = 2e-4 / 60.0

beta = (k_off + k_int) / k_on

# Geometry (from Java output)
tumour_cells_2D = 3285
numberOfVessels = 365

CELL_LENGTH = 1e-5
INTERSTITIAL_FRACTION = 0.4

# Calculate volumes
radius_cells = np.sqrt(tumour_cells_2D / np.pi)
radius_m = radius_cells * CELL_LENGTH
height_m = 2.0 * radius_m
tumourVolume_m3 = np.pi * radius_m**2 * height_m
V_ec_m3 = INTERSTITIAL_FRACTION * tumourVolume_m3
vesselVolumeInTumor_m3 = (CELL_LENGTH**2) * height_m
V_vasc_m3 = numberOfVessels * vesselVolumeInTumor_m3
V_cen_m3 = 0.458e-3

# Receptors
AVOGADRO = 6.022e23
RECEPTORS_PER_CELL = 5e5
height_cells = height_m / CELL_LENGTH
totalCells3D = tumour_cells_2D * height_cells
R_total_mol = totalCells3D * RECEPTORS_PER_CELL / AVOGADRO
R_T_tilde_mol_m3 = R_total_mol / V_ec_m3

# Initial conditions
dose_mol = 100e-9
hot_fraction = 0.1
C_cen0 = dose_mol / V_cen_m3
C_cen0_H = hot_fraction * C_cen0
gamma = beta / C_cen0

# =============================================================================
# LOAD JAVA DATA
# =============================================================================

java_file = "results/single_runs/pkStateVariables.csv"

if not Path(java_file).exists():
    print(f"ERROR: Java output not found: {java_file}")
    print("Run the Java simulation first!")
    exit(1)

java_data = np.loadtxt(java_file, delimiter=',')
print(f"✓ Loaded Java data: {java_file}")
print(f"  Shape: {java_data.shape}")
print(f"  Duration: {len(java_data)} hours = {len(java_data)/24:.1f} days")

# Extract Java compartments (concentrations in mol/m³)
java_t_hours = np.arange(len(java_data))
java_t_days = java_t_hours / 24.0

java_C_cen_H = java_data[:, 0]    # Central hot
java_C_v_H = java_data[:, 2]      # Vascular hot
java_C_ec_H = java_data[:, 4]     # Extracellular hot
java_C_b_H = java_data[:, 6]      # Bound hot
java_C_ic_H = java_data[:, 8]     # Intracellular hot

# =============================================================================
# COMPUTE ANALYTICAL SOLUTION
# =============================================================================

print("\n✓ Computing analytical solution...")

t_max_sec = len(java_data) * 3600  # Match Java duration
dt = 3600
t_sec = np.arange(0, t_max_sec, dt)
t_days = t_sec / (24 * 3600)

# Central
C_cen = C_cen0 * np.exp(-lambda_bio * t_sec)
C_cen_H = C_cen0_H * np.exp(-(lambda_bio + lambda_decay) * t_sec)

# Bound (QSS)
C_b = R_T_tilde_mol_m3 * C_cen / (C_cen + beta)
C_b_H = R_T_tilde_mol_m3 * C_cen_H / (C_cen + beta)

# Intracellular (numerical integration)
def integrand(s):
    return np.exp(k_rel * s) * np.exp(-lambda_bio * s) / (np.exp(-lambda_bio * s) + gamma)

integral_values = cumulative_trapezoid(integrand(t_sec), t_sec, initial=0)
C_ic = k_int * R_T_tilde_mol_m3 * np.exp(-k_rel * t_sec) * integral_values
C_ic_H = hot_fraction*np.exp(-lambda_decay * t_sec) * C_ic

# =============================================================================
# NUMERICAL ODE SOLUTION (Forward Euler for validation)
# =============================================================================

print("\n✓ Computing numerical ODE solution (Forward Euler)...")

# Time parameters (same as analytical)
dt_euler = 3600  # 1 hour time step (seconds)
n_steps = int(t_max_sec / dt_euler)

# Initialize arrays
t_euler = np.zeros(n_steps)
Euler_C_cen = np.zeros(n_steps)
Euler_C_cen_H = np.zeros(n_steps)
Euler_C_b = np.zeros(n_steps)
Euler_C_b_H = np.zeros(n_steps)
Euler_C_ic = np.zeros(n_steps)
Euler_C_ic_H = np.zeros(n_steps)

# Initial conditions (at t=0, immediately after injection)
Euler_C_cen[0] = C_cen0
Euler_C_cen_H[0] = C_cen0_H
Euler_C_ic[0] = 0.0
Euler_C_ic_H[0] = 0.0

# Forward Euler integration
for i in range(n_steps - 1):
    t_euler[i+1] = t_euler[i] + dt_euler
    
    # Current state
    Euler_C_cen_tot = Euler_C_cen[i]
    Euler_C_cen_hot = Euler_C_cen_H[i]
    Euler_C_cen_cold = Euler_C_cen_tot - Euler_C_cen_hot
    Euler_C_ic_tot = Euler_C_ic[i]
    Euler_C_ic_hot = Euler_C_ic_H[i]
    
    # Compute QSS values (bound concentrations)
    Euler_C_b_tot = R_T_tilde_mol_m3 * Euler_C_cen_tot / (Euler_C_cen_tot + beta)
    Euler_C_b_hot = R_T_tilde_mol_m3 * Euler_C_cen_hot / (Euler_C_cen_tot + beta)
    
    # Store for plotting
    Euler_C_b[i] = Euler_C_b_tot
    Euler_C_b_H[i] = Euler_C_b_hot
    
    # ODEs for total
    dC_cen_tot_dt = -lambda_bio * Euler_C_cen_tot
    dC_ic_tot_dt = k_int * Euler_C_b_tot - k_rel * Euler_C_ic_tot
    
    # ODEs for hot
    dC_cen_hot_dt = -lambda_bio * Euler_C_cen_hot - lambda_decay * Euler_C_cen_hot
    dC_ic_hot_dt = k_int * Euler_C_b_hot - k_rel * Euler_C_ic_hot - lambda_decay * Euler_C_ic_hot
    
    # Forward Euler step
    Euler_C_cen[i+1] = Euler_C_cen_tot + dC_cen_tot_dt * dt_euler
    Euler_C_cen_H[i+1] = Euler_C_cen_hot + dC_cen_hot_dt * dt_euler
    Euler_C_ic[i+1] = Euler_C_ic_tot + dC_ic_tot_dt * dt_euler
    Euler_C_ic_H[i+1] = Euler_C_ic_hot + dC_ic_hot_dt * dt_euler

# Fill last QSS values
Euler_C_b[-1] = R_T_tilde_mol_m3 * Euler_C_cen[-1] / (Euler_C_cen[-1] + beta)
Euler_C_b_H[-1] = R_T_tilde_mol_m3 * Euler_C_cen_H[-1] / (Euler_C_cen[-1] + beta)

print(f"  Euler integration complete: {n_steps} steps")




# =============================================================================
# DETAILED ODE TERM COMPARISON (Add after the Euler integration section)
# =============================================================================

print("\n" + "="*70)
print("ODE TERM-BY-TERM COMPARISON")
print("="*70)

# Compute individual terms for JAVA
java_term_kinton = k_int * java_C_b_H
java_term_krel = k_rel * java_C_ic_H
java_term_decay = lambda_decay * java_C_ic_H
java_dCic_dt = java_term_kinton - java_term_krel - java_term_decay

# Compute individual terms for EULER
euler_term_kinton = k_int * Euler_C_b_H
euler_term_krel = k_rel * Euler_C_ic_H
euler_term_decay = lambda_decay * Euler_C_ic_H
euler_dCic_dt = euler_term_kinton - euler_term_krel - euler_term_decay

# Print comparison at key time points
print("\nComparison at key hours:")
test_hours = [0, 1, 2, 3, 5, 10, 24, 48, 120]  # Hours after injection (day 5)

for h in test_hours:
    if h >= len(java_data):
        continue
        
    print(f"  k_int * C_b_H:")
    print(f"    Java:  {java_term_kinton[h]:.6e} mol/(m³·s)")
    print(f"    Euler: {euler_term_kinton[h]:.6e} mol/(m³·s)")
    ratio = java_term_kinton[h]/euler_term_kinton[h] if euler_term_kinton[h] > 1e-20 else None
    print(f"    Ratio: {ratio:.4f}" if ratio is not None else "    Ratio: N/A")
    
    print(f"  k_rel * C_ic_H:")
    print(f"    Java:  {java_term_krel[h]:.6e} mol/(m³·s)")
    print(f"    Euler: {euler_term_krel[h]:.6e} mol/(m³·s)")
    ratio = java_term_krel[h]/euler_term_krel[h] if euler_term_krel[h] > 1e-20 else None
    print(f"    Ratio: {ratio:.4f}" if ratio is not None else "    Ratio: N/A")
    
    print(f"  lambda_decay * C_ic_H:")
    print(f"    Java:  {java_term_decay[h]:.6e} mol/(m³·s)")
    print(f"    Euler: {euler_term_decay[h]:.6e} mol/(m³·s)")
    ratio = java_term_decay[h]/euler_term_decay[h] if euler_term_decay[h] > 1e-20 else None
    print(f"    Ratio: {ratio:.4f}" if ratio is not None else "    Ratio: N/A")
    
    print(f"  Net dC_ic_H/dt:")
    print(f"    Java:  {java_dCic_dt[h]:.6e} mol/(m³·s)")
    print(f"    Euler: {euler_dCic_dt[h]:.6e} mol/(m³·s)")
    ratio = java_dCic_dt[h]/euler_dCic_dt[h] if abs(euler_dCic_dt[h]) > 1e-20 else None
    print(f"    Ratio: {ratio:.4f}" if ratio is not None else "    Ratio: N/A")






    
# =============================================================================
# CONVERT TO AMOUNTS FOR PLOTTING
# =============================================================================

# Java amounts (mol → nmol)
java_A_cen_H = java_C_cen_H * V_cen_m3 * 1e9
java_A_v_H = java_C_v_H * V_vasc_m3 * 1e9
java_A_ec_H = java_C_ec_H * V_ec_m3 * 1e9
java_A_b_H = java_C_b_H * V_ec_m3 * 1e9
java_A_ic_H = java_C_ic_H * V_ec_m3 * 1e9
java_A_tumor_H = java_A_v_H + java_A_ec_H + java_A_b_H + java_A_ic_H
java_A_captive_H = java_A_b_H + java_A_ic_H
java_A_free_H = java_A_v_H + java_A_ec_H

# Analytical amounts (mol → nmol)
ana_A_cen = C_cen * V_cen_m3 * 1e9
ana_A_cen_H = C_cen_H * V_cen_m3 * 1e9
ana_A_v_H = C_cen_H * V_vasc_m3 * 1e9  # C_v = C_cen in QSS
ana_A_ec_H = C_cen_H * V_ec_m3 * 1e9   # C_ec = C_cen in QSS
ana_A_b = C_b * V_ec_m3 * 1e9
ana_A_b_H = C_b_H * V_ec_m3 * 1e9
ana_A_ic = C_ic * V_ec_m3 * 1e9
ana_A_ic_H = C_ic_H * V_ec_m3 * 1e9
ana_A_tumor_H = ana_A_v_H + ana_A_ec_H + ana_A_b_H + ana_A_ic_H
ana_A_captive_H = ana_A_b_H + ana_A_ic_H
ana_A_free_H = ana_A_v_H + ana_A_ec_H


# Convert numerical solutions to amounts for plotting
Euler_A_cen = Euler_C_cen * V_cen_m3 * 1e9
Euler_A_b = Euler_C_b * V_ec_m3 * 1e9
Euler_A_ic = Euler_C_ic * V_ec_m3 * 1e9
Euler_A_cen_H = Euler_C_cen_H * V_cen_m3 * 1e9
Euler_A_b_H = Euler_C_b_H * V_ec_m3 * 1e9
Euler_A_ic_H = Euler_C_ic_H * V_ec_m3 * 1e9

# =============================================================================
# PLOTTING
# =============================================================================

print("\n✓ Generating comparison plots...")

t_days_euler = t_euler / (24 * 3600)

fig, axes = plt.subplots(3, 1, figsize=(14, 12))

# Row 1: Central compartment
#ax = axes[0, 0]
ax = axes[0]
ax.plot(java_t_days, java_A_cen_H, 'r-', linewidth=2, label='Java (full)', alpha=0.7)
ax.plot(t_days, ana_A_cen_H, 'b--', linewidth=1.5, label='Analytical (reduced)')
ax.plot(t_days, ana_A_cen, 'b:', linewidth=1.5, label='Analytical (reduced, tot)')
ax.plot(t_days_euler, Euler_A_cen_H, 'g--', linewidth=2, label='Numerical (Euler)', alpha=0.7)
ax.plot(t_days_euler, Euler_A_cen, 'g:', linewidth=2, label='Numerical (Euler)', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Amount (nmol)')
ax.set_title('Central Compartment (Hot)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

#ax = axes[0, 1]
#ax.plot(java_t_days, err_cen * 100, 'k-', linewidth=1.5)
#ax.set_xlabel('Time (days)')
#ax.set_ylabel('Relative Error (%)')
#ax.set_title('Central: Error')
#ax.axhline(1, color='r', linestyle='--', alpha=0.5, label='1%')
#ax.axhline(5, color='orange', linestyle='--', alpha=0.5, label='5%')
#ax.legend()
#ax.grid(True, alpha=0.3)

# Row 2: Bound
#ax = axes[1, 0]
ax = axes[1]
ax.plot(java_t_days, java_A_b_H, 'r.', linewidth=2, label='Java (full)', alpha=0.7)
ax.plot(t_days, ana_A_b_H, 'b--', linewidth=1.5, label='Analytical (reduced)')
ax.plot(t_days, ana_A_b, 'b:', linewidth=1.5, label='Analytical (reduced, tot)')
ax.plot(t_days_euler, Euler_A_b_H, 'g--', linewidth=2, label='Numerical (Euler)', alpha=0.7)
ax.plot(t_days_euler, Euler_A_b, 'g:', linewidth=2, label='Numerical (Euler)', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Amount (nmol)')
ax.set_title('Bound (Hot)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

#ax = axes[1, 1]
#ax.plot(java_t_days, err_b * 100, 'k-', linewidth=1.5)
#ax.set_xlabel('Time (days)')
#ax.set_ylabel('Relative Error (%)')
#ax.set_title('Bound: Error')
#ax.axhline(1, color='r', linestyle='--', alpha=0.5, label='1%')
#ax.axhline(5, color='orange', linestyle='--', alpha=0.5, label='5%')
#ax.legend()
#ax.grid(True, alpha=0.3)

# Row 3: Intracellular
#ax = axes[2, 0]
ax = axes[2]
ax.plot(java_t_days, java_A_ic_H, 'r.', linewidth=2, label='Java (full)', alpha=0.7)
ax.plot(t_days, ana_A_ic_H, 'b--', linewidth=1.5, label='Analytical (reduced, Hot)')
#ax.plot(t_days, ana_A_ic, 'b:', linewidth=1.5, label='Analytical (reduced, Tot)')
ax.plot(t_days_euler, Euler_A_ic_H, 'g.', linewidth=2, label='Numerical (Euler, Hot)', alpha=0.7)
#ax.plot(t_days_euler, Euler_A_ic, 'g:', linewidth=2, label='Numerical (Euler, Tot)', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Amount (nmol)')
ax.set_title('Intracellular (Hot)')
#ax.legend()
ax.grid(True, alpha=0.3)
#ax.set_yscale('log')

#ax = axes[2, 1]
#ax.plot(java_t_days, err_ic * 100, 'k-', linewidth=1.5)
#ax.set_xlabel('Time (days)')
#ax.set_ylabel('Relative Error (%)')
#ax.set_title('Intracellular: Error')
#ax.axhline(1, color='r', linestyle='--', alpha=0.5, label='1%')
#ax.axhline(5, color='orange', linestyle='--', alpha=0.5, label='5%')
#ax.legend()
#ax.grid(True, alpha=0.3)

plt.suptitle(f'PK Comparison: Java Full Model vs Analytical Reduced Model\n' +
             f'{tumour_cells_2D} cells, frozen tumor, single injection',
             fontsize=13, fontweight='bold')

plt.tight_layout()

output_file = 'results/compare_models/pk_model_comparison.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Figure saved: {output_file}")







# =============================================================================
# PLOT ODE TERMS (Add as a new figure before plt.show())
# =============================================================================

print("\n✓ Generating ODE term comparison plots...")

fig_terms, axes_terms = plt.subplots(4, 1, figsize=(14, 16))

# Plot 1: Internalization term (k_int * C_b_H)
ax = axes_terms[0]
ax.plot(java_t_days, java_term_kinton * 1e9, 'r-', linewidth=2, label='Java', alpha=0.7)
ax.plot(t_days_euler, euler_term_kinton * 1e9, 'g--', linewidth=2, label='Euler', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Rate (nmol/(L·s))')
ax.set_title('Internalization Term: k_int × C_b_H')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# Plot 2: Release term (k_rel * C_ic_H)
ax = axes_terms[1]
ax.plot(java_t_days, java_term_krel * 1e9, 'r-', linewidth=2, label='Java', alpha=0.7)
ax.plot(t_days_euler, euler_term_krel * 1e9, 'g--', linewidth=2, label='Euler', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Rate (nmol/(L·s))')
ax.set_title('Release Term: k_rel × C_ic_H')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# Plot 3: Decay term (lambda_decay * C_ic_H)
ax = axes_terms[2]
ax.plot(java_t_days, java_term_decay * 1e9, 'r-', linewidth=2, label='Java', alpha=0.7)
ax.plot(t_days_euler, euler_term_decay * 1e9, 'g--', linewidth=2, label='Euler', alpha=0.7)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Rate (nmol/(L·s))')
ax.set_title('Decay Term: λ_decay × C_ic_H')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# Plot 4: Net derivative
ax = axes_terms[3]
ax.plot(java_t_days, java_dCic_dt * 1e9, 'r-', linewidth=2, label='Java', alpha=0.7)
ax.plot(t_days_euler, euler_dCic_dt * 1e9, 'g--', linewidth=2, label='Euler', alpha=0.7)
ax.axhline(0, color='k', linestyle=':', alpha=0.5)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Rate (nmol/(L·s))')
ax.set_title('Net Derivative: dC_ic_H/dt = k_int×C_b_H - k_rel×C_ic_H - λ×C_ic_H')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle(f'ODE Term Comparison: Java vs Euler\n' +
             f'Intracellular Hot Compartment Dynamics',
             fontsize=13, fontweight='bold')

plt.tight_layout()

output_file_terms = 'results/compare_models/pk_ode_terms_comparison.png'
plt.savefig(output_file_terms, dpi=300, bbox_inches='tight')
print(f"✓ ODE terms figure saved: {output_file_terms}")

# =============================================================================
# RATIO PLOTS (to see divergence clearly)
# =============================================================================

fig_ratio, axes_ratio = plt.subplots(2, 1, figsize=(14, 10))

# Plot 1: Ratio of C_b_H (should be ~1.0)
ax = axes_ratio[0]
ratio_Cb = java_C_b_H / (Euler_C_b_H + 1e-30)  # Avoid division by zero
ax.plot(java_t_days, ratio_Cb, 'b-', linewidth=2)
ax.axhline(1.0, color='r', linestyle='--', label='Perfect match')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Ratio (Java / Euler)')
ax.set_title('Bound Concentration Ratio: C_b_H(Java) / C_b_H(Euler)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim([0.9, 1.1])  # Zoom in to see small deviations

# Plot 2: Ratio of C_ic_H (should be ~1.0 but isn't!)
ax = axes_ratio[1]
ratio_Cic = java_C_ic_H / (Euler_C_ic_H + 1e-30)
ax.plot(java_t_days, ratio_Cic, 'b-', linewidth=2)
ax.axhline(1.0, color='r', linestyle='--', label='Perfect match')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Ratio (Java / Euler)')
ax.set_title('Intracellular Concentration Ratio: C_ic_H(Java) / C_ic_H(Euler)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

plt.suptitle('Concentration Ratios: Where Does Java Diverge?',
             fontsize=13, fontweight='bold')

plt.tight_layout()

output_file_ratio = 'results/compare_models/pk_ratio_comparison.png'
plt.savefig(output_file_ratio, dpi=300, bbox_inches='tight')
print(f"✓ Ratio figure saved: {output_file_ratio}")








plt.show()