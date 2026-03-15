#!/usr/bin/env python3
"""
Study the effect of binding affinity on cumulative activity

Vary k_off (unbinding rate) to explore affinity landscape.
Lower k_off = higher affinity (stronger binding)

Cumulative activity metric: Integral of captive radioligand (N_b^H + N_ic^H) over 30 days
This represents cumulative radiation exposure to the tumor.

Note: k_int sets a floor on benefit - once internalization becomes limiting,
further affinity increases provide diminishing returns.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import trapezoid
from matplotlib.ticker import ScalarFormatter
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
# FIXED PARAMETERS (from SimParams.java)
# =============================================================================

# Rate constants (converted to /min and nmol for numerical stability)
lambda_bio = 1.6e-4  # /min (biological clearance)
lambda_decay = 7.14e-5  # /min (Lu-177 decay)
#k_on = 0.046 * 1e6 * 1e-9 # m³/(nmol·min) (binding rate)
k_on = 0.0015 * 1e6 * 1e-9 # m³/(nmol·min) (binding rate)
# k_off will be varied
k_int = 0.001  # /min (internalization)
k_rel = 2e-4  # /min (release from cells)

# Geometry (for small tumor, ~3285 cells in 2D)
CELL_LENGTH = 1e-5  # m
INTERSTITIAL_FRACTION = 0.4
AVOGADRO = 6.022e23
RECEPTORS_PER_CELL = 3e5
TUMOR_CELLS_2D = 3285

radius_cells = np.sqrt(TUMOR_CELLS_2D / np.pi)
radius_m = radius_cells * CELL_LENGTH
height_m = 2.0 * radius_m
tumourVolume_m3 = np.pi * radius_m**2 * height_m
V_ec = INTERSTITIAL_FRACTION * tumourVolume_m3  # m³
V_cen = 0.5e-3  # m³

height_cells = height_m / CELL_LENGTH
totalCells3D = TUMOR_CELLS_2D * height_cells
R_total = totalCells3D * RECEPTORS_PER_CELL / AVOGADRO * 1e9  # nmol

# Initial conditions
dose_nmol = 100.0  # nmol total injection
#dose_nmol = 0.01  # nmol total injection
hot_fraction = 0.1
N_cen0 = dose_nmol
N_cen0_H = hot_fraction * dose_nmol

# Time array
t_max_days = 30
t_max_min = t_max_days * 24 * 60
t_min = np.arange(0, t_max_min, 60)  # Every hour
t_days = t_min / (24 * 60)

# =============================================================================
# DEFINE REDUCED PK MODEL
# =============================================================================

def solve_pk_model(k_off_value):
    """
    Solve reduced PK model for a given k_off value
    
    Returns:
        N_b_hot: Bound radioligand (nmol) over time
        N_ic_hot: Intracellular radioligand (nmol) over time
        captive_integral: Integral of (N_b + N_ic) over time (nmol·day)
    """
    
    # Compute beta for this k_off
    beta = (k_off_value + k_int) / k_on  # nmol/m³
    beta_scaled = V_cen * beta  # nmol
    
    def odes_reduced(y, t):
        """Reduced QSS PK model"""
        N_cen_H, N_cen_C, N_ic_H, N_ic_C = y
        N_cen = N_cen_H + N_cen_C
        
        # QSS for bound
        N_b_H = R_total * N_cen_H / (N_cen + beta_scaled)
        N_b_C = R_total * N_cen_C / (N_cen + beta_scaled)
        
        # ODEs
        dN_cen_H = -lambda_bio * N_cen_H - k_int * N_b_H - lambda_decay * N_cen_H
        dN_cen_C = -lambda_bio * N_cen_C - k_int * N_b_C + lambda_decay * N_cen_H
        dN_ic_H = k_int * N_b_H - k_rel * N_ic_H - lambda_decay * N_ic_H
        dN_ic_C = k_int * N_b_C - k_rel * N_ic_C + lambda_decay * N_ic_H
        
        return [dN_cen_H, dN_cen_C, dN_ic_H, dN_ic_C]
    
    # Integrate
    y0 = [N_cen0_H, N_cen0 - N_cen0_H, 0.0, 0.0]
    sol = odeint(odes_reduced, y0, t_min)
    
    N_cen_H = sol[:, 0]
    N_cen_C = sol[:, 1]
    N_ic_H = sol[:, 2]
    N_ic_C = sol[:, 3]
    N_cen = N_cen_H + N_cen_C
    
    # Compute N_b via QSS
    N_b_H = R_total * N_cen_H / (N_cen + beta_scaled)
    
    # Captive radioligand (bound + internalized)
    N_captive_H = N_b_H + N_ic_H
    
    # Integrate over time (convert time to days for the integral)
    CumAct = lambda_decay * trapezoid(N_captive_H, t_min) * 1e-9 * AVOGADRO / 1e9  # GBq·s
    return N_b_H, N_ic_H, CumAct

# =============================================================================
# SWEEP OVER k_off VALUES
# =============================================================================

print("\n" + "="*70)
print("AFFINITY STUDY: k_off vs Cumulative activity")
print("="*70)

# k_off range: Start below k_int (0.001 /min), go up to ~0.4 /min
# Use log spacing (doubling)
k_off_min = 0.0001  # /min (very high affinity)
k_off_max = 100*0.4     # /min (low affinity)

# Generate log-spaced values (doubling)
num_doublings = int(np.log2(k_off_max / k_off_min)) + 1
k_off_values = k_off_min * (2.0 ** np.arange(num_doublings))

# Also add k_int for reference
k_off_values = np.sort(np.append(k_off_values, k_int))
#k_off_values = np.sort(np.append(k_off_values, 0.386))
k_off_values = np.sort(np.append(k_off_values, 0.012))

print(f"\nScanning {len(k_off_values)} k_off values from {k_off_min} to {k_off_max:.2f} /min")
print(f"k_int = {k_int} /min (internalization rate - sets limiting timescale)")
print(f"k_on = {k_on:.2e} m³/(nmol·min)")
print(f"R_total = {R_total:.4e} nmol")

# Storage for results
CumAct_values = []
K_d_values = []  # Dissociation constant (k_off / k_on)

print("\nRunning simulations...")
for i, k_off in enumerate(k_off_values):
    N_b_H, N_ic_H, CumAct = solve_pk_model(k_off)
    CumAct_values.append(CumAct)
    K_d_values.append(k_off / k_on)  # nmol/m³
    
#    if i % 3 == 0 or k_off == k_int:  # Print every few values
    print(f"  k_off = {k_off:.4f} /min, K_d = {K_d_values[-1]:.2e} nmol/m³, "
      f"Cumulative Activity = {CumAct:.4f} GBq·s")

CumAct_values = np.array(CumAct_values)
K_d_values = np.array(K_d_values)

# Normalize cumulative activity to [0, 1] for easier interpretation
CumAct_norm = CumAct_values / CumAct_values.max()

print("\nSimulations complete!")

# =============================================================================
# ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS")
print("="*70)

# Find k_off that gives 95% of maximum cumulative activity
idx_95 = np.argmin(np.abs(CumAct_norm - 0.95))
k_off_95 = k_off_values[idx_95]
K_d_95 = K_d_values[idx_95]

print(f"\nMaximum cumulative activity: {CumAct_values.max():.4f} GBq·s")
print(f"  Achieved at k_off = {k_off_values[np.argmax(CumAct_values)]:.4f} /min")
print(f"  (K_d = {K_d_values[np.argmax(CumAct_values)]:.2e} nmol/m³)")

print(f"\n95% of maximum cumulative activity achieved at:")
print(f"  k_off = {k_off_95:.4f} /min")
print(f"  K_d = {K_d_95:.2e} nmol/m³")
print(f"  Affinity ratio (k_int/k_off) = {k_int/k_off_95:.2f}")

print(f"\nDiminishing returns:")
print(f"  When k_off << k_int ({k_int:.4f} /min), internalization limits additional improvement")
print(f"  Further affinity improvements provide minimal benefit")

# =============================================================================
# PLOTTING
# =============================================================================

print("\nCreating plots...")

fig, axes = plt.subplots(1, 1)

# Plot 1: Cumulative activity vs k_off
ax1 = axes
ax1.semilogx(k_off_values, CumAct_values, 'o-', linewidth=2, markersize=6)
ax1.axvline(k_int, color='red', linestyle='--', linewidth=1.5, 
            label=f'$k_{int}$ = {k_int} /min', alpha=0.7)
#ax1.axvline(0.386, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
ax1.axvline(0.012, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
#ax1.axhline(CumAct_values.max() * 0.95, color='gray', linestyle=':', 
#            linewidth=1, alpha=0.5)
ax1.set_xlabel('$k_{off}$ (unbinding rate, min$^{-1}$)')
ax1.set_ylabel('Cumulative Activity (GBq·s)')
ax1.set_title('Off rate vs Cumulative Activity', fontweight='bold')
ax1.legend
ax1.grid(True, alpha=0.3)
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

# Add annotation
#ax1.annotate('Internalization\nlimiting', 
#             xy=(k_int/3, CumAct_values.max()*0.8),
#              ha='center', color='red', alpha=0.7)

'''
# Plot 2: Cumulative activity vs K_d (dissociation constant)
ax2 = axes[1]
ax2.loglog(K_d_values, CumAct_values, 'o-', linewidth=2, markersize=6, color='tab:orange')
ax2.set_xlabel('K_d (dissociation constant, nmol/m³)', fontsize=12)
ax2.set_ylabel('Cumulative Activity\n\lambda_decay∫(N_b^H + N_ic^H) dt (Bq·s)', fontsize=12)
ax2.set_title('Affinity (K_d) vs Cumulative Activity', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)

# Add annotation for "high affinity" and "low affinity" regions
ax2.annotate('High affinity\n(diminishing returns)', 
             xy=(K_d_values.min()*3, CumAct_values.max()*0.5),
             fontsize=10, ha='left', color='darkgreen', alpha=0.7)
ax2.annotate('Low affinity\n(poor targeting)', 
             xy=(K_d_values.max()/3, CumAct_values.min()*2),
             fontsize=10, ha='right', color='darkred', alpha=0.7)
'''
plt.tight_layout()

# Save figure
import os
output_dir = 'results/compare_models'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'koff_vs_CumAct.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_file}")

# =============================================================================
# DETAILED EXAMPLE: Compare high vs medium vs low affinity
# =============================================================================

print("\n" + "="*70)
print("DETAILED COMPARISON: Three Affinity Levels")
print("="*70)

# Select three representative k_off values
k_off_high = 0.0001      # Very high affinity (k_off << k_int)
k_off_medium = k_int     # Medium affinity (k_off = k_int)
k_off_low = 0.1          # Low affinity (k_off >> k_int)

fig2, axes2 = plt.subplots(3, 1, figsize=(10, 10))

for i, (k_off, label, color) in enumerate([
    (k_off_high, 'High affinity (k_off << $k_{int}$)', 'tab:green'),
    (k_off_medium, 'Medium affinity (k_off = $k_{int}$)', 'tab:orange'),
    (k_off_low, 'Low affinity (k_off >> $k_{int}$)', 'tab:red')
]):
    
    N_b_H, N_ic_H, CumAct = solve_pk_model(k_off)
    N_captive = N_b_H + N_ic_H
    
    ax = axes2[i]
    ax.plot(t_days, N_b_H, label='Bound (N_b^H)', linewidth=2, color=color, linestyle='--')
    ax.plot(t_days, N_ic_H, label='Internalized (N_ic^H)', linewidth=2, color=color, linestyle=':')
    ax.plot(t_days, N_captive, label='Captive (N_b^H + N_ic^H)', linewidth=2.5, color=color)
    
    ax.set_xlabel('Time (days)', fontsize=11)
    ax.set_ylabel('Amount (nmol)', fontsize=11)
    ax.set_title(f'{label}\nk_off = {k_off:.4f} min$^{-1}$, '
                 f'K_d = {k_off/k_on:.2e} nmol/m³, Cumulative activity = {CumAct:.4f} GBq·s',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 30])
    
    print(f"\n{label}:")
    print(f"  k_off = {k_off:.4f} /min")
    print(f"  K_d = {k_off/k_on:.2e} nmol/m³")
    print(f"  Peak N_b^H = {N_b_H.max():.4e} nmol")
    print(f"  Peak N_ic^H = {N_ic_H.max():.4e} nmol")
    print(f"  Cumulative activity = {CumAct:.4f} GBq·s")

plt.tight_layout()
output_file2 = os.path.join(output_dir, 'affinity_comparison_timecourses_koff.png')
plt.savefig(output_file2, dpi=300, bbox_inches='tight')
print(f"\nTime course figure saved to: {output_file2}")

print("\n" + "="*70)
print("STUDY COMPLETE")
print("="*70 + "\n")

plt.show()
