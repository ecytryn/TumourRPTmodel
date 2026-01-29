#!/usr/bin/env python3
"""
Plot beta particle retention fraction as a function of tumor radius

Shows the uniform deposition sphere model compared to the old geometric approximation.
Uses similar matplotlib styling to pk_comparison_plot.py for consistency.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib as mpl

# Match styling from pk_comparison_plot.py
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

ELL = 1.0  # mm, mean range of Lu-177 beta particles

# =============================================================================
# BETA RETENTION MODEL
# =============================================================================

def sphere_intersection_volume(r, R, ell):
    """Volume of intersection between tumor and deposition sphere"""
    d = r
    
    if d >= R + ell:
        return 0.0
    
    if d <= abs(R - ell):
        return (4.0/3.0) * np.pi * min(R, ell)**3
    
    term1 = (R + ell - d)**2
    term2 = d**2 + 2*d*(ell + R) - 3*(R - ell)**2
    return (np.pi / (12*d)) * term1 * term2

def compute_retention_fraction(R, ell=ELL):
    """Compute retention fraction using uniform deposition sphere model"""
    V_ell = (4.0/3.0) * np.pi * ell**3
    
    def integrand(r):
        V_int = sphere_intersection_volume(r, R, ell)
        return (V_int / V_ell) * r**2
    
    integral, error = quad(integrand, 0, R, limit=100)
    f = (3.0 / R**3) * integral
    
    return f

def geometric_approximation(R, ell=ELL):
    """Old geometric approximation: R³/(R+ℓ)³"""
    return R**3 / (R + ell)**3

# =============================================================================
# GENERATE DATA
# =============================================================================

print("\n" + "="*70)
print("BETA RETENTION FRACTION PLOT")
print("="*70)
print(f"\nMean range: ℓ = {ELL} mm")

# Tumor radius range (clinically relevant: 0.5 - 3 mm)
R_values = np.linspace(0.1, 5.0, 100)

print(f"\nComputing retention fractions for {len(R_values)} tumor radii...")

# Compute retention fractions
retention_uniform = []
retention_geometric = []

for R in R_values:
    f_uniform = compute_retention_fraction(R, ELL)
    f_geom = geometric_approximation(R, ELL)
    
    retention_uniform.append(f_uniform)
    retention_geometric.append(f_geom)

retention_uniform = np.array(retention_uniform)
retention_geometric = np.array(retention_geometric)

print("  Computation complete")

# =============================================================================
# CREATE PLOT
# =============================================================================

print("\nCreating plot...")

fig, ax = plt.subplots(1, 1)

# Plot both models
ax.plot(R_values, retention_uniform, '-', 
       label='Uniform deposition sphere', linewidth=1.5, color='#2E86AB')
#ax.plot(R_values, retention_geometric, '--', 
#       label=r'Geometric $R^3/(R+\ell)^3$', linewidth=1.5, color='#A23B72')

# Mark mean range
ax.axvline(ELL, color='gray', linestyle=':', linewidth=0.8, alpha=0.6, zorder=1)
ax.text(ELL, 0.5, f'  $\ell$ = {ELL} mm', rotation=90, 
       va='center', fontsize=7, color='gray')

# Formatting
ax.set_xlabel('Tumor radius $R$ (mm)')
ax.set_ylabel('Fraction deposited in the tumour')
ax.set_title('Intratumoral energy deposition fraction', fontweight='bold')
#ax.legend(loc='lower right', framealpha=0.95)
ax.grid(True, alpha=0.3, linewidth=0.5)
ax.set_xlim([0, 5])
ax.set_ylim([0, 1.0])

# Add minor gridlines
ax.grid(True, which='minor', alpha=0.15, linewidth=0.3)
ax.minorticks_on()

plt.tight_layout()

# =============================================================================
# SAVE FIGURE
# =============================================================================

output_file = 'results/beta_retention_fraction.png'
import os
os.makedirs('results', exist_ok=True)

plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_file}")

# Also save as PDF for publication
output_pdf = 'results/beta_retention_fraction.pdf'
plt.savefig(output_pdf, bbox_inches='tight')
print(f"PDF saved to: {output_pdf}")

# =============================================================================
# PRINT COMPARISON STATISTICS
# =============================================================================

print("\n" + "="*70)
print("COMPARISON STATISTICS")
print("="*70)

# Relative error
relative_error = (retention_uniform - retention_geometric) / retention_geometric * 100

print(f"\nRelative error (Uniform - Geometric)/Geometric:")
print(f"  Mean:    {np.mean(np.abs(relative_error)):.1f}%")
print(f"  Maximum: {np.max(np.abs(relative_error)):.1f}%")
print(f"  At R=1mm: {relative_error[np.argmin(np.abs(R_values - 1.0))]:.1f}%")

# Sample values at key radii
print("\nRetention fractions at clinically relevant radii:")
print(f"{'R (mm)':<10} {'Uniform':<12} {'Geometric':<12} {'Difference':<10}")
print("-" * 50)

#for R_sample in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
for R_sample in [0.1, 0.333, 0.95]:
    idx = np.argmin(np.abs(R_values - R_sample))
    f_u = retention_uniform[idx]
    f_g = retention_geometric[idx]
    diff = (f_u - f_g) / f_g * 100
    print(f"{R_sample:<10.3f} {f_u:<12.4f} {f_g:<12.4f} {diff:<10.1f}%")

# Limiting behavior
print("\nLimiting behavior:")
print(f"  Small R (R=0.1mm): Uniform={retention_uniform[np.argmin(np.abs(R_values - 0.1))]:.6f}")
print(f"                     Theory: R³/ℓ³ = {(0.1/ELL)**3:.6f}")
print(f"  Large R (R=5.0mm): Uniform={retention_uniform[-1]:.6f}")
print(f"                     Theory: 1 - 3ℓ/R = {1 - 3*ELL/5.0:.6f}")

print("\n" + "="*70)
print("PLOT COMPLETE")
print("="*70 + "\n")

plt.show()
