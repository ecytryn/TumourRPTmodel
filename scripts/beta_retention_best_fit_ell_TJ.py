#!/usr/bin/env python3
"""
Fit the mean range parameter ELL of the uniform deposition sphere model
to TJ's beta retention lookup table data.

Fits over R = 0 to 2 mm, including the anchor point (0, 0).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import matplotlib as mpl
import os

mpl.rcParams.update({
    "figure.figsize": (3.35, 2.4),
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# =============================================================================
# TJ LOOKUP TABLE (from BetaRetention.java) + anchor point (0, 0)
# =============================================================================

R_VALUES_TJ = np.array([
    0.0000,
    0.1000, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000,
    0.7000, 0.8000, 0.9000, 1.0000, 1.2000, 1.4000,
    1.6000, 1.8000, 2.0000, 2.5000, 3.0000, 3.5000,
    4.0000, 4.5000, 5.0000
])

RETENTION_TJ = np.array([
    0.0000,
    0.1982, 0.2899, 0.3641, 0.4275, 0.4861, 0.5408,
    0.5921, 0.6420, 0.6898, 0.7342, 0.8137, 0.8760,
    0.9187, 0.9449, 0.9593, 0.9701, 0.9700, 0.9695,
    0.9700, 0.9714, 0.9734
])

# Restrict to fitting range 0 to 2 mm
FIT_MASK = R_VALUES_TJ <= 3.0
R_FIT = R_VALUES_TJ[FIT_MASK]
F_FIT = RETENTION_TJ[FIT_MASK]

# =============================================================================
# UNIFORM DEPOSITION SPHERE MODEL
# =============================================================================

def sphere_intersection_volume(r, R, ell):
    d = r
    if d >= R + ell:
        return 0.0
    if d <= abs(R - ell):
        return (4.0/3.0) * np.pi * min(R, ell)**3
    term1 = (R + ell - d)**2
    term2 = d**2 + 2*d*(ell + R) - 3*(R - ell)**2
    return (np.pi / (12*d)) * term1 * term2

def compute_retention_fraction(R, ell):
    """Retention fraction for a single R value."""
    if R == 0.0:
        return 0.0
    V_ell = (4.0/3.0) * np.pi * ell**3
    def integrand(r):
        return (sphere_intersection_volume(r, R, ell) / V_ell) * r**2
    integral, _ = quad(integrand, 0, R, limit=100)
    return (3.0 / R**3) * integral

def compute_retention_curve(R_array, ell):
    return np.array([compute_retention_fraction(R, ell) for R in R_array])

# =============================================================================
# FIT
# =============================================================================

def residual_sse(ell):
    """Sum of squared errors over fitting range."""
    predicted = compute_retention_curve(R_FIT, ell)
    return np.sum((predicted - F_FIT)**2)

print("Fitting ELL to TJ table data (R = 0 to 2 mm)...")
result = minimize_scalar(residual_sse, bounds=(0.05, 5.0), method='bounded')
ell_best = result.x
sse_best = result.fun

print(f"\nBest fit ELL = {ell_best:.4f} mm")
print(f"SSE          = {sse_best:.6f}")
print(f"RMSE         = {np.sqrt(sse_best / len(R_FIT)):.4f}")

# Compare with original ELL = 0.125 mm
sse_original = residual_sse(0.125)
print(f"\nFor comparison, ELL = 0.125 mm gives SSE = {sse_original:.6f}")

# =============================================================================
# PLOT
# =============================================================================

R_plot = np.linspace(0.001, 5.0, 300)
#R_plot = np.logspace(-2, np.log10(5.0), 300)
f_best  = compute_retention_curve(R_plot, ell_best)
f_orig  = compute_retention_curve(R_plot, 0.125)

fig, ax = plt.subplots()

#ax.plot(R_plot, f_orig, '-', color='#2E86AB', linewidth=1.5,
#        label=f'Uniform sphere ($\\ell$ = 0.125 mm)')
ax.plot(R_plot, R_plot, '--', color='#2E86AB', linewidth=1)
ax.plot(R_plot, R_plot**3, '--', color='#2E86AB', linewidth=1)

ax.plot(R_plot, f_best, '-', color='#2E86AB', linewidth=1.5,
        label=f'Best fit ($\\ell$ = {ell_best:.3f} mm)')
ax.plot(R_VALUES_TJ, RETENTION_TJ, 'o', color='#E84855',
        markersize=3.5, markeredgewidth=0.5, label='Monte Carlo model')
#ax.axvline(2.0, color='gray', linewidth=0.8, linestyle=':', alpha=0.6)
ax.text(2.05, 0.1, 'fit limit', fontsize=6, color='gray')

ax.set_xlabel('Tumour radius $R$ (mm)')
ax.set_ylabel('Fraction deposited in tumour')
ax.set_title('Best fit Uniform Sphere EDF to Monte Carlo EDF', fontweight='bold')
ax.legend(loc='lower right', framealpha=0.95)
ax.grid(True, alpha=0.3, linewidth=0.5)
ax.grid(True, which='minor', alpha=0.15, linewidth=0.3)
ax.minorticks_on()
#ax.set_xlim([0, 3.0])
#ax.set_ylim([0, 1.0])
ax.set_xlim([0.05, 5.0])
ax.set_ylim([0.01, 1.0])

ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()

os.makedirs('results', exist_ok=True)
plt.savefig('results/beta_retention_ell_fit.pdf', bbox_inches='tight')
plt.savefig('results/beta_retention_ell_fit.png', dpi=300, bbox_inches='tight')
print("\nSaved to results/beta_retention_ell_fit.pdf and .png")

plt.show()
