#!/usr/bin/env python3
"""
Plot beta particle retention fraction as a function of tumor radius

Compares:
  1. Uniform deposition sphere model (continuous curve)
  2. TJ lookup table values (discrete markers)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib as mpl
import os

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
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# =============================================================================
# CONFIGURATION
# =============================================================================

ELL = 0.3  # mm, mean range of Lu-177 beta particles

# TJ lookup table (from BetaRetention.java)
R_VALUES_TJ = np.array([
    0.1000, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000,
    0.7000, 0.8000, 0.9000, 1.0000, 1.2000, 1.4000,
    1.6000, 1.8000, 2.0000, 2.5000, 3.0000, 3.5000,
    4.0000, 4.5000, 5.0000
])

RETENTION_TJ = np.array([
    0.1982, 0.2899, 0.3641, 0.4275, 0.4861, 0.5408,
    0.5921, 0.6420, 0.6898, 0.7342, 0.8137, 0.8760,
    0.9187, 0.9449, 0.9593, 0.9701, 0.9700, 0.9695,
    0.9700, 0.9714, 0.9734
])

# =============================================================================
# UNIFORM DEPOSITION SPHERE MODEL
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
    integral, _ = quad(integrand, 0, R, limit=100)
    return (3.0 / R**3) * integral

# =============================================================================
# GENERATE DATA
# =============================================================================

R_values = np.linspace(0.01, 5.0, 300)
retention_uniform = np.array([compute_retention_fraction(R) for R in R_values])

# =============================================================================
# PLOT
# =============================================================================

fig, ax = plt.subplots(1, 1)

ax.plot(R_values, retention_uniform, '-',
        label='Uniform deposition sphere', linewidth=1.5, color='#2E86AB')

ax.plot(R_VALUES_TJ, RETENTION_TJ, 'o',
        label='TJ table', markersize=3.5, color='#E84855',
        markerfacecolor='#E84855', markeredgewidth=0.5)

ax.set_xlabel('Tumour radius $R$ (mm)')
ax.set_ylabel('Fraction deposited in tumour')
ax.set_title('Intratumoral energy deposition fraction', fontweight='bold')
ax.legend(loc='lower right', framealpha=0.95)
ax.grid(True, alpha=0.3, linewidth=0.5)
ax.grid(True, which='minor', alpha=0.15, linewidth=0.3)
ax.minorticks_on()
ax.set_xlim([0, 5.0])
ax.set_ylim([0, 1.0])

plt.tight_layout()

# =============================================================================
# SAVE
# =============================================================================

os.makedirs('results', exist_ok=True)
plt.savefig('results/beta_retention_fraction.png', dpi=300, bbox_inches='tight')
plt.savefig('results/beta_retention_fraction.pdf', bbox_inches='tight')
print("Saved to results/beta_retention_fraction.png and .pdf")

plt.show()
