#!/usr/bin/env python3
"""
Plot beta particle escape fraction as a function of tumor radius

For Lu-177 with beta range ~1 mm, shows how much energy is deposited 
vs escapes for different tumor sizes.
"""

import numpy as np
import matplotlib.pyplot as plt
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
# PARAMETERS
# =============================================================================

r_beta = 1.0  # mm (Lu-177 mean beta range in tissue)

# Tumor radius range
r_tumor = np.linspace(0, 2.0, 1000)  # mm

# =============================================================================
# CALCULATE DEPOSITION AND ESCAPE FRACTIONS
# =============================================================================

# Deposition fraction (energy absorbed in tumor)
f_deposit = (r_tumor / (r_tumor + r_beta))**3

# Escape fraction (energy that escapes tumor)
f_escape = 1 - f_deposit

# =============================================================================
# PLOTTING
# =============================================================================

fig, ax = plt.subplots()

# Plot both fractions
ax.plot(r_tumor, f_deposit * 100, linewidth=2.5, label='Deposited in tumour', color='tab:blue')
ax.plot(r_tumor, f_escape * 100, linewidth=2.5, label='Escaped from tumour', color='tab:red')

# Mark key points
key_radii = [0.1, 0.333, 0.95]
for r in key_radii:
    f_dep = (r / (r + r_beta))**3
    ax.plot(r, f_dep * 100, 'o', markersize=8, color='tab:blue', zorder=5)
    ax.axvline(r, color='gray', linestyle=':', linewidth=1, alpha=0.5)

#    ax.annotate(f'({r:.1f} mm, {f_dep*100:.1f}%)',
#               xy=(r, f_dep * 100),
#               xytext=(r, f_dep * 100),
 #              fontsize=9,
#               arrowprops=dict(arrowstyle='->', color='black', lw=1),
#               ha='right')
    

# Formatting
ax.set_xlabel('Tumor Radius (mm)')
ax.set_ylabel('Energy Fraction (%)')
ax.set_title(f'Beta Particle Energy Deposition',fontweight='bold')
ax.legend(loc='right')
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 2.0])
ax.set_ylim([0, 100])

# Add shaded regions to show regimes
#ax.axvspan(0, 0.5, alpha=0.1, color='red', label='_nolegend_')
#ax.axvspan(0.5, 2.0, alpha=0.1, color='blue', label='_nolegend_')

#ax.text(0.25, 95, 'Escape-dominated\n(small tumors)', 
#        ha='center', va='top', fontsize=10, color='darkred', style='italic')
#ax.text(1.25, 95, 'Deposition-dominated\n(larger tumors)', 
#        ha='center', va='top', fontsize=10, color='darkblue', style='italic')

plt.tight_layout()

# Save figure
import os
output_dir = 'results/compare_models'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'beta_escape_vs_radius.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_file}")

# =============================================================================
# PRINT KEY VALUES
# =============================================================================

print("\n" + "="*60)
print("BETA ESCAPE ANALYSIS")
print("="*60)
print(f"\nLu-177 beta range: {r_beta} mm")
print("\nKey tumor sizes:")
print(f"{'Radius (mm)':<15} {'Deposited':<15} {'Escaped':<15}")
print("-" * 45)

for r in [0.1, 0.333, 0.95]:
    f_dep = (r / (r + r_beta))**3
    f_esc = 1 - f_dep
    print(f"{r:<15.1f} {f_dep*100:<14.2f}% {f_esc*100:<14.2f}%")

print("\n" + "="*60)
print("INTERPRETATION")
print("="*60)
print(f"• For R << r_beta ({r_beta} mm): Most energy escapes")
print(f"• For R = r_beta:  Only {12.5:.1f}% deposited")
print(f"• For R >> r_beta: Most energy deposited")
print(f"\nSmall tumors (R < 0.5 mm) are poorly treated by RPT due to escape!")
print("="*60 + "\n")

plt.show()
