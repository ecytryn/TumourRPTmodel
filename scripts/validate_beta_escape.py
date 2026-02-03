#!/usr/bin/env python3
"""
Beta Escape Fraction Validation

Computes actual 3D dose deposition in spherical tumor and compares
to the R³/(R+1)³ geometric escape fraction approximation.

For a tumor of radius R with uniform activity, computes:
1. Radial dose profile d(r) via Monte Carlo integration
2. Average dose in tumor
3. Compares to geometric approximation

Theory:
-------
Dose at position r from uniform source:
    d(r) = ∫∫∫ k(|r - r'|) dV'
    
where k(s) is the beta particle kernel (exponential decay)
    k(s) = (A/λ³) exp(-s/λ)

Geometric escape fraction (current model):
    f_escape = R³/(R+1)³
    
Actual escape fraction:
    f_escape = 1 - (average dose in tumor)/(dose if no escape)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import tplquad
from scipy.spatial.distance import cdist

# Lu-177 beta particle parameters
LAMBDA_BETA = 0.25  # mm, mean range of Lu-177 betas (~0.2-0.3 mm)
E_BETA_AVG = 134    # keV, average beta energy for Lu-177

def beta_kernel(s, lambda_range=LAMBDA_BETA):
    """
    Beta particle dose kernel - exponential approximation
    
    k(s) = (A/λ³) exp(-s/λ)
    
    Normalized so that ∫∫∫ k(s) dV = 1 over all space
    """
    A = 1.0 / (4 * np.pi * lambda_range**2)  # Normalization constant
    return A * np.exp(-s / lambda_range)

def dose_at_point_mc(r, R, lambda_range=LAMBDA_BETA, n_samples=100000):
    """
    Compute dose at radial position r in tumor of radius R
    using Monte Carlo integration
    
    d(r) = ∫∫∫_{|r'|<R} k(|r - r'|) dV'
    
    Parameters:
    -----------
    r : float
        Radial position (distance from center), in mm
    R : float
        Tumor radius in mm
    lambda_range : float
        Beta particle mean range in mm
    n_samples : int
        Number of Monte Carlo samples
        
    Returns:
    --------
    dose : float
        Dose at position r (arbitrary units)
    """
    # Generate random points uniformly in tumor sphere
    # Use rejection sampling in spherical coordinates for uniformity
    
    phi = np.random.uniform(0, 2*np.pi, n_samples)
    cos_theta = np.random.uniform(-1, 1, n_samples)
    theta = np.arccos(cos_theta)
    radius = R * np.random.uniform(0, 1, n_samples)**(1/3)  # r^(1/3) for uniform volume sampling
    
    # Convert to Cartesian (source points on z-axis, dose point at (r,0,0))
    x_source = radius * np.sin(theta) * np.cos(phi)
    y_source = radius * np.sin(theta) * np.sin(phi)
    z_source = radius * np.cos(theta)
    
    # Dose point at (r, 0, 0)
    x_dose = r
    y_dose = 0.0
    z_dose = 0.0
    
    # Compute distances
    distances = np.sqrt((x_source - x_dose)**2 + (y_source - y_dose)**2 + (z_source - z_dose)**2)
    
    # Apply beta kernel
    kernel_values = beta_kernel(distances, lambda_range)
    
    # Monte Carlo estimate: (Volume of tumor / N) * sum of kernel values
    tumor_volume = (4/3) * np.pi * R**3
    dose = (tumor_volume / n_samples) * np.sum(kernel_values)
    
    return dose

def compute_dose_profile(R, n_points=50, lambda_range=LAMBDA_BETA, n_samples=50000):
    """
    Compute radial dose profile d(r) for tumor of radius R
    
    Returns:
    --------
    r_values : array
        Radial positions from 0 to R
    doses : array
        Dose at each radial position
    """
    r_values = np.linspace(0, R, n_points)
    doses = []
    
    print(f"Computing dose profile for R={R:.2f} mm...")
    for i, r in enumerate(r_values):
        dose = dose_at_point_mc(r, R, lambda_range, n_samples)
        doses.append(dose)
        if (i+1) % 10 == 0:
            print(f"  Progress: {i+1}/{n_points} points")
    
    return r_values, np.array(doses)

def geometric_escape_fraction(R):
    """
    Current model's geometric escape fraction
    
    f_escape = R³/(R+1)³
    
    Note: R should be in units of lambda_beta for this formula
    """
    return R**3 / (R + 1)**3

def compute_average_dose(r_values, doses):
    """
    Compute volume-weighted average dose in tumor
    
    <d> = ∫₀^R d(r) · 4πr² dr / [(4/3)πR³]
        = 3/R³ ∫₀^R d(r) r² dr
    """
    # Use trapezoidal rule with volume weighting
    r_squared = r_values**2
    integrand = doses * r_squared
    
    avg_dose = 3 / r_values[-1]**3 * np.trapz(integrand, r_values)
    return avg_dose

def plot_comparison(R_values_mm, lambda_range=LAMBDA_BETA):
    """
    Create comprehensive comparison plot:
    1. Dose profiles for different tumor sizes
    2. Escape fraction vs tumor size
    3. Comparison to geometric approximation
    """
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Panel A: Dose profiles for different tumor sizes
    ax1 = axes[0, 0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(R_values_mm)))
    
    dose_profiles = []
    avg_doses = []
    center_doses = []
    
    for i, R_mm in enumerate(R_values_mm):
        r_vals, doses = compute_dose_profile(R_mm, n_points=30, lambda_range=lambda_range, n_samples=30000)
        dose_profiles.append((r_vals, doses))
        
        # Normalize to center dose for comparison
        center_dose = doses[0]
        center_doses.append(center_dose)
        normalized_doses = doses / center_dose
        
        avg_dose = compute_average_dose(r_vals, doses)
        avg_doses.append(avg_dose)
        
        # Plot normalized profile
        ax1.plot(r_vals, normalized_doses, color=colors[i], linewidth=2, 
                label=f'R={R_mm:.1f} mm')
    
    ax1.set_xlabel('Radial Position (mm)', fontsize=12)
    ax1.set_ylabel('Normalized Dose d(r)/d(0)', fontsize=12)
    ax1.set_title('A. Radial Dose Profiles', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1.1])
    
    # Panel B: Average dose vs tumor size
    ax2 = axes[0, 1]
    avg_doses_normalized = np.array(avg_doses) / np.array(center_doses)
    
    ax2.plot(R_values_mm, avg_doses_normalized, 'o-', linewidth=2, markersize=8, 
            label='Monte Carlo', color='darkblue')
    
    # Geometric approximation
    R_normalized = np.array(R_values_mm) / lambda_range
    geom_retention = geometric_escape_fraction(R_normalized)
    ax2.plot(R_values_mm, geom_retention, 's--', linewidth=2, markersize=6,
            label='Geometric: R³/(R+1)³', color='red', alpha=0.7)
    
    ax2.set_xlabel('Tumor Radius (mm)', fontsize=12)
    ax2.set_ylabel('Dose Retention Fraction', fontsize=12)
    ax2.set_title('B. Average Dose vs Tumor Size', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1.1])
    
    # Panel C: Escape fraction comparison
    ax3 = axes[1, 0]
    escape_mc = 1 - avg_doses_normalized
    escape_geom = 1 - geom_retention
    
    ax3.plot(R_values_mm, escape_mc * 100, 'o-', linewidth=2, markersize=8,
            label='Monte Carlo', color='darkblue')
    ax3.plot(R_values_mm, escape_geom * 100, 's--', linewidth=2, markersize=6,
            label='Geometric R³/(R+1)³', color='red', alpha=0.7)
    
    ax3.set_xlabel('Tumor Radius (mm)', fontsize=12)
    ax3.set_ylabel('Escape Fraction (%)', fontsize=12)
    ax3.set_title('C. Beta Escape Fraction', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # Panel D: Relative error
    ax4 = axes[1, 1]
    relative_error = (escape_geom - escape_mc) / escape_mc * 100
    
    ax4.plot(R_values_mm, relative_error, 'o-', linewidth=2, markersize=8, color='darkgreen')
    ax4.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax4.fill_between(R_values_mm, -10, 10, alpha=0.2, color='green')
    
    ax4.set_xlabel('Tumor Radius (mm)', fontsize=12)
    ax4.set_ylabel('Relative Error (%)', fontsize=12)
    ax4.set_title('D. Geometric Approximation Error', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.text(0.95, 0.95, '±10% band', transform=ax4.transAxes,
            ha='right', va='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    return fig, (avg_doses_normalized, geom_retention, escape_mc, escape_geom)

def main():
    """Main analysis"""
    
    print("=" * 60)
    print("Beta Escape Fraction Validation")
    print("=" * 60)
    print(f"Lu-177 beta range: {LAMBDA_BETA} mm")
    print()
    
    # Test range of tumor sizes
    R_values_mm = np.array([0.5, 1.0, 2.0, 3.0, 5.0])  # mm
    
    print("Tumor sizes to test:")
    for R in R_values_mm:
        R_norm = R / LAMBDA_BETA
        print(f"  R = {R:.1f} mm ({R_norm:.1f} λ)")
    print()
    
    # Generate comparison plots
    fig, (retention_mc, retention_geom, escape_mc, escape_geom) = plot_comparison(R_values_mm)
    
    # Print numerical results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"{'R (mm)':<10} {'Retention (MC)':<18} {'Retention (Geom)':<18} {'Error (%)':<10}")
    print("-" * 60)
    
    for i, R in enumerate(R_values_mm):
        error = (retention_geom[i] - retention_mc[i]) / retention_mc[i] * 100
        print(f"{R:<10.1f} {retention_mc[i]:<18.3f} {retention_geom[i]:<18.3f} {error:<10.1f}")
    
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)
    
    # Assess approximation quality
    errors = np.abs((retention_geom - retention_mc) / retention_mc * 100)
    max_error = np.max(errors)
    mean_error = np.mean(errors)
    
    print(f"Maximum relative error: {max_error:.1f}%")
    print(f"Mean relative error: {mean_error:.1f}%")
    
    if max_error < 10:
        print("✓ Geometric approximation is GOOD (< 10% error)")
    elif max_error < 20:
        print("⚠ Geometric approximation is ACCEPTABLE (10-20% error)")
    else:
        print("✗ Geometric approximation has SIGNIFICANT error (> 20%)")
    
    # Save figure
    output_file = "beta_escape_validation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
