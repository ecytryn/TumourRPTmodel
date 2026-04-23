#!/usr/bin/env python3
"""
Beta Escape Fraction Validation - CORRECTED

Computes actual beta retention fraction and compares to R³/(R+1)³ approximation.

Key calculation:
1. Compute d(r) = dose at radius r from all sources in sphere
2. Integrate over tumor: E_in = ∫₀ᴿ d(r) · 4πr² dr
3. Total emission: E_total = (4π/3)R³ · ε (uniform source)
4. Retention: f = E_in / E_total

This should match the geometric approximation R³/(R+1)³ where R is in units of λ.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import time

# Lu-177 beta particle parameters
LAMBDA_BETA = 0.25  # mm, mean range
EPSILON = 1.0       # Energy per decay (arbitrary units, cancels out in ratio)

def beta_kernel(s, lambda_range=LAMBDA_BETA):
    """
    Beta particle dose kernel - exponential approximation
    
    k(s) = (A/λ³) exp(-s/λ)
    
    Represents energy deposited at distance s from a point source
    """
    A = 1.0 / (4 * np.pi * lambda_range**2)
    return A * np.exp(-s / lambda_range) * EPSILON

def dose_at_point_cylindrical(r, R, lambda_range=LAMBDA_BETA):
    """
    Compute dose at radius r using cylindrical symmetry
    
    d(r) = ∫∫ k(√[(x-r)² + z²]) · 2πz · ρ₀ dx dz
    
    where ρ₀ = 1 is the uniform activity density in the tumor
    
    Integration domain: Quarter-circle in x-z plane
    """
    def integrand(z, x):
        distance = np.sqrt((x - r)**2 + z**2)
        kernel_val = beta_kernel(distance, lambda_range)
        return kernel_val * 2 * np.pi * z  # Revolution factor
    
    result, error = dblquad(
        integrand,
        -R, R,                           # x limits
        lambda x: 0,                     # z lower
        lambda x: np.sqrt(R**2 - x**2)   # z upper (circle)
    )
    
    return result, error

def compute_retention_fraction(R, n_radial_points=30, lambda_range=LAMBDA_BETA):
    """
    Compute the fraction of emitted energy that is deposited back in the tumor
    
    Returns:
    --------
    retention : float
        Fraction of energy retained (0 to 1)
    r_values : array
        Radial positions used
    doses : array
        Dose d(r) at each radius
    """
    print(f"\nComputing retention for R = {R:.2f} mm...")
    
    # Step 1: Compute d(r) at various radii
    r_values = np.linspace(0, R, n_radial_points)
    doses = []
    
    for i, r in enumerate(r_values):
        dose, error = dose_at_point_cylindrical(r, R, lambda_range)
        doses.append(dose)
        
        if (i+1) % 10 == 0:
            print(f"  Progress: {i+1}/{n_radial_points} radial points")
    
    doses = np.array(doses)
    
    # Step 2: Integrate d(r) over tumor volume using spherical shells
    # E_in = ∫₀ᴿ d(r) · 4πr² dr
    
    # Integrand: d(r) · r²
    integrand = doses * r_values**2
    
    # Integrate using trapezoidal rule
    E_in_tumor = 4 * np.pi * np.trapz(integrand, r_values)
    
    # Step 3: Total energy emitted from tumor
    # E_total = (4π/3) R³ · ρ₀ · ε
    # where ρ₀ = 1 (uniform activity density), ε = energy per decay
    
    E_total = (4.0/3.0) * np.pi * R**3 * EPSILON
    
    # Step 4: Retention fraction
    retention = E_in_tumor / E_total
    
    print(f"  E_in_tumor = {E_in_tumor:.6e}")
    print(f"  E_total = {E_total:.6e}")
    print(f"  Retention = {retention:.6f}")
    
    return retention, r_values, doses

def geometric_retention(R, lambda_range=LAMBDA_BETA):
    """
    Geometric retention approximation: R³/(R+1)³
    
    where R is in units of λ (normalized radius)
    """
    R_normalized = R / lambda_range
    return R_normalized**3 / (R_normalized + 1)**3

def plot_validation(R_values_mm, lambda_range=LAMBDA_BETA):
    """
    Create comprehensive validation figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(R_values_mm)))
    
    results = []
    
    # Compute for each tumor size
    for i, R in enumerate(R_values_mm):
        print(f"\n{'='*70}")
        print(f"TUMOR SIZE: R = {R:.2f} mm (R/λ = {R/lambda_range:.2f})")
        print(f"{'='*70}")
        
        retention, r_vals, doses = compute_retention_fraction(R, n_radial_points=40)
        
        # Geometric approximation
        retention_geom = geometric_retention(R, lambda_range)
        
        results.append({
            'R': R,
            'retention': retention,
            'retention_geom': retention_geom,
            'r_vals': r_vals,
            'doses': doses
        })
        
        print(f"  Cylindrical retention: {retention:.6f}")
        print(f"  Geometric R³/(R+1)³:   {retention_geom:.6f}")
        print(f"  Relative error: {abs(retention - retention_geom)/retention_geom * 100:.2f}%")
        
        # Panel A: Normalized dose profiles
        ax = axes[0, 0]
        normalized_doses = doses / doses[0] if doses[0] > 0 else doses
        ax.plot(r_vals, normalized_doses, color=colors[i], linewidth=2, 
               label=f'R={R:.2f} mm')
    
    # Finalize Panel A
    ax = axes[0, 0]
    ax.set_xlabel('Radial Position r (mm)', fontsize=12)
    ax.set_ylabel('Normalized Dose d(r)/d(0)', fontsize=12)
    ax.set_title('A. Radial Dose Profiles', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])
    
    # Panel B: Retention vs Tumor Size
    ax = axes[0, 1]
    
    R_array = np.array([res['R'] for res in results])
    ret_array = np.array([res['retention'] for res in results])
    ret_geom_array = np.array([res['retention_geom'] for res in results])
    
    ax.plot(R_array, ret_array, 'o-', linewidth=2.5, markersize=10,
           label='Cylindrical Integration', color='darkblue')
    ax.plot(R_array, ret_geom_array, 's--', linewidth=2.5, markersize=8,
           label='Geometric R³/(R+1)³', color='red', alpha=0.7)
    
    # Add fine geometric curve
    R_fine = np.linspace(0.05, max(R_values_mm), 200)
    ret_geom_fine = geometric_retention(R_fine, lambda_range)
    ax.plot(R_fine, ret_geom_fine, '--', linewidth=1.5, color='red', alpha=0.3)
    
    ax.set_xlabel('Tumor Radius R (mm)', fontsize=12)
    ax.set_ylabel('Retention Fraction', fontsize=12)
    ax.set_title('B. Beta Retention vs Tumor Size', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    ax.set_xlim([0, max(R_values_mm)])
    
    # Add vertical line at λ
    ax.axvline(lambda_range, color='gray', linestyle=':', alpha=0.5, linewidth=2)
    ax.text(lambda_range, 0.5, f'  λ = {lambda_range} mm', rotation=90,
           va='center', fontsize=10, color='gray')
    
    # Panel C: Relative Error
    ax = axes[1, 0]
    
    relative_error = (ret_array - ret_geom_array) / ret_geom_array * 100
    
    ax.plot(R_array, relative_error, 'o-', linewidth=2.5, markersize=10,
           color='darkgreen')
    ax.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.fill_between(R_array, -10, 10, alpha=0.2, color='green', 
                    label='±10% band')
    ax.fill_between(R_array, -5, 5, alpha=0.3, color='green',
                    label='±5% band')
    
    ax.set_xlabel('Tumor Radius R (mm)', fontsize=12)
    ax.set_ylabel('Relative Error (%)', fontsize=12)
    ax.set_title('C. Geometric Approximation Error\n(Cyl - Geom)/Geom × 100%', 
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Panel D: Summary Table
    ax = axes[1, 1]
    ax.axis('off')
    
    # Create summary text
    summary_lines = [
        "VALIDATION SUMMARY",
        "="*50,
        f"Beta particle range (λ): {lambda_range} mm",
        f"",
        "Retention Fraction Results:",
        "",
        f"{'R (mm)':<10} {'R/λ':<8} {'Cylindrical':<12} {'Geometric':<12} {'Error (%)':<10}",
        "-"*50
    ]
    
    for res in results:
        R = res['R']
        R_norm = R / lambda_range
        ret_cyl = res['retention']
        ret_geom = res['retention_geom']
        error = (ret_cyl - ret_geom) / ret_geom * 100
        
        summary_lines.append(
            f"{R:<10.2f} {R_norm:<8.2f} {ret_cyl:<12.4f} {ret_geom:<12.4f} {error:<10.2f}"
        )
    
    summary_lines.extend([
        "",
        "="*50,
        "CONCLUSIONS:",
        "",
    ])
    
    max_error = np.max(np.abs(relative_error))
    mean_error = np.mean(np.abs(relative_error))
    
    summary_lines.append(f"Maximum error: {max_error:.2f}%")
    summary_lines.append(f"Mean error: {mean_error:.2f}%")
    summary_lines.append("")
    
    if max_error < 10:
        summary_lines.append("✓ Geometric approximation is EXCELLENT")
        summary_lines.append("  (< 10% error across all sizes)")
    elif max_error < 20:
        summary_lines.append("✓ Geometric approximation is GOOD")
        summary_lines.append("  (< 20% error)")
    else:
        summary_lines.append("⚠ Geometric approximation has errors")
        summary_lines.append("  (> 20% in some regimes)")
    
    summary_lines.extend([
        "",
        "Physical interpretation:",
        f"• Small R (< λ): Most betas escape",
        f"• Large R (> λ): Most betas retained",
        f"• Transition around R ≈ λ = {lambda_range} mm"
    ])
    
    summary_text = "\n".join(summary_lines)
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top', family='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
    
    plt.tight_layout()
    return fig, results

def main():
    """Main validation"""
    
    print("="*70)
    print("BETA RETENTION FRACTION VALIDATION")
    print("="*70)
    print("\nValidating geometric approximation: R³/(R+1)³")
    print(f"Lu-177 beta range: λ = {LAMBDA_BETA} mm")
    print("\nComputing energy retention for various tumor sizes...")
    
    # Test range of tumor sizes (in mm)
    # Include sizes both smaller and larger than λ
    R_values_mm = np.array([0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    
    print(f"\nTumor sizes to test (λ = {LAMBDA_BETA} mm):")
    for R in R_values_mm:
        print(f"  R = {R:.2f} mm (R/λ = {R/LAMBDA_BETA:.2f})")
    
    # Run validation
    fig, results = plot_validation(R_values_mm, LAMBDA_BETA)
    
    # Save figure
    output_file = "beta_retention_validation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    print(f"\n{'='*70}")
    print(f"Figure saved: {output_file}")
    print(f"{'='*70}")
    
    plt.show()

if __name__ == "__main__":
    main()
