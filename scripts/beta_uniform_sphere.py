#!/usr/bin/env python3
"""
Beta Retention - Uniform Deposition Sphere Method

Physical model:
- Each decay deposits energy uniformly in sphere of radius ℓ (mean free path)
- Computes fraction deposited inside tumor by integrating sphere-sphere overlaps

f_deposit = (1/V_tumor) ∫_tumor [V_int(y) / V_ℓ] dy

where V_int(y) is the volume of intersection between:
  - Tumor sphere (radius R)
  - Deposition sphere (radius ℓ, centered at y)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def sphere_intersection_volume(r, R, ell):
    """
    Volume of intersection between two spheres
    
    Sphere 1: Radius R, centered at origin (tumor)
    Sphere 2: Radius ℓ, centered at distance r from origin (deposition)
    
    Parameters:
    -----------
    r : float
        Distance from tumor center to deposition center
    R : float
        Tumor radius
    ell : float
        Deposition sphere radius (mean free path)
        
    Returns:
    --------
    V_int : float
        Intersection volume
    """
    d = r  # Distance between centers
    
    # Case 1: No overlap
    if d >= R + ell:
        return 0.0
    
    # Case 2: One sphere completely inside the other
    if d <= abs(R - ell):
        return (4.0/3.0) * np.pi * min(R, ell)**3
    
    # Case 3: Partial overlap
    # Formula: V = π/(12d) * (R+ℓ-d)² * [d² + 2d(ℓ+R) - 3(R-ℓ)²]
    term1 = (R + ell - d)**2
    term2 = d**2 + 2*d*(ell + R) - 3*(R - ell)**2
    V_int = (np.pi / (12*d)) * term1 * term2
    
    return V_int

def compute_retention_fraction(R, ell, n_points=100):
    """
    Compute beta retention fraction using uniform deposition sphere model
    
    f = (3/R³) ∫₀ᴿ [V_int(r) / V_ℓ] r² dr
    
    Parameters:
    -----------
    R : float
        Tumor radius (mm)
    ell : float
        Mean free path / deposition radius (mm)
        
    Returns:
    --------
    f : float
        Retention fraction (0 to 1)
    """
    V_tumor = (4.0/3.0) * np.pi * R**3
    V_ell = (4.0/3.0) * np.pi * ell**3
    
    # Integrand: [V_int(r) / V_ℓ] * r²
    def integrand(r):
        V_int = sphere_intersection_volume(r, R, ell)
        return (V_int / V_ell) * r**2
    
    # Integrate from 0 to R
    integral, error = quad(integrand, 0, R, limit=n_points)
    
    # f = (3/R³) * integral
    f = (3.0 / R**3) * integral
    
    return f

def geometric_retention(R, ell):
    """
    Geometric approximation: R³/(R+ℓ)³
    """
    return R**3 / (R + ell)**3

def plot_comparison(R_values, ell, title_suffix=""):
    """
    Compare uniform deposition model to geometric approximation
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Compute retention for each tumor size
    retentions_uniform = []
    retentions_geom = []
    
    print(f"\nMean free path ℓ = {ell:.2f} mm")
    print(f"\n{'R (mm)':<10} {'R/ℓ':<8} {'Uniform':<12} {'Geometric':<12} {'Diff (%)':<10}")
    print("-" * 60)
    
    for R in R_values:
        # Uniform deposition model
        f_uniform = compute_retention_fraction(R, ell)
        retentions_uniform.append(f_uniform)
        
        # Geometric approximation
        f_geom = geometric_retention(R, ell)
        retentions_geom.append(f_geom)
        
        # Difference
        diff = (f_uniform - f_geom) / f_geom * 100 if f_geom > 0 else 0
        
        print(f"{R:<10.2f} {R/ell:<8.2f} {f_uniform:<12.4f} {f_geom:<12.4f} {diff:<10.2f}")
    
    retentions_uniform = np.array(retentions_uniform)
    retentions_geom = np.array(retentions_geom)
    
    # Panel 1: Retention vs Tumor Size
    ax = axes[0]
    
    ax.plot(R_values, retentions_uniform, 'o-', linewidth=2.5, markersize=10,
           label='Uniform Deposition Sphere', color='darkblue')
    ax.plot(R_values, retentions_geom, 's--', linewidth=2.5, markersize=8,
           label='Geometric R³/(R+ℓ)³', color='red', alpha=0.7)
    
    # Fine geometric curve
    R_fine = np.linspace(0.05, max(R_values), 200)
    f_geom_fine = geometric_retention(R_fine, ell)
    ax.plot(R_fine, f_geom_fine, '--', linewidth=1.5, color='red', alpha=0.3)
    
    ax.set_xlabel('Tumor Radius R (mm)', fontsize=12)
    ax.set_ylabel('Retention Fraction', fontsize=12)
    ax.set_title(f'Beta Retention vs Tumor Size\n(ℓ = {ell} mm){title_suffix}', 
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    ax.set_xlim([0, max(R_values)])
    
    # Mark ℓ
    ax.axvline(ell, color='gray', linestyle=':', alpha=0.5, linewidth=2)
    ax.text(ell, 0.5, f'  ℓ = {ell} mm', rotation=90,
           va='center', fontsize=10, color='gray')
    
    # Panel 2: Relative Difference
    ax = axes[1]
    
    relative_diff = (retentions_uniform - retentions_geom) / retentions_geom * 100
    
    ax.plot(R_values, relative_diff, 'o-', linewidth=2.5, markersize=10,
           color='darkgreen')
    ax.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.fill_between(R_values, -10, 10, alpha=0.2, color='green', label='±10%')
    ax.fill_between(R_values, -5, 5, alpha=0.3, color='green', label='±5%')
    
    ax.set_xlabel('Tumor Radius R (mm)', fontsize=12)
    ax.set_ylabel('Relative Difference (%)', fontsize=12)
    ax.set_title('Geometric Approximation Error\n(Uniform - Geom)/Geom × 100%',
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, max(R_values)])
    
    plt.tight_layout()
    
    # Summary statistics
    max_error = np.max(np.abs(relative_diff))
    mean_error = np.mean(np.abs(relative_diff))
    
    print(f"\n{'='*60}")
    print(f"SUMMARY:")
    print(f"  Maximum error: {max_error:.2f}%")
    print(f"  Mean error: {mean_error:.2f}%")
    
    if max_error < 10:
        print(f"  ✓ Geometric approximation is EXCELLENT (< 10% error)")
    elif max_error < 20:
        print(f"  ✓ Geometric approximation is GOOD (< 20% error)")
    else:
        print(f"  ⚠ Geometric approximation has significant error (> 20%)")
    
    return fig, (retentions_uniform, retentions_geom)

def test_different_ell_values():
    """
    Test multiple values of ℓ to see which gives best match
    """
    R_values = np.array([0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    ell_values = [0.2, 0.25, 0.3, 0.5, 1.0]
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, ell in enumerate(ell_values):
        ax = axes[idx]
        
        retentions_uniform = [compute_retention_fraction(R, ell) for R in R_values]
        retentions_geom = [geometric_retention(R, ell) for R in R_values]
        
        relative_diff = (np.array(retentions_uniform) - np.array(retentions_geom)) / np.array(retentions_geom) * 100
        max_error = np.max(np.abs(relative_diff))
        
        ax.plot(R_values, retentions_uniform, 'o-', linewidth=2, label='Uniform Sphere')
        ax.plot(R_values, retentions_geom, 's--', linewidth=2, label='Geometric', alpha=0.7)
        
        ax.set_xlabel('R (mm)', fontsize=10)
        ax.set_ylabel('Retention', fontsize=10)
        ax.set_title(f'ℓ = {ell} mm (max error: {max_error:.1f}%)', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1.05])
    
    axes[-1].axis('off')
    
    plt.tight_layout()
    return fig

def main():
    """Main execution"""
    
    print("="*60)
    print("BETA RETENTION - UNIFORM DEPOSITION SPHERE METHOD")
    print("="*60)
    
    # Test with different ℓ values
    print("\n1. Testing different mean free path values:")
    fig_test = test_different_ell_values()
    plt.savefig("beta_retention_ell_comparison.png", dpi=300, bbox_inches='tight')
    print("   Saved: beta_retention_ell_comparison.png")
    
    # Detailed analysis with ℓ = 1.0 mm (your model value)
    print("\n2. Detailed analysis with ℓ = 1.0 mm:")
    R_values = np.linspace(0.1, 3.0, 30)
    fig_main, (f_uniform, f_geom) = plot_comparison(R_values, ell=1.0, title_suffix=" (Your Model)")
    plt.savefig("beta_retention_uniform_sphere.png", dpi=300, bbox_inches='tight')
    print("   Saved: beta_retention_uniform_sphere.png")
    
    print(f"\n{'='*60}")
    print("INTERPRETATION:")
    print(f"{'='*60}")
    print("The uniform deposition sphere model:")
    print("  • Assumes each decay deposits uniformly in sphere of radius ℓ")
    print("  • Computes geometric overlap between tumor and deposition spheres")
    print("  • More physically motivated than exponential kernel")
    print("")
    print("Comparison to R³/(R+ℓ)³:")
    print("  • Shows how well the simple geometric formula approximates")
    print("  • Can tune ℓ to match your simulation results")
    print("  • Both capture correct physics (small R → escape, large R → retain)")
    
    plt.show()

if __name__ == "__main__":
    main()
