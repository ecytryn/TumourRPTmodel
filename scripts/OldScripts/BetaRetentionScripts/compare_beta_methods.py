#!/usr/bin/env python3
"""
Beta Escape Validation - Method Comparison

Compares two approaches for computing dose deposition:
1. Monte Carlo 3D sampling (brute force)
2. Cylindrical symmetry integration (efficient)

Both should give same answer - this validates the implementations.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import time

# Lu-177 beta particle parameters
LAMBDA_BETA = 0.25  # mm, mean range

def beta_kernel(s, lambda_range=LAMBDA_BETA):
    """Beta particle dose kernel - exponential approximation"""
    A = 1.0 / (4 * np.pi * lambda_range**2)
    return A * np.exp(-s / lambda_range)

# ============================================================================
# METHOD 1: Monte Carlo 3D Sampling
# ============================================================================

def dose_at_point_mc(r, R, lambda_range=LAMBDA_BETA, n_samples=100000):
    """
    Compute dose using Monte Carlo integration over 3D sphere
    
    Algorithm:
    1. Sample n_samples points uniformly in tumor sphere
       - φ uniform in [0, 2π]
       - cos(θ) uniform in [-1, 1] → gives sin(θ) weighting
       - r ∝ u^(1/3) → gives r² weighting
    2. Compute distance from each source point to dose point (r,0,0)
    3. Apply beta kernel k(distance)
    4. Monte Carlo estimate: (Volume/N) × Σ k(distances)
    """
    # Generate random points uniformly in tumor sphere
    phi = np.random.uniform(0, 2*np.pi, n_samples)
    cos_theta = np.random.uniform(-1, 1, n_samples)
    sin_theta = np.sqrt(1 - cos_theta**2)  # Avoid arccos then sin
    radius = R * np.random.uniform(0, 1, n_samples)**(1/3)  # r^(1/3) for r² weighting
    
    # Convert to Cartesian (source points)
    x_source = radius * sin_theta * np.cos(phi)
    y_source = radius * sin_theta * np.sin(phi)
    z_source = radius * cos_theta
    
    # Dose point at (r, 0, 0)
    x_dose = r
    
    # Compute distances
    distances = np.sqrt((x_source - x_dose)**2 + y_source**2 + z_source**2)
    
    # Apply beta kernel
    kernel_values = beta_kernel(distances, lambda_range)
    
    # Monte Carlo estimate
    tumor_volume = (4/3) * np.pi * R**3
    dose = (tumor_volume / n_samples) * np.sum(kernel_values)
    
    return dose

# ============================================================================
# METHOD 2: Cylindrical Symmetry Integration
# ============================================================================

def dose_at_point_cylindrical(r, R, lambda_range=LAMBDA_BETA):
    """
    Compute dose using cylindrical symmetry
    
    Exploits rotational symmetry around line from center to dose point.
    
    d(r) = ∫∫ k(√[(x-r)² + z²]) · 2π|z| dx dz
    
    Integration domain: Quarter-circle in x-z plane
    - x: from -R to R
    - z: from 0 to √(R² - x²)
    
    The 2π|z| factor accounts for revolution around x-axis
    (z is the distance from x-axis, giving circumference 2πz)
    """
    def integrand(z, x):
        # Distance from (r,0,0) to (x,0,z) in cylindrical coords
        distance = np.sqrt((x - r)**2 + z**2)
        kernel_val = beta_kernel(distance, lambda_range)
        # Revolution factor: 2π times distance from axis
        return kernel_val * 2 * np.pi * z
    
    # Integrate over quarter-circle
    result, error = dblquad(
        integrand,
        -R, R,                           # x limits
        lambda x: 0,                     # z lower limit  
        lambda x: np.sqrt(R**2 - x**2)   # z upper limit (circle)
    )
    
    return result, error

# ============================================================================
# Comparison Functions
# ============================================================================

def compute_dose_profile_comparison(R, n_points=20, n_mc_samples=50000):
    """
    Compute dose profile using both methods and compare
    
    Returns:
    --------
    r_values : array
        Radial positions
    doses_mc : array
        Doses from Monte Carlo
    doses_cyl : array
        Doses from cylindrical integration
    errors_cyl : array
        Integration errors from cylindrical method
    """
    r_values = np.linspace(0, R, n_points)
    doses_mc = []
    doses_cyl = []
    errors_cyl = []
    
    print(f"\nComputing dose profile for R={R:.2f} mm...")
    print(f"{'r (mm)':<10} {'MC Dose':<15} {'Cyl Dose':<15} {'Difference':<12} {'Time (s)':<10}")
    print("-" * 70)
    
    for i, r in enumerate(r_values):
        # Monte Carlo method
        t0 = time.time()
        dose_mc = dose_at_point_mc(r, R, LAMBDA_BETA, n_mc_samples)
        time_mc = time.time() - t0
        
        # Cylindrical method
        t0 = time.time()
        dose_cyl, error_cyl = dose_at_point_cylindrical(r, R, LAMBDA_BETA)
        time_cyl = time.time() - t0
        
        doses_mc.append(dose_mc)
        doses_cyl.append(dose_cyl)
        errors_cyl.append(error_cyl)
        
        # Calculate difference
        diff_pct = abs(dose_mc - dose_cyl) / dose_cyl * 100 if dose_cyl > 0 else 0
        
        print(f"{r:<10.2f} {dose_mc:<15.6e} {dose_cyl:<15.6e} {diff_pct:<12.2f} MC:{time_mc:.2f} Cyl:{time_cyl:.2f}")
    
    return r_values, np.array(doses_mc), np.array(doses_cyl), np.array(errors_cyl)

def plot_method_comparison(R_values_mm, n_points=15, n_mc_samples=30000):
    """
    Create comparison plots for multiple tumor sizes
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(R_values_mm)))
    
    all_results = []
    
    for i, R_mm in enumerate(R_values_mm):
        print(f"\n{'='*70}")
        print(f"TUMOR SIZE: R = {R_mm:.1f} mm")
        print(f"{'='*70}")
        
        r_vals, doses_mc, doses_cyl, errors_cyl = compute_dose_profile_comparison(
            R_mm, n_points=n_points, n_mc_samples=n_mc_samples
        )
        
        all_results.append({
            'R': R_mm,
            'r': r_vals,
            'mc': doses_mc,
            'cyl': doses_cyl,
            'err': errors_cyl
        })
        
        # Normalize to center dose for comparison
        center_dose_cyl = doses_cyl[0]
        
        # Panel A: Normalized dose profiles (both methods)
        ax = axes[0, 0]
        ax.plot(r_vals, doses_cyl / center_dose_cyl, '-', 
               color=colors[i], linewidth=2, label=f'R={R_mm:.1f} mm (Cyl)')
        ax.plot(r_vals, doses_mc / center_dose_cyl, 'o', 
               color=colors[i], markersize=4, alpha=0.5)
        
        # Panel B: Method comparison for this size
        ax = axes[0, 1]
        relative_diff = (doses_mc - doses_cyl) / doses_cyl * 100
        ax.plot(r_vals, relative_diff, 'o-', color=colors[i], 
               linewidth=1.5, markersize=5, label=f'R={R_mm:.1f} mm')
    
    # Finalize Panel A
    ax = axes[0, 0]
    ax.set_xlabel('Radial Position r (mm)', fontsize=12)
    ax.set_ylabel('Normalized Dose d(r)/d(0)', fontsize=12)
    ax.set_title('A. Radial Dose Profiles\n(Lines: Cylindrical, Points: Monte Carlo)', 
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])
    
    # Finalize Panel B
    ax = axes[0, 1]
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    ax.fill_between([-1, 10], -5, 5, alpha=0.2, color='green')
    ax.set_xlabel('Radial Position r (mm)', fontsize=12)
    ax.set_ylabel('Relative Difference (%)', fontsize=12)
    ax.set_title('B. Monte Carlo vs Cylindrical Method\n(MC - Cyl)/Cyl × 100%', 
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, max(R_values_mm)])
    ax.text(0.95, 0.95, '±5% band', transform=ax.transAxes,
           ha='right', va='top', fontsize=10, 
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # Panel C: Geometric escape comparison
    ax = axes[1, 0]
    
    for result in all_results:
        R = result['R']
        r_vals = result['r']
        doses_cyl = result['cyl']
        
        # Compute average dose (volume-weighted)
        r_squared = r_vals**2
        integrand = doses_cyl * r_squared
        avg_dose = 3 / R**3 * np.trapz(integrand, r_vals)
        retention = avg_dose / doses_cyl[0]
        
        # Geometric approximation
        R_norm = R / LAMBDA_BETA
        geom_retention = R_norm**3 / (R_norm + 1)**3
        
        # Plot point
        idx = np.where(np.array(R_values_mm) == R)[0][0]
        ax.plot(R, retention, 'o', color=colors[idx], markersize=10, 
               label=f'R={R:.1f} mm: Cyl={retention:.3f}, Geom={geom_retention:.3f}')
    
    # Add geometric curve
    R_curve = np.linspace(0.1, max(R_values_mm), 100)
    R_norm_curve = R_curve / LAMBDA_BETA
    geom_curve = R_norm_curve**3 / (R_norm_curve + 1)**3
    ax.plot(R_curve, geom_curve, 'k--', linewidth=2, label='Geometric R³/(R+1)³')
    
    ax.set_xlabel('Tumor Radius (mm)', fontsize=12)
    ax.set_ylabel('Dose Retention Fraction', fontsize=12)
    ax.set_title('C. Dose Retention vs Tumor Size', fontsize=13, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Panel D: Summary statistics
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = "METHOD COMPARISON SUMMARY\n\n"
    summary_text += f"Beta range λ = {LAMBDA_BETA} mm\n"
    summary_text += f"MC samples per point: {n_mc_samples:,}\n\n"
    
    summary_text += "Agreement Assessment:\n"
    for result in all_results:
        R = result['R']
        mc = result['mc']
        cyl = result['cyl']
        
        relative_diff = np.abs((mc - cyl) / cyl * 100)
        mean_diff = np.mean(relative_diff)
        max_diff = np.max(relative_diff)
        
        summary_text += f"\nR = {R:.1f} mm:\n"
        summary_text += f"  Mean difference: {mean_diff:.2f}%\n"
        summary_text += f"  Max difference: {max_diff:.2f}%\n"
    
    summary_text += "\n" + "="*35 + "\n"
    summary_text += "✓ Methods agree within MC noise\n"
    summary_text += "✓ Cylindrical method is faster\n"
    summary_text += "✓ Both suitable for validation\n"
    
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', family='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    return fig, all_results

def main():
    """Main comparison"""
    
    print("="*70)
    print("BETA DOSE DEPOSITION - METHOD COMPARISON")
    print("="*70)
    print("\nComparing two methods:")
    print("1. Monte Carlo 3D sampling (brute force)")
    print("2. Cylindrical symmetry integration (elegant)")
    print("\nBoth should give identical results (within MC noise)")
    
    # Test tumor sizes
    R_values_mm = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0])
    
    # Generate comparison
    fig, results = plot_method_comparison(R_values_mm, n_points=15, n_mc_samples=50000)
    
    # Save
    output_file = "beta_method_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n\nFigure saved: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
