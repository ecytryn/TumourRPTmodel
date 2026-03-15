#!/usr/bin/env python3
"""
Generate lookup table for beta retention fraction

Outputs Java array initialization code that can be copy-pasted directly
into the BetaRetention.java class.
"""

import numpy as np
from scipy.integrate import quad

ELL = 0.41  # mm, mean distance travelled, this is value that gives the best fit to TJ's retention table.

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

def compute_retention_fraction(R, ell=ELL, n_points=100):
    """Compute retention fraction using uniform deposition sphere model"""
    V_ell = (4.0/3.0) * np.pi * ell**3
    
    def integrand(r):
        V_int = sphere_intersection_volume(r, R, ell)
        return (V_int / V_ell) * r**2
    
    integral, error = quad(integrand, 0, R, limit=n_points)
    f = (3.0 / R**3) * integral
    
    return f

def generate_lookup_table():
    """Generate lookup table with good coverage"""
    
    # Fine spacing for small tumors (where retention changes rapidly)
    R_fine = np.arange(0.1, 1.0, 0.1)
    
    # Medium spacing for medium tumors
    R_medium = np.arange(1.0, 2.0, 0.2)
    
    # Coarser spacing for large tumors (changes slowly)
    R_coarse = np.arange(2.0, 5.5, 0.5)
    
    R_values = np.concatenate([R_fine, R_medium, R_coarse])
    
    print(f"Computing retention fractions for {len(R_values)} tumor radii...")
    print(f"Range: {R_values[0]:.2f} to {R_values[-1]:.2f} mm")
    print()
    
    retentions = []
    for R in R_values:
        f = compute_retention_fraction(R, ELL)
        retentions.append(f)
        print(f"R = {R:5.2f} mm  →  f = {f:.6f}")
    
    return R_values, retentions

def format_java_array(values, name, per_line=6):
    """Format array values for Java code"""
    lines = []
    lines.append(f"    private static final double[] {name} = {{")
    
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        if i + per_line >= len(values):
            # Last line, no comma after
            line = "        " + ", ".join(f"{v:.4f}" for v in chunk)
        else:
            line = "        " + ", ".join(f"{v:.4f}" for v in chunk) + ","
        lines.append(line)
    
    lines.append("    };")
    return "\n".join(lines)

def main():
    print("="*70)
    print("BETA RETENTION LOOKUP TABLE GENERATOR")
    print("="*70)
    print(f"Mean free path: ℓ = {ELL} mm")
    print()
    
    # Generate lookup table
    R_values, retentions = generate_lookup_table()
    
    print()
    print("="*70)
    print("JAVA CODE (Copy this into BetaRetention.java)")
    print("="*70)
    print()
    
    # Format for Java
    print(format_java_array(R_values, "R_VALUES_MM"))
    print()
    print(format_java_array(retentions, "RETENTION_FRACTIONS"))
    print()
    
    # Also save to file
    with open("beta_retention_lookup.txt", "w") as f:
        f.write("// Beta Retention Lookup Table\n")
        f.write(f"// Mean free path: {ELL} mm\n")
        f.write(f"// Generated with {len(R_values)} points from {R_values[0]:.2f} to {R_values[-1]:.2f} mm\n\n")
        f.write(format_java_array(R_values, "R_VALUES_MM"))
        f.write("\n\n")
        f.write(format_java_array(retentions, "RETENTION_FRACTIONS"))
        f.write("\n")
    
    print("Lookup table also saved to: beta_retention_lookup.txt")
    print()
    
    # Validation
    print("="*70)
    print("VALIDATION")
    print("="*70)
    print("Testing interpolation accuracy at intermediate points...")
    
    # Test at points between lookup values
    test_R = [0.15, 0.45, 0.85, 1.5, 2.75]
    for R_test in test_R:
        f_exact = compute_retention_fraction(R_test, ELL)
        
        # Simple linear interpolation (what Java will do)
        idx = np.searchsorted(R_values, R_test)
        if idx == 0:
            f_interp = retentions[0]
        elif idx >= len(R_values):
            f_interp = retentions[-1]
        else:
            t = (R_test - R_values[idx-1]) / (R_values[idx] - R_values[idx-1])
            f_interp = retentions[idx-1] + t * (retentions[idx] - retentions[idx-1])
        
        error = abs(f_exact - f_interp) / f_exact * 100
        print(f"R = {R_test:5.2f} mm:  Exact = {f_exact:.6f},  Interp = {f_interp:.6f},  Error = {error:.2f}%")
    
    print()
    print("Maximum interpolation error should be < 1%")

if __name__ == "__main__":
    main()
