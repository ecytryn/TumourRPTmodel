package TumorRPT;

/**
 * Intratumoural energy deposition fraction calculator
 * 
 * Uses uniform deposition sphere model with pre-computed lookup table
 * for computational efficiency. In principle, this could be replaced by
 * an analytical form (the integrals end up being piecewise polynommial).
 * 
 * Physical model:
 * - Each Lu-177 decay deposits energy uniformly in sphere of radius ℓ (mean range)
 * - Computes fraction of energy retained within tumor using geometric overlap calculation
 * - See appendix for full mathematical derivation
 * 
 * @author Generated for TumorRPT model
 */
public class BetaRetention {
    
    /** Mean free path for Lu-177 beta particles in tissue (mm) */
    private static final double ELL_MM = 1.0;
    
    /** 
     * Lookup table: tumor radii in mm
     * Range: 0.1 to 5.0 mm with finer resolution for small tumors
     */
    private static final double[] R_VALUES_MM = {
        0.1000, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000,
        0.7000, 0.8000, 0.9000, 1.0000, 1.2000, 1.4000,
        1.6000, 1.8000, 2.0000, 2.5000, 3.0000, 3.5000,
        4.0000, 4.5000, 5.0000
    };
    
    /**
     * Lookup table: retention fractions corresponding to R_VALUES_MM
     * Computed by
     * Taylor J. McColl, Nouran R. R. Zaid, Carlos F. Uribe, Arman Rahmim
     * A unified multi-scale S-value framework enabling pharmacokinetic-pharmacodynamic 
     * model coupling in radiopharmaceutical therapy
     * Proc. SNMMI Annual Meeting, 2026 
     */
    private static final double[] RETENTION_FRACTIONS_TJ = {
		0.1982, 0.2899, 0.3641, 0.4275, 0.4861, 0.5408, 
		0.5921, 0.642, 0.6898, 0.7342, 0.8137, 0.876, 
		0.9187, 0.9449, 0.9593, 0.9701, 0.97, 0.9695, 
		0.97, 0.9714, 0.9734
	};

    /**
     * Lookup table: retention fractions corresponding to R_VALUES_MM
     * Computed using uniform deposition sphere model with ℓ = 1.0 mm
     */
    private static final double[] RETENTION_FRACTIONS_UNIFORM_DEPOSITION = {
        0.0010, 0.0080, 0.0270, 0.0640, 0.1250, 0.2072,
        0.2875, 0.3579, 0.4179, 0.4688, 0.5493, 0.6096,
        0.6561, 0.6929, 0.7227, 0.7770, 0.8137, 0.8400,
        0.8599, 0.8753, 0.8878
    };

    /**
     * Lookup table: retention fractions corresponding to R_VALUES_MM
     * A blend of the uniform sphere at low R and TJ's table at higher R
     */
    private static final double[] RETENTION_FRACTIONS_BLEND = {
		0.011956, 0.03932, 0.083183, 0.1447778, 0.2253056, 0.3184,
		0.405956, 0.484167, 0.55385, 0.61624, 0.725567, 0.8168,
		0.889522, 0.9449, 0.9593, 0.9701, 0.97, 0.9695,
		0.97, 0.9714, 0.9734
    };

	/**
	 * Pick which table to use.
	*/
    private static final double[] RETENTION_FRACTIONS = RETENTION_FRACTIONS_BLEND;
        
    /**
     * Get beta retention fraction for given tumor radius
     * 
     * Uses linear interpolation between tabulated values for efficiency.
     * Extrapolates (with warning) for radii outside lookup table range.
     * 
     * @param R_mm Tumor radius in millimeters
     * @return Retention fraction (0 to 1)
     */
    public static double getRetentionFraction(double R_mm) {
        
        // Validate input
        if (R_mm < 0) {
            throw new IllegalArgumentException("Tumor radius cannot be negative: " + R_mm);
        }
        
        if (R_mm == 0) {
            return 0.0;  // No retention for zero-size tumor
        }
        
        // Check bounds and extrapolate if needed
        if (R_mm <= R_VALUES_MM[0]) {
            // Small tumor extrapolation: f ≈ (R/ℓ)³ for R << ℓ
            // But table starts at 0.1 mm, so just use first value
            return RETENTION_FRACTIONS[0] * Math.pow(R_mm / R_VALUES_MM[0], 3);
        }
        
        if (R_mm >= R_VALUES_MM[R_VALUES_MM.length - 1]) {
            // Large tumor extrapolation: f → 1 - 3ℓ/R for R >> ℓ
            if (R_mm > 10.0 && SimParams.VERBOSE_ON) {
                System.out.println("Warning: BetaRetention extrapolating for large R = " + R_mm + " mm");
            }
            // Use asymptotic formula
            return 1.0 - 3.0 * ELL_MM / R_mm;
        }
        
        // Linear interpolation
        for (int i = 0; i < R_VALUES_MM.length - 1; i++) {
            if (R_mm >= R_VALUES_MM[i] && R_mm <= R_VALUES_MM[i + 1]) {
                // Interpolation parameter
                double t = (R_mm - R_VALUES_MM[i]) / (R_VALUES_MM[i + 1] - R_VALUES_MM[i]);
                
                // Linear interpolation
                double f = RETENTION_FRACTIONS[i] + t * (RETENTION_FRACTIONS[i + 1] - RETENTION_FRACTIONS[i]);
                
                return f;
            }
        }
        
        // Should never reach here, but return last value as fallback
        return RETENTION_FRACTIONS[RETENTION_FRACTIONS.length - 1];
    }
    
    /**
     * Get beta escape fraction (complement of retention)
     * 
     * @param R_mm Tumor radius in millimeters
     * @return Escape fraction (0 to 1)
     */
    public static double getEscapeFraction(double R_mm) {
        return 1.0 - getRetentionFraction(R_mm);
    }
    
    /**
     * Get mean free path used in this model
     * 
     * @return Mean free path in millimeters
     */
    public static double getMeanFreePath() {
        return ELL_MM;
    }
    
    /**
     * Test method to verify lookup table interpolation
     */
    public static void main(String[] args) {
        System.out.println("Beta Retention Fraction Lookup Table Test");
        System.out.println("==========================================");
        System.out.println("Mean free path: " + ELL_MM + " mm");
        System.out.println();
        
        // Test at lookup table points
        System.out.println("Exact table values:");
        for (int i = 0; i < R_VALUES_MM.length; i++) {
            double R = R_VALUES_MM[i];
            double f = getRetentionFraction(R);
            System.out.printf("R = %5.2f mm  →  f = %.4f%n", R, f);
        }
        
        System.out.println();
        System.out.println("Interpolated values:");
        
        // Test interpolation
        double[] testR = {0.15, 0.45, 0.85, 1.5, 2.75, 3.75};
        for (double R : testR) {
            double f = getRetentionFraction(R);
            System.out.printf("R = %5.2f mm  →  f = %.4f%n", R, f);
        }
        
        System.out.println();
        System.out.println("Boundary cases:");
        
        // Test boundaries
        System.out.printf("R = %5.2f mm  →  f = %.4f (small R extrapolation)%n", 
                         0.05, getRetentionFraction(0.05));
        System.out.printf("R = %5.2f mm  →  f = %.4f (large R extrapolation)%n",
                         7.0, getRetentionFraction(7.0));
        System.out.printf("R = %5.2f mm  →  f = %.4f (very large R)%n",
                         20.0, getRetentionFraction(20.0));
    }
}
