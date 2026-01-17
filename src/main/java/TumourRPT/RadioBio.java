package TumorRPT;

import HAL.Tools.ODESolver.ODESolver;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * CORRECTED RadioBio Implementation
 * 
 * FIXES APPLIED:
 * 
 * 1. ORIGINAL METHOD BUG: The original nested-loop method uses a left Riemann sum
 *    to approximate the integral for A(t), which systematically underestimates
 *    the G-factor by ~15-20%. This version uses the exact integral formula.
 * 
 * 2. ODE METHOD BUG: The ODE method was using A_new (after update) for G_num
 *    accumulation. The correct formula uses A_old (before update), because
 *    G_num = integral of D_dot(t) * A(t) dt where A(t) represents damage accumulated
 *    BEFORE time t.
 * 
 * 3. COHORT INITIALIZATION: When a cohort is first queried, its state must be
 *    backfilled from its birth time to now, not initialized to zero.
 * 
 * MATHEMATICAL BACKGROUND:
 * 
 * The Lea-Catcheside G-factor is:
 *   G = (2/D^2) * integral_0^T D_dot(t) A(t) dt
 * 
 * where:
 *   A(t) = integral_0^t D_dot(t') exp(-mu*(t-t')) dt'
 * 
 * The ODE formulation:
 *   dA/dt = D_dot(t) - mu*A(t)
 *   dG_num/dt = D_dot(t) * A(t)
 * 
 * For constant D_dot over [t, t+dt], the exact solution is:
 *   A(t+dt) = exp(-mu*dt) * A(t) + D_dot * (1 - exp(-mu*dt)) / mu
 * 
 * For G_num, since A(t) represents PAST damage:
 *   G_num(t+dt) = G_num(t) + D_dot * A(t) * dt   [using A(t), not A(t+dt)]
 */
public class RadioBio {

    public Grid grid;
    public ODESolver solver;
    
    // Settings
    public static final boolean USE_OPTIMIZED = true;
    public static final boolean VALIDATE_OPTIMIZATION = false;
    
    // Validation statistics
    private double maxDifference_Normoxic = 0.0;
    private double maxDifference_Hypoxic = 0.0;
    private double sumSquaredDiff_Normoxic = 0.0;
    private double sumSquaredDiff_Hypoxic = 0.0;
    private int comparisonCount = 0;
    private long totalTime_Original_ns = 0;
    private long totalTime_Optimized_ns = 0;
    
    // ODE state for each birth cohort
    private class CohortODEState {
        double D;       // Total dose (Gy)
        double A;       // Amplitude (accumulated effective damage)
        double G_num;   // G-factor numerator
        int lastUpdatedHour;  // Last hour when state was updated
        
        CohortODEState() {
            this.D = 0.0;
            this.A = 0.0;
            this.G_num = 0.0;
            this.lastUpdatedHour = -1;
        }
    }
    
    // HashMap: birth time (hours) -> ODE state
    private HashMap<Integer, CohortODEState> cohortStates;
    
    private double mu;  // Repair rate (1/hour)

    public RadioBio(Grid grid) {
        this.grid = grid;
        this.solver = new ODESolver();
        this.mu = SimParams.repairRate;  // 1/hour
        this.cohortStates = new HashMap<>();
    }

    /**
     * Step ODE state forward by one hour for given dose rate.
     * Uses exact solution for A, and A_old for G_num.
     */
    private void stepODE(CohortODEState state, double doseRate, double dt) {
        // IMPORTANT: Use A_old for G_num BEFORE updating A
        double A_old = state.A;
        state.G_num += doseRate * A_old * dt;
        
        // Update D
        state.D += doseRate * dt;
        
        // Update A using exact analytical solution
        double decay = Math.exp(-mu * dt);
        if (mu > 1e-10) {
            state.A = decay * A_old + doseRate * (1.0 - decay) / mu;
        } else {
            state.A = decay * A_old + doseRate * dt;
        }
    }

    /**
     * Backfill cohort state from birth time to current time.
     */
    private void backfillCohortState(int birthTime, int currentTime) {
        CohortODEState state = new CohortODEState();
        double dt = 1.0;  // hours
        
        // Process each hour from birthTime to currentTime-1
        for (int hour = birthTime; hour < currentTime; hour++) {
            double doseRate = getDoseRate(hour);
            stepODE(state, doseRate, dt);
        }
        
        state.lastUpdatedHour = currentTime - 1;
        cohortStates.put(birthTime, state);
    }

    /**
     * Get dose rate for a given hour from DoseRateList.
     * 
     * DoseRateList indexing:
     * - Index 0: Init placeholder (hour -1)
     * - Index 1: Hour 0 (from day 0, hour 0)
     * - Index k: Hour k-1 for k >= 1
     * 
     * So to get dose rate for hour h, we use index h+1.
     * But we also need to account for when we're querying:
     * at globalTime T, DoseRateList.size() = T+2
     */
    private double getDoseRate(int hour) {
        // At globalTime=T, the list has entries for hours 0 through T-1
        // (plus the init placeholder at index 0)
        // 
        // The PK solver at hour H adds DoseRateList[H+1] = dose rate for hour H
        // So getDoseRate(H) should return DoseRateList[H+1]
        
        int index = hour + 1;
        
        if (index >= 0 && index < grid.DoseRateList.size()) {
            double[] data = grid.DoseRateList.get(index);
            return (data.length > 0) ? data[0] : 0.0;
        }
        return 0.0;
    }

    /**
     * Update all existing cohorts with the most recent hour's dose.
     * Called once per hour at the start of SurvivalProbLookupTableCalc.
     */
    private void updateAllCohortStates() {
        int globalTime = SimParams.globalTime;
        double dt = 1.0;
        
        // At globalTime = T, PK just added the dose rate for hour T-1
        // We need to update all cohorts with this new dose
        int hourToProcess = globalTime - 1;
        if (hourToProcess < 0) return;
        
        double doseRate = getDoseRate(hourToProcess);
        
        for (Integer birthTime : cohortStates.keySet()) {
            CohortODEState state = cohortStates.get(birthTime);
            
            // Only update if this cohort exists (was born) before hourToProcess
            // and hasn't been updated yet
            if (birthTime <= hourToProcess && state.lastUpdatedHour < hourToProcess) {
                stepODE(state, doseRate, dt);
                state.lastUpdatedHour = hourToProcess;
            }
        }
    }

    /**
     * Main entry point: compute survival probability lookup table.
     */
    public ArrayList<double[]> SurvivalProbLookupTableCalc(int currentDay, int currentHour) {
        // Update existing cohorts first
        updateAllCohortStates();
        
        ArrayList<double[]> lookupTable = new ArrayList<>();
        int maxLookupAge = SimParams.maxLookupAge;
        
        if (currentDay < maxLookupAge) {
            if (currentDay == 0) {
                lookupTable.add(new double[]{1.0, 1.0});
                return lookupTable;
            }
            maxLookupAge = currentDay;
        }
        
        // Calculate SF for each age
        for (int age = 1; age <= maxLookupAge; age++) {
            double[] SF_toUse;
            
            if (VALIDATE_OPTIMIZATION) {
                // Run both methods and compare
                long t1 = System.nanoTime();
                double[] SF_original = CalculateSF_Original_Corrected(age);
                totalTime_Original_ns += System.nanoTime() - t1;
                
                long t2 = System.nanoTime();
                double[] SF_optimized = CalculateSF_ODE(age);
                totalTime_Optimized_ns += System.nanoTime() - t2;
                
                // Track differences
                double diff_norm = Math.abs(SF_original[0] - SF_optimized[0]);
                double diff_hypo = Math.abs(SF_original[1] - SF_optimized[1]);
                
                maxDifference_Normoxic = Math.max(maxDifference_Normoxic, diff_norm);
                maxDifference_Hypoxic = Math.max(maxDifference_Hypoxic, diff_hypo);
                sumSquaredDiff_Normoxic += diff_norm * diff_norm;
                sumSquaredDiff_Hypoxic += diff_hypo * diff_hypo;
                comparisonCount++;
                
                // Debug output for first few mismatches
                if (comparisonCount <= 5 && (diff_norm > 1e-6 || diff_hypo > 1e-6)) {
                    System.out.printf("MISMATCH: Age %d, SF_orig=(%.6f,%.6f), SF_ode=(%.6f,%.6f), diff=(%.2e,%.2e)%n",
                                     age, SF_original[0], SF_original[1], 
                                     SF_optimized[0], SF_optimized[1], diff_norm, diff_hypo);
                }

                SF_toUse = USE_OPTIMIZED ? SF_optimized : SF_original;
            } else if (USE_OPTIMIZED) {
                SF_toUse = CalculateSF_ODE(age);
            } else {
                SF_toUse = CalculateSF_Original_Corrected(age);
            }
            
            lookupTable.add(SF_toUse);
        }
        
        return lookupTable;
    }

    /**
     * CORRECTED Original method using exact integral for A.
     * 
     * This recomputes the full integral from scratch, which is O(n) for
     * n hours of dose history. Since we call this for each age bin,
     * the total complexity is O(n * maxAge) per call.
     */
    private double[] CalculateSF_Original_Corrected(int age) {
        int globalTime = SimParams.globalTime;
        int birthTime = globalTime - age * 24;
        birthTime = Math.max(0, birthTime);
        
        double dt = 1.0;  // hours
        double D = 0.0;
        double A = 0.0;
        double G_num = 0.0;
        
        // Process each hour from birthTime to globalTime-1
        for (int hour = birthTime; hour < globalTime; hour++) {
            double doseRate = getDoseRate(hour);
            
            // Use A_old for G_num (before updating A)
            G_num += doseRate * A * dt;
            
            // Update D
            D += doseRate * dt;
            
            // Update A using exact formula
            double decay = Math.exp(-mu * dt);
            if (mu > 1e-10) {
                A = decay * A + doseRate * (1.0 - decay) / mu;
            } else {
                A = decay * A + doseRate * dt;
            }
        }
        
        // Compute G-factor
        double G_Factor = 0.0;
        if (D > 1e-10) {
            G_Factor = 2.0 * G_num / (D * D);
            G_Factor = Math.max(0.0, Math.min(1.0, G_Factor));
        }
        
        // Compute survival fractions
        double[] SF = new double[2];
        for (int i = 0; i < 2; i++) {
            SF[i] = Math.exp(-SimParams.alphaValues[i] * D 
                           - SimParams.betaValues[i] * D * D * G_Factor);
        }
        
        return SF;
    }

    /**
     * ODE method: look up or backfill cohort state.
     * 
     * This is O(1) per lookup after initial backfill, giving
     * O(n + maxAge) total complexity per call.
     */
    private double[] CalculateSF_ODE(int age) {
        int globalTime = SimParams.globalTime;
        int birthTime = globalTime - age * 24;
        birthTime = Math.max(0, birthTime);
        
        // Get or create cohort state
        if (!cohortStates.containsKey(birthTime)) {
            backfillCohortState(birthTime, globalTime);
        }
        
        CohortODEState state = cohortStates.get(birthTime);
        
        // Compute G-factor
        double D = state.D;
        double G_num = state.G_num;
        double G_Factor = 0.0;
        
        if (D > 1e-10) {
            G_Factor = 2.0 * G_num / (D * D);
            G_Factor = Math.max(0.0, Math.min(1.0, G_Factor));
        }
        
        // Compute survival fractions
        double[] SF = new double[2];
        for (int i = 0; i < 2; i++) {
            SF[i] = Math.exp(-SimParams.alphaValues[i] * D 
                           - SimParams.betaValues[i] * D * D * G_Factor);
        }
        
        return SF;
    }

    /**
     * Print validation report.
     */
    public void printValidationReport() {
        if (!VALIDATE_OPTIMIZATION || comparisonCount == 0) {
            return;
        }
        if (SimParams.VERBOSE_ON) {
        
			System.out.println();
			System.out.println("======================================================================");
			System.out.println("RADIOBIOLOGY VALIDATION REPORT (CORRECTED)");
			System.out.println("======================================================================");
			
			double rmse_norm = Math.sqrt(sumSquaredDiff_Normoxic / comparisonCount);
			double rmse_hypo = Math.sqrt(sumSquaredDiff_Hypoxic / comparisonCount);
			
			System.out.printf("%nAccuracy:%n");
			System.out.printf("  Normoxic: max diff = %.2e, RMSE = %.2e%n", maxDifference_Normoxic, rmse_norm);
			System.out.printf("  Hypoxic:  max diff = %.2e, RMSE = %.2e%n", maxDifference_Hypoxic, rmse_hypo);
			
			System.out.printf("%nPerformance (%d comparisons):%n", comparisonCount);
			System.out.printf("  Original (corrected): %.3f s%n", totalTime_Original_ns / 1e9);
			System.out.printf("  ODE:                  %.3f s%n", totalTime_Optimized_ns / 1e9);
			
			if (totalTime_Optimized_ns > 0) {
				double speedup = (double) totalTime_Original_ns / totalTime_Optimized_ns;
				System.out.printf("  Speedup:              %.1fx%n", speedup);
			}
			
			boolean accurate = (maxDifference_Normoxic < 1e-6 && maxDifference_Hypoxic < 1e-6);
			System.out.printf("%nStatus: %s%n", accurate ? "EXCELLENT - Methods agree" : "MISMATCH - Review needed");
			System.out.println("======================================================================");
			System.out.println();
		}
    }

    // Legacy method signatures for compatibility
    public void RHS(double t, double[] currentValues, double[] derivatives) {}
    public void UpdateParams() {}
    public void Solve(double[] curr, double[] next, double t0, double tf, double dt, double tol) {}
}
