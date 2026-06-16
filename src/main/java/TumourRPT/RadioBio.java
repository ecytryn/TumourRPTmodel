package TumourRPT;

import HAL.Tools.ODESolver.ODESolver;

import java.util.ArrayList;

/**
 * RadioBio: Radiation biology calculations using ODE-based Lea-Catcheside model
 * 
 * APPROACH:
 * Pre-allocate ODE state for every hour of the simulation. Each hour represents
 * a birth cohort - all cells born in that hour share the same radiation exposure history.
 * 
 * MATHEMATICAL BACKGROUND:
 * 
 * The Lea-Catcheside G-factor accounts for sub-lethal damage repair:
 *   G = (2/D²) * ∫₀ᵀ Ḋ(t) A(t) dt
 * 
 * where A(t) is the accumulated effective damage:
 *   A(t) = ∫₀ᵗ Ḋ(t') exp(-μ(t-t')) dt'
 * 
 * The ODE formulation:
 *   dA/dt = D_dot(t) - mu*A(t)
 *   dG_num/dt = D_dot(t) * A(t)
 * 
 * Survival fraction:
 *   SF = exp(-α*D - β*G*D²)
 *
 * The 2/D^2 in the definition of G causes unnecessary trouble for small D so I'm reformulating the SF: 
 *   G = 2 G_num and
 *   SF = exp(-α*D - β*G) 
 * This is totally equivalent - I've just cancelled a D^2 on top and bottom.
 * 
 * IMPLEMENTATION:
 * - hourStates[h] stores {D, A, G_num} for cells born at hour h
 * - Every simulation hour: update all existing hour states with current dose rate
 * - Cells query their SF based on birthTime (hour index)
 */
public class RadioBio {

    public Grid grid;
    
    /**
     * ODE state for one hour cohort (all cells born in same hour)
     */
    private class HourState {
        double D = 0.0;       // Total dose (Gy)
        double A = 0.0;       // Amplitude (accumulated effective damage)
        double G_num = 0.0;   // G-factor numerator
    }
    
    private HourState[] hourStates;  // Index = birth hour
    private double mu;               // Repair rate (1/hour)
    private int maxSimulationHours;

    public RadioBio(Grid grid) {
        this.grid = grid;
        this.mu = SimParams.repairRate;  // 1/hour
        
        // Pre-allocate for maximum simulation length
        this.maxSimulationHours = 275 * 24;  // 250 days worth
        this.hourStates = new HourState[maxSimulationHours];
        
        // Initialize all states to zero
        for (int i = 0; i < maxSimulationHours; i++) {
            hourStates[i] = new HourState();
        }
    }

    /**
     * Step ODE state forward by one hour for given dose rate.
     * Uses exact analytical solution for A(t).
     * 
     * @param state Hour state to update
     * @param doseRate Current dose rate (Gy/h)
     * @param dt Time step (hours)
     */
	private void stepODE(HourState state, double doseRate, double dt) {
		// D is exact (linear accumulation)
		state.D += doseRate * dt;
		
		// G_num and A need careful handling
		// Subdivide for accuracy when A is changing rapidly
		int nSubsteps = 10;
		double dt_sub = dt / nSubsteps;
		
		for (int i = 0; i < nSubsteps; i++) {
			double A_old = state.A;
			
			// G_num accumulates doseRate * A(t)
			state.G_num += doseRate * A_old * dt_sub;
			
			// Update A using analytical solution for this substep
			double decay = Math.exp(-mu * dt_sub);
			if (mu > 1e-10) {
				state.A = decay * A_old + doseRate * (1.0 - decay) / mu;
			} else {
				state.A = decay * A_old + doseRate * dt_sub;
			}
		}
	}
/**
	// This version of stepODE uses dt= 1 hr which is too big during the peak after an injeciton
    private void stepODE(HourState state, double doseRate, double dt) {
        // IMPORTANT: Use A_old for G_num BEFORE updating A
        double A_old = state.A;
        state.G_num += doseRate * A_old * dt;
        
        // Update total dose
        state.D += doseRate * dt;
        
        // Update A using exact analytical solution
        double decay = Math.exp(-mu * dt);
        if (mu > 1e-10) {
            state.A = decay * A_old + doseRate * (1.0 - decay) / mu;
        } else {
            // For very small mu, use Euler approximation
            state.A = decay * A_old + doseRate * dt;
        }
    }
*/

    /**
     * Update all hour states with the most recent dose rate.
     * Called once per hour at the start of each simulation hour.
     */
    public void updateAllStates() {
        int globalTime = SimParams.globalTime;
        if (globalTime == 0) return;
        
        // Get dose rate for the hour that just completed
        double doseRate = getDoseRate(globalTime - 1);
        double dt = 1.0;  // hours
        
        // Update all hours that represent existing cell cohorts
        // (All hours from 0 to globalTime-1)
        for (int hour = 0; hour < globalTime; hour++) {
            stepODE(hourStates[hour], doseRate, dt);
        }
    }

    /**
     * Get dose rate for a given hour from DoseRateList.
     * 
     * DoseRateList indexing:
     * - Index 0: Init placeholder (hour -1)
     * - Index 1: Hour 0
     * - Index h+1: Hour h
     * 
     * @param hour The hour to query (0-based)
     * @return Dose rate in Gy/h, or 0.0 if not available
     */
    private double getDoseRate(int hour) {
        int index = hour + 1;
        
        if (index >= 0 && index < grid.DoseRateList.size()) {
            double[] data = grid.DoseRateList.get(index);
            return (data.length > 0) ? data[0] : 0.0;
        }
        return 0.0;
    }

    /**
     * Get the ODE state (D, A, G_num) for a given birth cohort.
     * Useful for debugging: lets you see whether dose is accumulating
     * and whether A and G_num are responding.
     * 
     * @param birthTime Hour when cell was born
     * @return double[3] = {D, A, G_num}, or {0,0,0} if birthTime out of bounds
     */
    public double[] getCohortState(int birthTime) {
        if (birthTime < 0 || birthTime >= maxSimulationHours) {
            System.err.printf("ERROR: birthTime %d out of bounds [0, %d) - returning zeros for D, A, G_num\n", 
                             birthTime, maxSimulationHours);
            return new double[] {0.0, 0.0, 0.0};
        }
        HourState state = hourStates[birthTime];
        return new double[] { state.D, state.A, state.G_num };
    }

    /**
     * Calculate survival fraction for a cell based on its birth time.
     * 
     * @param birthTime Hour when cell was born
     * @param cellType Cell type (NORMAL or HYPOXIC)
     * @return Survival probability [0, 1]
     */
    public double calculateSF(int birthTime, int cellType) {
        // Get the ODE state for this birth cohort
        if (birthTime < 0 || birthTime >= maxSimulationHours) {
            System.err.printf("ERROR: birthTime %d out of bounds [0, %d)\n", 
                             birthTime, maxSimulationHours);
            return 1.0;
        }
        
        HourState state = hourStates[birthTime];
        
        // Compute G-factor
        double D = state.D;
        double G_num = state.G_num;
        double G = 2.0 * G_num;
//		This commented section is the original G*D^2 formulation
//        double G = 0.0;
//        if (D > 1e-10) {
//            G = 2.0 * G_num / (D * D);
//            G = Math.max(0.0, Math.min(1.0, G));
//        }
        
        // Get radiosensitivity parameters for cell type
        double alpha, beta;
        if (cellType == SimParams.HYPOXIC) {
            alpha = SimParams.ALPHA_HYPOXIC;
            beta = SimParams.BETA_HYPOXIC;
        } else {
            alpha = SimParams.ALPHA_NORMAL;
            beta = SimParams.BETA_NORMAL;
        }
        
        // Calculate survival fraction using LQ model with G-factor
        double SF = Math.exp(-alpha * D - beta * G);
//		This commented section is the original G*D^2 formulation
//        double SF = Math.exp(-alpha * D - beta * G * D * D);
        
        return SF;
    }

    /**
     * Print diagnostic information about radiation states.
     * Useful for debugging dose accumulation issues.
     */
    public void printDiagnostics(int currentDay, int currentHour) {
        int globalTime = currentDay * 24 + currentHour;
        
        System.out.printf("\n=== RADIOBIO DIAGNOSTICS Day %d Hour %d ===\n", currentDay, currentHour);
        System.out.printf("Global time: %d hours\n", globalTime);
        
        // Sample a few birth cohorts
        int[] sampleHours = {0, globalTime/2, Math.max(0, globalTime-24)};
        
        for (int hour : sampleHours) {
            if (hour < globalTime) {
                HourState state = hourStates[hour];
                double age = (globalTime - hour) / 24.0;
                
                double G = 2.0 * state.G_num;
//		This commented section is the original G*D^2 formulation
//                double G = (state.D > 1e-10) ? 2.0 * state.G_num / (state.D * state.D) : 0.0;
                double SF_norm = calculateSF(hour, SimParams.NORMAL);
                double SF_hypo = calculateSF(hour, SimParams.HYPOXIC);
                
                System.out.printf("  Hour %d (age %.1f days): D=%.4e Gy, G=%.4f, SF_norm=%.4f, SF_hypo=%.4f\n",
                                 hour, age, state.D, G, SF_norm, SF_hypo);
            }
        }
        
        System.out.println("===========================================\n");
    }

    // Legacy method signatures for compatibility
    public void RHS(double t, double[] currentValues, double[] derivatives) {}
    public void UpdateParams() {}
    public void Solve(double[] curr, double[] next, double t0, double tf, double dt, double tol) {}
}
