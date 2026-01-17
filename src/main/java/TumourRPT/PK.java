package TumorRPT;

import HAL.Tools.ODESolver.ODESolver;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * Pharmacokinetic model for radioligand therapy using quasi-steady-state approximation
 * 
 * APPROACH:
 * - Fast compartments (C_v, C_int, C_b) assumed in QSS with C_BP
 * - Only slow compartments (C_BP, C_intern) integrated as ODEs
 * - Tracks "hot" (radioactive) and "cold" (decayed) ligand separately
 * 
 * UNITS IN THIS FILE (consistent with SimParams):
 * - Time: seconds (s)
 * - Concentration: mol/m³
 * - Amount: moles (mol)
 * - Volume: m³
 * - Rate constants: 1/s or m³/(mol·s)
 * 
 * NOTE: This file integrates with SimParams.java which uses SI units throughout.
 * No conversion factors should appear in the math - everything is already SI.
 */
public class PK {

    public Grid grid;
    public ODESolver solver;
    
    // Current tumor geometry (updated hourly by Grid)
    private double V_ec;           // Extracellular volume (m³)
    private double R_total;        // Total receptors (mol)
    private double V_v;            // Total vessel volume (m³)
    private double tumorVolume;    // Total tumor volume (m³)
    private int numVessels;        // Vessel count

    public PK(Grid grid){
        this.grid = grid;
        this.solver = new ODESolver();
        
        // Initialize with default values (will be updated by updateGeometry)
        this.V_ec = 0;
        this.R_total = 0;
        this.V_v = 0;
        this.tumorVolume = 0;
        this.numVessels = 0;
    }

	// Getters for reporting values
	public double getV_ec() { 
		return V_ec; 
	}
	
	public double getV_v() { 
		return V_v; 
	}
	
	public double getR_total() { 
		return R_total; 
	}
	
	public double getTumorVolume() { 
		return tumorVolume; 
	}

    /**
     * Right-hand side of ODE system (slow compartments only)
     * 
     * State vector:
     * [0] C_cen_hot    - Central compartment hot concentration (mol/m³)
     * [1] C_cen_cold   - Central compartment cold concentration (mol/m³)
     * [2] C_ic_hot     - Intracellular hot concentration (mol/m³)
     * [3] C_ic_cold    - Intracellular cold concentration (mol/m³)
     * [4] A_blob       - Test variable for decay validation (mol)
     * 
     * Fast compartments (C_v, C_ec, C_b) computed via QSS functions
     */
    public void RHS_MINIPBPK(double t, double[] currentValues, double[] derivatives){

        double C_cen_hot = currentValues[0];    // mol/m³
        double C_cen_cold = currentValues[1];   // mol/m³
        double C_ic_hot = currentValues[2];     // mol/m³
        double C_ic_cold = currentValues[3];    // mol/m³
        double A_blob = currentValues[4];       // mol (test variable)

        // Compute fast compartments in QSS
        double C_v_hot    = calc_C_v_QSS(C_cen_hot);
        double C_v_cold   = calc_C_v_QSS(C_cen_cold);
        double C_ec_hot   = calc_C_ec_QSS(C_cen_hot);
        double C_ec_cold  = calc_C_ec_QSS(C_cen_cold);
        double C_b_hot    = calc_C_b_QSS(C_cen_hot, C_cen_hot + C_cen_cold);
        double C_b_cold   = calc_C_b_QSS(C_cen_cold, C_cen_hot + C_cen_cold);

        // === ODE System (all in SI units) ===
        
        // Central compartment (hot)
        // dC_cen_hot/dt = -clearance - decay --- flux between cen and vasc is net zero in QSS so that term is omitted
        derivatives[0] = -SimParams.LAMBDA_BIO * C_cen_hot 
                        - SimParams.LAMBDA_DECAY * C_cen_hot;
        
        // Central compartment (cold) 
        // dC_cen_cold/dt = -clearance + decay(from hot)
        derivatives[1] = -SimParams.LAMBDA_BIO * C_cen_cold 
                        + SimParams.LAMBDA_DECAY * C_cen_hot;
        
        // Intracellular (hot)
        // dC_ic_hot/dt = internalization - release - decay
        derivatives[2] = SimParams.K_INT * C_b_hot 
                        - SimParams.K_REL * C_ic_hot 
                        - SimParams.LAMBDA_DECAY * C_ic_hot;
        
        // Intracellular (cold)
        // dC_ic_cold/dt = internalization - release + decay(from hot)
        derivatives[3] = SimParams.K_INT * C_b_cold 
                        - SimParams.K_REL * C_ic_cold 
                        + SimParams.LAMBDA_DECAY * C_ic_hot;
        
        // Test blob (just decays)
        derivatives[4] = -SimParams.LAMBDA_DECAY * A_blob;
    }

    /**
     * Computes bound ligand concentration via QSS approximation
     * 
     * From: k_on*C_ec*R_free = (k_off + k_int)*C_b
     * With: R_free = R_total - C_b*V_ec
     * 
     * Solution: C_b = R_total/V_ec * C_cen / (C_cen_total + beta)
     * where beta = (k_off + k_int)/k_on
     * 
     * Uses stored geometry (updated hourly by Grid.updatePKGeometry)
     * 
     * @param C_cen Single species central concentration (mol/m³)
     * @param C_cen_total Total (hot+cold) central concentration (mol/m³)
     * @return Bound concentration C_b (mol/m³)
     */
    private double calc_C_b_QSS(double C_cen, double C_cen_total) {
        // Use stored geometry (already in SI units)
        double R_T_tilde = R_total / V_ec;  // mol/m³ (receptor concentration)
        double beta = (SimParams.K_OFF + SimParams.K_INT) / SimParams.K_ON;  // mol/m³
       
        return R_T_tilde * C_cen / (C_cen_total + beta);
    }

    /**
     * Vessel concentration QSS (C_v ≈ C_cen)
     * @param C_cen Central concentration (mol/m³)
     * @return Vessel concentration (mol/m³)
     */
    private double calc_C_v_QSS(double C_cen) {
        return C_cen;  // Fast equilibration assumption
    }
    
    /**
     * Extracellular concentration QSS (C_ec ≈ C_cen)
     * @param C_cen Central concentration (mol/m³)
     * @return Extracellular concentration (mol/m³)
     */
    private double calc_C_ec_QSS(double C_cen) {
        return C_cen;  // Fast equilibration assumption
    }

    /**
     * Main PK integration step - advances ODEs and computes dose rate
     * 
     * @param t_0 Start time (s)
     * @param t_f End time (s)
     * @param currentDay Current day number
     * @param currentHour Current hour (0-23)
     * @return 0 on success
     */
    public int DoseRateCalc(double t_0, double t_f, int currentDay, int currentHour){
        
        // Initialize on first call
        if (currentDay == 0 && currentHour == 0){
            this.grid.DoseRateList.add(new double[]{0.0});
            return 0;
        }

        // Get current state from grid
        double[] currentValues = this.grid.PKStateVariables.get(this.grid.PKStateVariables.size()-1);

        // Extract ODE variables (slow compartments only)
        double[] ODE_currentValues = {
            currentValues[0],  // C_cen_hot
            currentValues[1],  // C_cen_cold
            currentValues[8],  // C_ic_hot
            currentValues[9],  // C_ic_cold
            currentValues[10]   // A_blob
//            currentValues[2],  // C_ic_hot
//            currentValues[3],  // C_ic_cold
//            currentValues[4]   // A_blob
        };

        double[] ODE_nextValues = new double[5];

        // Integrate using HAL's Runge45 solver
        // Note: t_0 and t_f are now in seconds (not hours!)
        double result = this.solver.Runge45(
            this::RHS_MINIPBPK, 
            ODE_currentValues, 
            ODE_nextValues, 
            t_0, 
            t_f, 
            SimParams.TIME_STEP * 1e-5,  // Initial dt (much smaller than time step)
            1e-5  // Error tolerance
        );

        // Extract integrated values
        double C_cen_hot = ODE_nextValues[0];
        double C_cen_cold = ODE_nextValues[1];
        double C_ic_hot = ODE_nextValues[2];
        double C_ic_cold = ODE_nextValues[3];
        double A_blob = ODE_nextValues[4];
        
        // Compute QSS values at new time
        double C_v_hot    = calc_C_v_QSS(C_cen_hot);
        double C_v_cold   = calc_C_v_QSS(C_cen_cold);
        double C_ec_hot   = calc_C_ec_QSS(C_cen_hot);
        double C_ec_cold  = calc_C_ec_QSS(C_cen_cold);
        double C_b_hot    = calc_C_b_QSS(C_cen_hot, C_cen_hot + C_cen_cold);
        double C_b_cold   = calc_C_b_QSS(C_cen_cold, C_cen_hot + C_cen_cold);

        // Optional logging
        if (SimParams.ENABLE_PBPK_LOGGING) {
            logPBPKState(currentDay, currentHour, 
                        C_cen_hot, C_cen_cold, 
                        C_v_hot, C_v_cold,
                        C_ec_hot, C_ec_cold, 
                        C_b_hot, C_b_cold,
                        C_ic_hot, C_ic_cold);
        }
                
        // Store complete state (slow + fast variables)
        double[] nextValues = {
            C_cen_hot, C_cen_cold,
            C_v_hot, C_v_cold,
            C_ec_hot, C_ec_cold,
            C_b_hot, C_b_cold,
            C_ic_hot, C_ic_cold,
            A_blob
        };
        
        this.grid.PKStateVariables.add(nextValues);
        
        // === COMPUTE DOSE RATE ===
        // Use stored geometry (already in SI units)
        double V_v_current = this.V_v;
        double V_ec_current = this.V_ec;
        
        // Total hot activity in tumor (decays/s)
        double N_hot = (
            C_v_hot * V_v_current +
            C_ec_hot * V_ec_current +
            C_b_hot * V_ec_current +
            C_ic_hot * V_ec_current
        );  // mol/m³ * m³ = mol
        
        // Convert to Bq (decays/s per particle)
        double activityBq = N_hot * SimParams.LAMBDA_DECAY * SimParams.AVOGADRO;  // mol* 1/s * (1/mol) = particles/s = Bq
        
        // Energy emission rate (J/s) - before escape
		double energyRate_emitted = activityBq * SimParams.E_BETA_LU177;  // Bq * J = J/s

        // β particle escape correction
        // Lu-177 beta range ~1 mm in tissue
        double beta_range_m = 1.0e-3;  // m
		// Assume the tumour cells are in a sphere of some effective radius
        double R_tumor_eff = Math.pow(3.0 * tumorVolume / (4.0 * Math.PI), 1.0/3.0);
        double f_deposit = Math.pow(R_tumor_eff / (R_tumor_eff + beta_range_m), 3.0);
		// Spread the energy over a sphere of radius R_tumor_eff + beta_range_m
        double energyRate = energyRate_emitted * f_deposit;
                
        // Tumor mass (assuming density = 1000 kg/m³)
        double tumorMass = tumorVolume * 1000.0;  // m³ * kg/m³ = kg
        
        // Dose rate in Gy/s
        double doseRate_Gy_per_s = energyRate / tumorMass;  // (J/s) / kg = Gy/s
        
        // Convert to Gy/h for compatibility with old outputs
        double doseRate_Gy_per_h = doseRate_Gy_per_s * 3600.0;

        // DIAGNOSTIC: Print every 10 days
        if (currentDay % 10 == 0 && currentHour == 0) {
            System.out.printf("\n=== PK DOSE DIAGNOSTIC (Day %d) ===\n", currentDay);
            System.out.printf("Tumor: V=%.3f mm³, M=%.3f mg, R=%.3f mm\n",
                             tumorVolume*1e9, tumorMass*1e6, R_tumor_eff*1e3);
            System.out.printf("RL: N=%.3e mol, A=%.3e Bq\n", N_hot, activityBq);
            System.out.printf("Energy: %.3e J/s emitted, %.3e J/s deposited (%.1f%%)\n",
                             energyRate_emitted, energyRate, f_deposit*100);
            System.out.printf("Dose rate: %.3e Gy/h\n", doseRate_Gy_per_h);
            System.out.printf("=====================================\n\n");
        }
                
        this.grid.DoseRateList.add(new double[]{doseRate_Gy_per_h});
        
        return 0;
    }

    /**
     * Update tumor geometry (called hourly by Grid)
     * 
     * @param newTumorVolume Updated tumor volume (m³)
     * @param newReceptorMoles Updated receptor moles (mol)
     * @param newNumVessels Updated vessel count
     */
    public void updateGeometry(double newTumorVolume, double newReceptorMoles, int newNumVessels, double tumorHeight_m) {    
		this.tumorVolume = newTumorVolume;
        this.V_ec = SimParams.computeInterstitialVolume(newTumorVolume);
        this.R_total = newReceptorMoles;
    	double vesselVolumeInTumor = Math.pow(SimParams.CELL_LENGTH, 2) * tumorHeight_m;
	    this.V_v = newNumVessels * vesselVolumeInTumor;
        this.numVessels = newNumVessels;
    }

    /**
     * Log PBPK state to CSV file for debugging
     */
    private void logPBPKState(int day, int hour,
                             double C_cen_hot, double C_cen_cold,
                             double C_v_hot, double C_v_cold,
                             double C_ec_hot, double C_ec_cold,
                             double C_b_hot, double C_b_cold,
                             double C_ic_hot, double C_ic_cold) {
        File logFile = new File("Results/PBPK_log.csv");

        try (FileWriter fw = new FileWriter(logFile, true);
             BufferedWriter bw = new BufferedWriter(fw);
             PrintWriter out = new PrintWriter(bw)) {
    
            // Header if file is empty
            if (logFile.length() == 0) {
                out.println("day,hour,time_days,C_cen_hot,C_cen_cold,C_v_hot,C_v_cold," +
                           "C_ec_hot,C_ec_cold,C_b_hot,C_b_cold,C_ic_hot,C_ic_cold");
            }
    
            double timeDays = day + hour / 24.0;
            
            // Convert to nmol/L for logging (easier to read)
            double conv = 1e9;  // mol/m³ → nmol/L (factor of 1e9)
            
            out.printf("%.0f,%.0f,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e%n",
                (double) day, (double) hour, timeDays,
                C_cen_hot * conv, C_cen_cold * conv, 
                C_v_hot * conv, C_v_cold * conv, 
                C_ec_hot * conv, C_ec_cold * conv, 
                C_b_hot * conv, C_b_cold * conv,
                C_ic_hot * conv, C_ic_cold * conv
            );
    
        } catch (IOException e) {
            System.err.println("Error writing PBPK log: " + e.getMessage());
        }
    }
}
