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
 * - Fast compartments (N_v, N_int, N_b) assumed in QSS with N_BP
 * - Only slow compartments (N_BP, N_intern) integrated as ODEs
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
     * [0] N_cen_hot    - Central compartment hot amount (mol)
     * [1] N_cen_cold   - Central compartment cold amount (mol)
     * [2] N_ic_hot     - Intracellular hot amount (mol)
     * [3] N_ic_cold    - Intracellular cold amount (mol)
     * [4] A_blob       - Test variable for decay validation (mol)
     * 
     * Fast compartments (N_v, N_ec, N_b) computed via QSS functions
     */
    public void RHS_MINIPBPK(double t, double[] currentValues, double[] derivatives){

        double N_cen_hot = currentValues[0];    // mol
        double N_cen_cold = currentValues[1];   // mol
        double N_ic_hot = currentValues[2];     // mol
        double N_ic_cold = currentValues[3];    // mol
        double A_blob = currentValues[4];       // mol (test variable)

        // Compute fast compartments in QSS
        double N_v_hot    = calc_N_v_QSS(N_cen_hot);
        double N_v_cold   = calc_N_v_QSS(N_cen_cold);
        double N_ec_hot   = calc_N_ec_QSS(N_cen_hot);
        double N_ec_cold  = calc_N_ec_QSS(N_cen_cold);
        double N_b_hot    = calc_N_b_QSS(N_cen_hot, N_cen_hot + N_cen_cold);
        double N_b_cold   = calc_N_b_QSS(N_cen_cold, N_cen_hot + N_cen_cold);

        // === ODE System (all in SI units) ===
        
        // Central compartment (hot)
        // dN_cen_hot/dt = -clearance - decay 
        derivatives[0] = - SimParams.LAMBDA_BIO * N_cen_hot 
        				 - SimParams.K_INT * N_b_hot // This term appears here because of the QSS
                         - SimParams.LAMBDA_DECAY * N_cen_hot;
        
        // Central compartment (cold) 
        // dN_cen_cold/dt = -clearance + decay(from hot)
        derivatives[1] = - SimParams.LAMBDA_BIO * N_cen_cold 
        				 - SimParams.K_INT * N_b_cold // This term appears here because of the QSS
                         + SimParams.LAMBDA_DECAY * N_cen_hot;
        
        // Intracellular (hot)
        // dN_ic_hot/dt = internalization - release - decay
        derivatives[2] = SimParams.K_INT * N_b_hot 
                        - SimParams.K_REL * N_ic_hot 
                        - SimParams.LAMBDA_DECAY * N_ic_hot;
        
        // Intracellular (cold)
        // dN_ic_cold/dt = internalization - release + decay(from hot)
        derivatives[3] = SimParams.K_INT * N_b_cold 
                        - SimParams.K_REL * N_ic_cold 
                        + SimParams.LAMBDA_DECAY * N_ic_hot;
        
        // Test blob (just decays)
        derivatives[4] = -SimParams.LAMBDA_DECAY * A_blob;
    }

    /**
     * Computes bound ligand via QSS approximation
     * 
     * Uses stored geometry (updated hourly by Grid.updatePKGeometry)
     * 
     * @param N_cen Single species central amount (mol)
     * @param N_cen_total Total (hot+cold) central amount (mol)
     * @return Bound amount N_b (mol)
     */
    private double calc_N_b_QSS(double N_cen, double N_cen_total) {
        // Use stored geometry (already in SI units)
        double beta = SimParams.V_CENTRAL * (SimParams.K_OFF + SimParams.K_INT) / SimParams.K_ON;  // mol
       
        return R_total * N_cen / (N_cen_total + beta);
    }

    /**
     * Vessel amount QSS (N_v ≈ N_cen*V_v/V_cen)
     * @param N_cen Central amount (mol)
     * @return Vessel amount (mol)
     */
    private double calc_N_v_QSS(double N_cen) {
        return N_cen * V_v / SimParams.V_CENTRAL;  // Fast equilibration assumption
    }
    
    /**
     * Extracellular amount QSS (N_ec ≈ N_cen)
     * @param N_cen Central amount (mol)
     * @return Extracellular amount (mol)
     */
    private double calc_N_ec_QSS(double N_cen) {
        return N_cen * V_ec / SimParams.V_CENTRAL;  // Fast equilibration assumption
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
            currentValues[0],  // N_cen_hot
            currentValues[1],  // N_cen_cold
            currentValues[8],  // N_ic_hot
            currentValues[9],  // N_ic_cold
            currentValues[10]   // A_blob
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
        double N_cen_hot = ODE_nextValues[0];
        double N_cen_cold = ODE_nextValues[1];
        double N_ic_hot = ODE_nextValues[2];
        double N_ic_cold = ODE_nextValues[3];
        double A_blob = ODE_nextValues[4];
        
        // Compute QSS values at new time
        double N_v_hot    = calc_N_v_QSS(N_cen_hot);
        double N_v_cold   = calc_N_v_QSS(N_cen_cold);
        double N_ec_hot   = calc_N_ec_QSS(N_cen_hot);
        double N_ec_cold  = calc_N_ec_QSS(N_cen_cold);
        double N_b_hot    = calc_N_b_QSS(N_cen_hot, N_cen_hot + N_cen_cold);
        double N_b_cold   = calc_N_b_QSS(N_cen_cold, N_cen_hot + N_cen_cold);

        // Optional logging
        if (SimParams.ENABLE_PBPK_LOGGING) {
            logPBPKState(currentDay, currentHour, 
                        N_cen_hot, N_cen_cold, 
                        N_v_hot, N_v_cold,
                        N_ec_hot, N_ec_cold, 
                        N_b_hot, N_b_cold,
                        N_ic_hot, N_ic_cold);
        }
                
        // Store complete state (slow + fast variables)
        double[] nextValues = {
            N_cen_hot, N_cen_cold,
            N_v_hot, N_v_cold,
            N_ec_hot, N_ec_cold,
            N_b_hot, N_b_cold,
            N_ic_hot, N_ic_cold,
            A_blob
        };
        
        this.grid.PKStateVariables.add(nextValues);
        
        // === COMPUTE DOSE RATE ===
        // Use stored geometry (already in SI units)
        double V_v_current = this.V_v;
        double V_ec_current = this.V_ec;
        
        // Total hot activity in tumor (decays/s)
        double N_hot = (
            N_v_hot  +
            N_ec_hot +
            N_b_hot  +
            N_ic_hot
        );  // mol
        
        // Convert to Bq (decays/s per particle)
        double activityBq = N_hot * SimParams.LAMBDA_DECAY * SimParams.AVOGADRO;  // mol* 1/s * (1/mol) = particles/s = Bq
        
        // Energy emission rate (J/s) - before escape
		double energyRate_emitted = activityBq * SimParams.E_BETA_LU177;  // Bq * J = J/s

        // β particle escape correction
        // Lu-177 beta range ~1 mm in tissue
        double beta_range_m = 1.0e-3;  // m
		// Assume the tumour cells are in the extruded cylinder shape (sphere gives the same formula)
		// with an effective radius
        double R_tumor_eff = Math.pow(3.0 * tumorVolume / (4.0 * Math.PI), 1.0/3.0);
        double f_deposit = Math.pow(R_tumor_eff / (R_tumor_eff + beta_range_m), 3.0);
		// Spread the energy over a sphere of radius R_tumor_eff + beta_range_m
        double energyRate = energyRate_emitted * f_deposit;
                
        // Tumor mass (assuming density = 1000 kg/m³)
        double tumorMass = tumorVolume * 1000.0;  // m³ * kg/m³ = kg
        
        // Dose rate in Gy/s
        double doseRate_Gy_per_s = (tumorMass > 1e-15) 
                                   ? energyRate / tumorMass 
                                   : 0.0;
        
        // Convert to Gy/h for compatibility with old outputs
        double doseRate_Gy_per_h = doseRate_Gy_per_s * 3600.0;

        // DIAGNOSTIC: Print every 10 days
        if (currentDay % 10 == 0 && currentHour == 0 && SimParams.VERBOSE_ON) {
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
                             double N_cen_hot, double N_cen_cold,
                             double N_v_hot, double N_v_cold,
                             double N_ec_hot, double N_ec_cold,
                             double N_b_hot, double N_b_cold,
                             double N_ic_hot, double N_ic_cold) {
        File logFile = new File("Results/PBPK_log.csv");

        try (FileWriter fw = new FileWriter(logFile, true);
             BufferedWriter bw = new BufferedWriter(fw);
             PrintWriter out = new PrintWriter(bw)) {
    
            // Header if file is empty
            if (logFile.length() == 0) {
                out.println("day,hour,time_days,N_cen_hot,N_cen_cold,N_v_hot,N_v_cold," +
                           "N_ec_hot,N_ec_cold,N_b_hot,N_b_cold,N_ic_hot,N_ic_cold");
            }
    
            double timeDays = day + hour / 24.0;
            
            // Convert to nmol for logging (easier to read)
            double conv = 1e9;  // mol → nmol (factor of 1e9)
            
            out.printf("%.0f,%.0f,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e%n",
                (double) day, (double) hour, timeDays,
                N_cen_hot * conv, N_cen_cold * conv, 
                N_v_hot * conv, N_v_cold * conv, 
                N_ec_hot * conv, N_ec_cold * conv, 
                N_b_hot * conv, N_b_cold * conv,
                N_ic_hot * conv, N_ic_cold * conv
            );
    
        } catch (IOException e) {
            System.err.println("Error writing PBPK log: " + e.getMessage());
        }
    }
}
