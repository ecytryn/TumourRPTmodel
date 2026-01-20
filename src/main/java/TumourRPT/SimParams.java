package TumorRPT;

// File I/O for writing reports
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

// Timestamps
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

// Array formatting
import java.util.Arrays;

/**
 * All simulation parameters in SI units
 * 
 * UNITS USED THROUGHOUT CODE:
 * - Time: seconds (s)
 * - Length: meters (m)  
 * - Amount: moles (mol)
 * - Concentration: mol/m³
 * - Volume: m³
 * - Pressure: Pascals (Pa)
 * - Mass: kilograms (kg)
 * - Energy: Joules (J)
 * 
 * NO UNIT CONVERSIONS NEEDED IN SIMULATION CODE!
 * All parameters are already in the correct SI units.
 */
public class SimParams {

	// Define the output string using the purpose of the run	
	private static final String EXPERIMENT_NAME_RAW = "";  // ← Set this for real runs
	private static final String EXPERIMENT_DESCRIPTION_RAW = "";
//	private static final String EXPERIMENT_NAME_RAW = "RPTonNormoxicTumour";  // ← Set this for real runs
//	private static final String EXPERIMENT_DESCRIPTION_RAW = "This is the first run of a two-run experiment. Here, I hit the tumour with a dose of RPT at Day 5 when the tumour is mostly normoxic.";
//	private static final String EXPERIMENT_NAME_RAW = "RPTonHypoxicTumour";  // ← Set this for real runs
//	private static final String EXPERIMENT_DESCRIPTION_RAW = "This is the second run of a two-run experiment. Here, I hit the tumour with a dose of RPT at Day 35 when a large portion of the tumour is hypoxic.";

	// Auto-default to "testing_run" if empty
	public static final String EXPERIMENT_NAME = 
		EXPERIMENT_NAME_RAW.isEmpty() ? "testing_run" : EXPERIMENT_NAME_RAW;
	
	public static final String EXPERIMENT_DESCRIPTION = 
		EXPERIMENT_DESCRIPTION_RAW.isEmpty() 
			? "Quick test/debug run" 
			: EXPERIMENT_DESCRIPTION_RAW;

    // Output directories (set by Main at startup)
    public static String OUTPUT_DIR_BASE = "";
    public static String OUTPUT_DIR_TUMOUR_IMAGES = "";
    public static String OUTPUT_DIR_OXYGEN_IMAGES = "";
    
    // =================================================================
    // RUNTIME CONFIGURATION
    // =================================================================
    public static final boolean VERBOSE_ON = false;  // Set true for reporting to the terminal
    public static final boolean EXPORT_OX_IMAGES = true;  // Set true for exporting oxygen images
    public static final boolean EXPORT_CONSUMP_IMAGES = false;  // Set true for exporting consumption images
    public static final boolean EXPORT_TUMOUR_OX_IMAGES = true; // Set true for exporting tumour/oxygen images
    public static final boolean PLOT_LIVE_IMAGES = false;  // Set true to pop up images during the run
    public static final boolean ENABLE_PBPK_LOGGING = false;   // Set true for PK debugging
	public static final boolean FREEZE_TUMOR = false;          // Set true to test PK without cell dynamics    
	public static final boolean USE_CALIBRATED_BC = true;  // true --> use a BC for the O2 PDE that is calibrated to average O2 levels in healthy tissue
    
    // =================================================================
    // PHYSICAL CONSTANTS
    // =================================================================
    public static final double AVOGADRO = 6.022e23;  // molecules/mol
    public static final double MMHG_TO_PA = 133.322; // Pa/mmHg
    
    // =================================================================
    // DOMAIN & DISCRETIZATION
    // =================================================================
    public static final double CELL_LENGTH = 1e-5;     // m (10 μm)
    public static final int GRID_SIZE = 400;            // cells per side
    public static final double DOMAIN_SIZE = GRID_SIZE * CELL_LENGTH;  // m (4 mm)
    public static final double TIME_STEP = 3600.0;      // s (1 hour)
    public static final double CELL_CYCLE = 86400.0;    // s (1 day)
	public static int globalTime = 0;  // Current simulation time in hours (updated during run)    
	// These are just rescalings that get used in a couple places.
	public static final double CELL_CYCLE_IN_HRS = CELL_CYCLE / 3600.0;  // 86400 s → 24 hours
	public static final double TIME_STEP_IN_HRS = TIME_STEP / 3600.0;    // 3600 s → 1 hour
    // =================================================================
    // PHARMACOKINETICS - Rate Constants
    // =================================================================
    // All in 1/s (per second) - written as per minute values with conversion factors to SI units
    public static final double LAMBDA_BIO = 1.6e-4 / 60.0;    // 1/s (biological clearance)
    public static final double LAMBDA_DECAY = 7.14e-5 / 60.0; // 1/s (Lu-177 physical decay)
    public static final double K_ON = 0.046 * 1e6 / 60.0;     // m³/(mol·s) (binding) (0.046 lit/nmol/min)
    public static final double K_OFF = 0.368 / 60.0;          // 1/s (unbinding)
    public static final double K_INT = 0.001 / 60.0;          // 1/s (internalization)
    public static final double K_REL = 2e-4 / 60.0;           // 1/s (release from cells)
    
    // =================================================================
    // PHARMACOKINETICS - Volumes
    // =================================================================
    public static final double V_CENTRAL = 0.458e-3;          // m³ (central/arterial, 0.458 L)
    public static final double INTERSTITIAL_FRACTION = 0.4;   // fraction of tumor that is interstitial
    
    // Per-vessel parameters
    public static final double VESSEL_VOLUME = Math.pow(CELL_LENGTH, 2) * DOMAIN_SIZE;  // m³ (cross-section × height)
    public static final double VESSEL_FLOW = 5e-9 / 1000.0;   // m³/s (5 nL/s → m³/s)
    public static final double VESSEL_PS = (6.0/5.0) * VESSEL_FLOW;  // m³/s (permeability-surface)
    
    // Vessel influence for PK calculations
	public static final int VESSEL_INFLUENCE_RADIUS = 10;  // cells (~100 μm)
	
    // =================================================================
    // CELL BIOLOGY
    // =================================================================
    public static final double RECEPTORS_PER_CELL = 5e5;      // receptors/cell
    public static final double RECEPTORS_PER_CELL_MOL = RECEPTORS_PER_CELL / AVOGADRO;  // mol/cell
    
    public static final double INITIAL_TUMOR_RADIUS = 333e-6;  // m (1 mm)
    public static final double INITIAL_TUMOR_RADIUS_CELLS = INITIAL_TUMOR_RADIUS/CELL_LENGTH;  // m (1 mm)
    
    // Cell removal times
    public static final double NECROTIC_REMOVAL_TIME = 10.0 * 86400.0;  // s (10 days)
    public static final double APOPTOTIC_REMOVAL_TIME = 2.0 * 86400.0;  // s (2 days)
	public static final double ApopRemovalTime = APOPTOTIC_REMOVAL_TIME / 86400.0;  // Convert s→days
	public static final double APOP_REMOVAL_PROB_PER_HOUR = // <-- needs to be fixed if time step is ever changed!!
		(1.0 / APOPTOTIC_REMOVAL_TIME) * 3600.0;  // (1/s) * (s/hour) = 1/hour 
	public static final double strictRemovalCoeff = 1.0;


	public static final double DIVISON_PROB_MAX = 0.5;  // per day
	public static final int[] divHood = HAL.Util.CircleHood(true, 1);  // immediate neighbors    

    // =================================================================
    // OXYGEN DIFFUSION
    // =================================================================
    public static final double D_O2 = 2e-5 * 1e-4;            // m²/s (2e-5 cm²/s)
    
    // Oxygen partial pressures
    public static final double P_O2_VESSEL = 100.0 * MMHG_TO_PA;   // Pa (100 mmHg)
    public static final double P_O2_VEIN = 30.0 * MMHG_TO_PA;      // Pa (30 mmHg)
    public static final double P_O2_HYPOXIC = 10.0 * MMHG_TO_PA;   // Pa (10 mmHg)
    public static final double P_O2_NECROTIC = 0.5 * MMHG_TO_PA;   // Pa (0.5 mmHg)
    
    // Henry's constant for O2
    public static final double HENRY_O2 = 1.3e-3 / 760.0;     // mol/(L·atm) = 1.71e-6 mol/(L·Pa)
    // Convert to mol/(m³·Pa): multiply by 1000 L/m³
    public static final double HENRY_O2_SI = HENRY_O2 * 1000.0;  // mol/(m³·Pa)
    
    // Oxygen concentration thresholds
    public static final double C_O2_VESSEL = P_O2_VESSEL * HENRY_O2_SI;   // mol/m³
    public static final double C_O2_VEIN = P_O2_VEIN * HENRY_O2_SI;       // mol/m³
    public static final double C_O2_HYPOXIC = P_O2_HYPOXIC * HENRY_O2_SI; // mol/m³
    public static final double C_O2_NECROTIC = P_O2_NECROTIC * HENRY_O2_SI; // mol/m³

	//Consumption rate constants (1/s) - used in oxygen PDE: -C(x,y)·u    
    public static final double CONSUMPTION_HEALTHY = 1500.0/3600;   // 
    public static final double CONSUMPTION_NORMAL = 7500.0/3600;    // 
    public static final double CONSUMPTION_HYPOXIC = 3000.0/3600;   // 
    public static final double CONSUMPTION_NECROTIC = 1.0/3600;     // 
    public static final double CONSUMPTION_APOPTOTIC = 1.0/3600;    // 
    public static final double CONSUMPTION_VESSEL = 0.0;       // no consumption

	// Convert your individual constants to array:
	public static final double[] CELLS_CONSUMPTION_RATE_LIST = {
		CONSUMPTION_HEALTHY,   // 0
		CONSUMPTION_NORMAL,    // 1
		CONSUMPTION_HYPOXIC,   // 2
		CONSUMPTION_NECROTIC,  // 3
		CONSUMPTION_APOPTOTIC, // 4
		CONSUMPTION_VESSEL     // 5
	};    

    // =================================================================
    // RADIOBIOLOGY
    // =================================================================
//    public static final double ALPHA_NORMAL = 0.15;           // Gy⁻¹
//    public static final double ALPHA_HYPOXIC = 0.0107;        // Gy⁻¹ (10× less sensitive)
//    public static final double BETA_NORMAL = 0.048;           // Gy⁻²
//    public static final double BETA_HYPOXIC = 0.0024;         // Gy⁻² (10× less sensitive)
//    public static final double REPAIR_RATE = 0.3 / 3600.0;    // 1/s (0.3 per hour)

	public static final double ALPHA_NORMAL = 0.9 * 0.15;           // Gy⁻¹
    public static final double ALPHA_HYPOXIC = 0.06;          // Gy⁻¹ (0.15/2.5, OER=2.5)
    public static final double BETA_NORMAL = 0.048;           // Gy⁻²
    public static final double BETA_HYPOXIC = 0.019;          // Gy⁻² (0.048/2.5, OER=2.5)
    public static final double REPAIR_RATE = 0.65 / 3600.0;   // 1/s (0.65 per hour)    

    public static final double E_BETA_LU177 = 2.14e-14;       // J (average beta energy)

	public static final double repairRate = REPAIR_RATE * 3600.0;  // Convert 1/s → 1/hour
	
	public static final double[] alphaValues = {
		ALPHA_NORMAL,   // index 0: normoxic
		ALPHA_HYPOXIC   // index 1: hypoxic
	};
	
	public static final double[] betaValues = {
		BETA_NORMAL,    // index 0: normoxic
		BETA_HYPOXIC    // index 1: hypoxic
	};
	
	public static final int maxLookupAge = 200;  // hours
    
    // =================================================================
    // INJECTION PROTOCOL (defaults, can be overridden)
    // =================================================================
    public static final double DOSE_PER_INJECTION = 100e-9;   // mol (100 nmol)
    public static final double HOT_FRACTION = 0.1;            // 10% hot, 90% cold
    public static final int[] DEFAULT_INJECTION_DAYS = {5, 35, 65, 95};  // days
    
    // =================================================================
    // CELL TYPES (indices)
    // =================================================================
    public static final int HEALTHY = 0;
    public static final int NORMAL = 1;
    public static final int HYPOXIC = 2;
    public static final int NECROTIC = 3;
    public static final int APOPTOTIC = 4;
    public static final int VESSEL = 5;
    public static final int NUM_CELL_TYPES = 6;
    	
    // Vessel density configuration
    // Possibilities: "uniform", "sparse", "dense", "heterogeneous", etc.
    // Corresponds to CSV files in src/main/resources/vasculature/
    // Use scripts/GenerateUniformVessels2.py to create new csv files.
	public static final String VESSEL_DENSITY_CONFIG = "uniform";  // Default configuration
    // Optional: Numerical density for validation/reporting - not implemented
//    public static final double VESSEL_DENSITY_PER_MM2 = 100.0;  // vessels/mm² (for uniform)

    // =================================================================
    // DERIVED/COMPUTED PARAMETERS
    // =================================================================
    
    /**
     * Compute total tumor volume from cell count
     * @param cellCount Number of cells in 2D cross-section
     * @return Volume in m³
     */
    public static double computeTumorVolume(int cellCount) {
        double radius = Math.sqrt(cellCount / Math.PI) * CELL_LENGTH;
        double height = 2.0 * radius;  // cylindrical extrusion
        return Math.PI * radius * radius * height;
    }
    
    /**
     * Compute total receptors in tumor (accounting for 3D extrusion)
     * @param cellCount2D Number of cells in 2D cross-section
     * @param numVessels Number of vessels (to subtract)
     * @return Total receptors in moles
     */
	public static double computeReceptorMoles(int tumorCells2D, int numVessels) {
		// Parameter already excludes vessels - no subtraction needed!
		if (tumorCells2D <= 0) {
			return 0.0;
		}
		
		double radius = Math.sqrt(tumorCells2D / Math.PI);  // in cell lengths
		double height = 2.0 * radius;						// cylindrical extension
		double totalCells3D = tumorCells2D * height;
		return totalCells3D * RECEPTORS_PER_CELL_MOL;
	}    
    /**
     * Compute interstitial volume from given tumour volume
     * @param domainVolume Total tumour volume in m³
     * @return Interstitial volume in m³
     */
    public static double computeInterstitialVolume(double tumourVolume) {
        return INTERSTITIAL_FRACTION * tumourVolume;
    }

	/**
	 * Update global time from current day and hour
	 * @param currentDay Current simulation day
	 * @param currentHour Current hour (0-23)
	 */
	public static void updateGlobalTime(int currentDay, int currentHour) {
		globalTime = currentDay * 24 + currentHour;
	}
    // =================================================================
	// VISUALIZATION
    // =================================================================
/**
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFF69B4,  // 1: Normal (hot pink)
		0xFF4169E1,  // 2: Hypoxic (royal blue)
		0xFF8B4513,  // 3: Necrotic (saddle brown)
		0xFFFFFF00,  // 4: Apoptotic (yellow)
		0xFFFF0000   // 5: Vessel (red)
	};
*/
	public static final int[] COLORLIST = {
			0x88888888,  // 0: Healthy (gray)
			0xFF90EE90,  // 1: Normal tumor (light green)
			0xFF228B22,  // 2: Hypoxic (forest green)
			0xFF006400,  // 3: Necrotic (dark green)
			0xFF800080,  // 4: Apoptotic (purple)
			0xFFFF6347   // 5: Vessel (tomato red)
	};	

	/**
	 * Export an md file for human readable parameter summary
	 */
	public static void generateParameterReport(String filepath) throws IOException {
		try (PrintWriter out = new PrintWriter(new FileWriter(filepath))) {
			out.println("# Simulation Parameters Report");
			out.println();
			out.println("**Experiment:** " + EXPERIMENT_NAME);
			out.println("**Date:** " + java.time.LocalDateTime.now());
			out.println("**Description:** " + EXPERIMENT_DESCRIPTION);
			out.println();
			
			out.println("## Injection Protocol");
			out.println("| Parameter | Value |");
			out.println("|-----------|-------|");
			out.println("| Days | " + java.util.Arrays.toString(DEFAULT_INJECTION_DAYS) + " |");
			out.println("| Dose | " + (DOSE_PER_INJECTION * 1e9) + " nmol |");
			out.println("| Hot fraction | " + (HOT_FRACTION * 100) + "% |");
			out.println();
			
			out.println("## Domain");
			out.println("| Parameter | Value |");
			out.println("|-----------|-------|");
			out.println("| Grid size | " + GRID_SIZE + " cells |");
			out.println("| Cell length | " + (CELL_LENGTH * 1e6) + " μm |");
			out.println("| Domain size | " + (DOMAIN_SIZE * 1e3) + " mm |");
			out.println();
			
			out.println("## Pharmacokinetics");
			out.println("| Parameter | Value | Description |");
			out.println("|-----------|-------|-------------|");
			out.println("| λ_bio | " + (LAMBDA_BIO * 3600) + " hr⁻¹ | Biological clearance |");
			out.println("| λ_decay | " + (LAMBDA_DECAY * 3600) + " hr⁻¹ | Lu-177 decay |");
			out.println("| k_on | " + (K_ON / 1e6 * 60) + " L/(nmol·min) | Binding rate |");
			out.println("| k_off | " + (K_OFF * 60) + " min⁻¹ | Unbinding rate |");
			out.println("| k_int | " + (K_INT * 60) + " min⁻¹ | Internalization |");
			out.println();
			
			out.println("## Radiobiology");
			out.println("| Cell Type | α (Gy⁻¹) | β (Gy⁻²) |");
			out.println("|-----------|----------|----------|");
			out.println("| Normoxic | " + ALPHA_NORMAL + " | " + BETA_NORMAL + " |");
			out.println("| Hypoxic | " + ALPHA_HYPOXIC + " | " + BETA_HYPOXIC + " |");
			out.println("| Repair rate | " + (REPAIR_RATE * 3600) + " hr⁻¹ | |");
			out.println();
			
			out.println("## Oxygen");
			out.println("| Threshold | Value (mmHg) |");
			out.println("|-----------|--------------|");
			out.println("| Vessel | " + (P_O2_VESSEL / MMHG_TO_PA) + " |");
			out.println("| Hypoxic | " + (P_O2_HYPOXIC / MMHG_TO_PA) + " |");
			out.println("| Necrotic | " + (P_O2_NECROTIC / MMHG_TO_PA) + " |");
			
			System.out.println("Parameter report generated: " + filepath);
		}
	}
   

   
	/**
	 * Export key simulation parameters to CSV for reproducibility
	 */
	public static void exportParametersToCSV(String filepath) throws IOException {
		try (PrintWriter out = new PrintWriter(new FileWriter(filepath))) {
			out.println("parameter,value,units,description");
			
			// Experiment metadata
			out.println("experiment_name," + EXPERIMENT_NAME + ",,Experiment identifier");
			out.println("experiment_date," + java.time.LocalDateTime.now() + ",,Run timestamp");
			
			// Injection protocol
			out.println("injection_days," + java.util.Arrays.toString(DEFAULT_INJECTION_DAYS) + ",days,");
			out.println("dose_per_injection," + DOSE_PER_INJECTION * 1e9 + ",nmol,");
			out.println("hot_fraction," + HOT_FRACTION + ",,");
			
			// Grid
			out.println("grid_size," + GRID_SIZE + ",cells,");
			out.println("cell_length," + CELL_LENGTH * 1e6 + ",um,");
			out.println("domain_size," + DOMAIN_SIZE * 1e3 + ",mm,");
			
			// PK
			out.println("lambda_bio," + LAMBDA_BIO * 3600 + ",1/hr,Biological clearance");
			out.println("lambda_decay," + LAMBDA_DECAY * 3600 + ",1/hr,Lu-177 decay");
			out.println("k_on," + K_ON / 1e6 * 60 + ",L/(nmol*min),Binding rate");
			out.println("k_off," + K_OFF * 60 + ",1/min,Unbinding rate");
			out.println("k_int," + K_INT * 60 + ",1/min,Internalization rate");
			out.println("v_central," + V_CENTRAL * 1e3 + ",L,Central volume");
			
			// Radiobiology
			out.println("alpha_normoxic," + ALPHA_NORMAL + ",Gy^-1,");
			out.println("beta_normoxic," + BETA_NORMAL + ",Gy^-2,");
			out.println("alpha_hypoxic," + ALPHA_HYPOXIC + ",Gy^-1,");
			out.println("beta_hypoxic," + BETA_HYPOXIC + ",Gy^-2,");
			out.println("repair_rate," + REPAIR_RATE * 3600 + ",1/hr,");
			
			// Oxygen
			out.println("p_o2_vessel," + P_O2_VESSEL / MMHG_TO_PA + ",mmHg,");
			out.println("p_o2_hypoxic," + P_O2_HYPOXIC / MMHG_TO_PA + ",mmHg,");
			out.println("p_o2_necrotic," + P_O2_NECROTIC / MMHG_TO_PA + ",mmHg,");
			out.println("d_o2," + D_O2 * 1e4 + ",cm^2/s,");
			
			// Cell biology
			out.println("receptors_per_cell," + RECEPTORS_PER_CELL + ",,");
			out.println("cell_cycle," + CELL_CYCLE / 3600 + ",hours,");
			
			System.out.println("Parameters exported to: " + filepath);
		}
	}


}
