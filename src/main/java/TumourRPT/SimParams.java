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
 * - Concentration: mol/m^3
 * - Volume: m^3
 * - Pressure: Pascals (Pa)
 * - Mass: kilograms (kg)
 * - Energy: Joules (J)
 * 
 * NO UNIT CONVERSIfONS NEEDED IN SIMULATION CODE!
 * All parameters are already in the correct SI units.
 */
public class SimParams {
    // =================================================================
    // Default parameters that can be overwritten for different numerical experiments
    // =================================================================

	public static String EXPERIMENT_NAME = "default";
	public static String EXPERIMENT_DESCRIPTION = "Default run";
	public static int[] INJECTION_SCHEDULE = {45};
	public static double INITIAL_TUMOR_RADIUS = 333e-6;
	public static double ALPHA_HYPOXIC = 0.06;
	public static double BETA_HYPOXIC = 0.019;
    // Vessel density configuration
    // Use scripts/GenerateVessels/GenerateUniformVessels.py to create new csv files.
    // Possibilities: 40, 50, 60 - these refer to the R_REPEL parameter used in GenerateUniformVessels.py
    // Density of vessels is "inversely" proportional to R_REPEL
    // Corresponds to CSV files in src/main/resources/vasculature/
	// R_REPEL = 40 --> 616.1 vessels/mm^2 
	// R_REPEL = 50 --> 604.8 vessels/mm^2 
	// R_REPEL = 60 --> 591.9 vessels/mm^2 
	public static String VESSEL_DENSITY_CONFIG = "50";  // Default configuration

	
	// Add setter method
	public static void setExperiment(String name, String desc, int[] injections, 
									  double radius, double alphaH, double betaH) {
		EXPERIMENT_NAME = name;
		EXPERIMENT_DESCRIPTION = desc;
		INJECTION_SCHEDULE = injections;
		INITIAL_TUMOR_RADIUS = radius;
		ALPHA_HYPOXIC = alphaH;
		BETA_HYPOXIC = betaH;
	}

    // Output directories (set by Main at startup)
    public static String OUTPUT_DIR_BASE = "";
    public static String OUTPUT_DIR_TUMOUR_IMAGES = "";
    public static String OUTPUT_DIR_OXYGEN_IMAGES = "";
	public static String OUTPUT_DIR_SF_IMAGES = "";
    
    // =================================================================
    // RUNTIME CONFIGURATION
    // =================================================================
	
	// For generating manuscript figures
    public static boolean EXPORT_OX_IMAGES = true;  // Set true for exporting oxygen images
    public static boolean EXPORT_TUMOUR_OX_IMAGES = true; // Set true for exporting tumour/oxygen images

	// for testing and debugging
    public static boolean EXPORT_CONSUMP_IMAGES = false;  // Set true for exporting consumption images
	public static boolean EXPORT_SF_IMAGES = false;  // Set true for SF visualization    
    public static final boolean VERBOSE_ON = false;  // Set true for reporting to the terminal
	public static final boolean FREEZE_TUMOR = false;          // Set true to test PK without cell dynamics    

	// Likely to forever remain at these setting
	public static final boolean PLOT_LIVE_IMAGES = false;  // Set true to pop up images during the run
	public static final boolean USE_CALIBRATED_BC = true;  // true --> use a BC for the O2 PDE that is calibrated to average O2 levels in healthy tissue
    public static final boolean ENABLE_PBPK_LOGGING = false;   // Set true for PK debugging
    
    // =================================================================
    // PHYSICAL CONSTANTS
    // =================================================================
    public static final double AVOGADRO = 6.022e23;  // molecules/mol
    public static final double MMHG_TO_PA = 133.322; // Pa/mmHg
    
    // =================================================================
    // DOMAIN & DISCRETIZATION
    // =================================================================
    public static final double CELL_LENGTH = 1e-5;     // m (10 um)
    public static final int GRID_SIZE = 400;            // cells per side
    public static final double DOMAIN_SIZE = GRID_SIZE * CELL_LENGTH;  // m (4 mm)
    public static final double TIME_STEP = 3600.0;      // s (1 hour)
    public static final double CELL_CYCLE = 86400.0;    // s (1 day)
	public static int globalTime = 0;  // Current simulation time in hours (updated during run)    
	// These are just rescalings that get used in a couple places.
	public static final double CELL_CYCLE_IN_HRS = CELL_CYCLE / 3600.0;  // 86400 s -> 24 hours
	public static final double TIME_STEP_IN_HRS = TIME_STEP / 3600.0;    // 3600 s -> 1 hour
    // =================================================================
    // PHARMACOKINETICS - Rate Constants
    // =================================================================
    // All in 1/s (per second) - written as per minute values with conversion factors to SI units
    public static final double LAMBDA_BIO = 1.6e-4 / 60.0;    // 1/s (biological clearance)
    public static final double LAMBDA_DECAY = 7.14e-5 / 60.0; // 1/s (Lu-177 physical decay)
    public static final double K_ON = 0.046 * 1e6 / 60.0;     // m^3/(mol s) (binding) (0.046 lit/nmol/min)
    public static final double K_OFF = 0.368 / 60.0;          // 1/s (unbinding)
    public static final double K_INT = 0.001 / 60.0;          // 1/s (internalization)
    public static final double K_REL = 2e-4 / 60.0;           // 1/s (release from cells)
    
    // =================================================================
    // PHARMACOKINETICS - Volumes
    // =================================================================
    public static final double V_CENTRAL = 0.458e-3;          // m^3 (central/arterial, 0.458 L)
    public static final double INTERSTITIAL_FRACTION = 0.4;   // fraction of tumor that is interstitial
    
    // Per-vessel parameters
    public static final double VESSEL_VOLUME = Math.pow(CELL_LENGTH, 2) * DOMAIN_SIZE;  // m^3 (cross-section × height)
    public static final double VESSEL_FLOW = 5e-9 / 1000.0;   // m^3/s (5 nL/s -> m^3/s)
    public static final double VESSEL_PS = (6.0/5.0) * VESSEL_FLOW;  // m^3/s (permeability-surface)
    
    // Vessel influence for PK calculations
	public static final int VESSEL_INFLUENCE_RADIUS = 10;  // cells (~100 um)
	
    // =================================================================
    // CELL BIOLOGY
    // =================================================================
    public static final double RECEPTORS_PER_CELL = 4e5;      // receptors/cell
    public static double RECEPTORS_PER_CELL_MOL = RECEPTORS_PER_CELL / AVOGADRO;  // mol/cell

	public static double getInitialTumorRadiusCells() {
		return INITIAL_TUMOR_RADIUS / CELL_LENGTH;
	}    
    // Cell removal times
    public static final double NECROTIC_REMOVAL_TIME = 10.0 * 86400.0;  // s (10 days)
    public static final double APOPTOTIC_REMOVAL_TIME = 2.0 * 86400.0;  // s (2 days)
	public static final double ApopRemovalTime = APOPTOTIC_REMOVAL_TIME / 86400.0;  // Convert s->days
	public static final double APOP_REMOVAL_PROB_PER_HOUR = // <-- needs to be fixed if time step is ever changed!!
		(1.0 / APOPTOTIC_REMOVAL_TIME) * 3600.0;  // (1/s) * (s/hour) = 1/hour 
	public static final double strictRemovalCoeff = 1.0;


	public static final double DIVISION_PROB_MAX = 0.5;  // per day
	public static final int[] divHood = HAL.Util.CircleHood(true, 1);  // immediate neighbors    

    // =================================================================
    // OXYGEN DIFFUSION
    // =================================================================
    public static final double D_O2 = 2e-5 * 1e-4;            // m^2/s (2e-5 cm^2/s)
    
    // Oxygen partial pressures
    public static final double P_O2_VESSEL = 100.0 * MMHG_TO_PA;   // Pa (100 mmHg)
    public static final double P_O2_VEIN = 30.0 * MMHG_TO_PA;      // Pa (30 mmHg)
    public static final double P_O2_HYPOXIC = 10.0 * MMHG_TO_PA;   // Pa (10 mmHg)
    public static final double P_O2_NECROTIC = 0.5 * MMHG_TO_PA;   // Pa (0.5 mmHg)
    
    // Henry's constant for O2
    public static final double HENRY_O2 = 1.3e-3 / 760.0;     // mol/(L atm) = 1.71e-6 mol/(L Pa)
    // Convert to mol/(m^3 Pa): multiply by 1000 L/m^3
    public static final double HENRY_O2_SI = HENRY_O2 * 1000.0;  // mol/(m^3 Pa)
    
    // Oxygen concentration thresholds
    public static final double C_O2_VESSEL = P_O2_VESSEL * HENRY_O2_SI;   // mol/m^3
    public static final double C_O2_VEIN = P_O2_VEIN * HENRY_O2_SI;       // mol/m^3
    public static final double C_O2_HYPOXIC = P_O2_HYPOXIC * HENRY_O2_SI; // mol/m^3
    public static final double C_O2_NECROTIC = P_O2_NECROTIC * HENRY_O2_SI; // mol/m^3

	//Consumption rate constants (1/s) - used in oxygen PDE: -C(x,y)u    
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

	public static final double ALPHA_NORMAL = 0.15;           // Gy^(-1)
    public static final double BETA_NORMAL = 0.048;           // Gy^(-2)
/** These are now set in Main.
    public static final double ALPHA_HYPOXIC = 0.06;          // Gy^(-1) (0.15/2.5, OER=2.5)
    public static final double BETA_HYPOXIC = 0.019;          // Gy^(-2) (0.048/2.5, OER=2.5)
*/

    public static final double REPAIR_RATE = 0.65 / 3600.0;   // 1/s (0.65 per hour)    

    public static final double E_BETA_LU177 = 2.14e-14;       // J (average beta energy)

	public static final double repairRate = REPAIR_RATE * 3600.0;  // Convert 1/s -> 1/hour
	
	public static final double[] alphaValues = {
		ALPHA_NORMAL,   // index 0: normoxic
		ALPHA_HYPOXIC   // index 1: hypoxic
	};
	
	public static final double[] betaValues = {
		BETA_NORMAL,    // index 0: normoxic
		BETA_HYPOXIC    // index 1: hypoxic
	};
	
	public static final int maxLookupAge = 200;  // days
    
    // =================================================================
    // INJECTION PROTOCOL
    // =================================================================
    // **FOR SINGLE RUNS:** Change these parameters to run individual experiments
    // **FOR PARAMETER SWEEPS:** Use IntervalSkewSweep.java which programmatically 
    //                          varies inter-injection interval and dose distribution
    //                           or DoseReceptorSweep.java (injected amount and receptor density)
    // Main.java uses these values for single-run mode
    
    public static final double DOSE_PER_INJECTION = 100e-9;   // mol (100 nmol per injection is baseline)
    public static final double HOT_FRACTION = 0.1;            // Fraction that is radioactive (0.1 = 10%)
/**	Moved to Main
    public static final int[] INJECTION_SCHEDULE = {5};  // Days to inject (e.g., {5, 35, 65, 95})
*/    
    // Examples for different treatment schedules:
    // Single dose:      {35}
    // 2 doses:          {5, 35}
    // 4 doses:          {5, 35, 65, 95}
    // Weekly x 4:       {7, 14, 21, 28}
    
    // Simulation length (days to simulate after last injection)
    public static final int DAYS_AFTER_LAST_INJECTION = 30;
    
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

    // =================================================================
    // DERIVED/COMPUTED PARAMETERS
    // =================================================================
    
    /**
     * Compute total tumor volume from cell count
     * @param cellCount Number of cells in 2D cross-section
     * @return Volume in m^3
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
     * @param domainVolume Total tumour volume in m^3
     * @return Interstitial volume in m^3
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
// My original colour-blind unfriendly colour scheme
public static final int[] COLORLIST = {
			0x88888888,  // 0: Healthy (gray)
			0xFF90EE90,  // 1: Normal tumor (light green)
			0xFF228B22,  // 2: Hypoxic (forest green)
			0xFF006400,  // 3: Necrotic (dark green)
//			0xFF800080,  // 4: Apoptotic (purple)
			0xFFD0A000,  // 4: Apoptotic (orange)
			0xFFFF6347,  // 5: Vessel (tomato red)
			0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};	
	
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFF90EE90,  // 1: Normal (light green - RGB 0.56, 0.93, 0.56)
		0xFF3CB371,  // 2: Hypoxic (medium sea green - RGB 0.24, 0.70, 0.44)
		0xFF006400,  // 3: Necrotic (dark green - RGB 0, 0.39, 0)
		0xFFDA70D6,  // 4: Apoptotic (orchid/bright purple - RGB 0.85, 0.44, 0.84)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFF98FB98,  // 1: Normal (pale green - RGB 0.60, 0.98, 0.60)
		0xFF2E8B57,  // 2: Hypoxic (sea green - RGB 0.18, 0.55, 0.34)
		0xFF0B4F0B,  // 3: Necrotic (very dark green - RGB 0.04, 0.31, 0.04)
		0xFFBF40BF,  // 4: Apoptotic (medium-bright purple - RGB 0.75, 0.25, 0.75)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};

	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFF7FFFD4,  // 1: Normal (aquamarine - RGB 0.50, 1.0, 0.83)
		0xFF20B2AA,  // 2: Hypoxic (light sea green - RGB 0.13, 0.70, 0.67)
		0xFF00695C,  // 3: Necrotic (teal green - RGB 0, 0.41, 0.36)
		0xFFFF1493,  // 4: Apoptotic (deep pink/magenta - RGB 1.0, 0.08, 0.58)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};

	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFF66FF66,  // 1: Normal (RGB 0.4, 1.0, 0.4 - bright light green)
		0xFF00CC00,  // 2: Hypoxic (RGB 0, 0.8, 0 - strong green)
		0xFF004D00,  // 3: Necrotic (RGB 0, 0.3, 0 - dark forest green)
		0xFFB266FF,  // 4: Apoptotic (RGB 0.7, 0.4, 1.0 - bright purple)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};

	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFD2B48C,  // 1: Normal tumor (light tan)
		0xFFB87333,  // 2: Hypoxic (copper)
		0xFF654321,  // 3: Necrotic (dark brown)
		0xFFB8A8C8,  // 4: Apoptotic (soft lavender)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};


	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFFE87C,  // 1: Normal tumor (light yellow)
		0xFFFF8C00,  // 2: Hypoxic (dark orange)
		0xFF8B4513,  // 3: Necrotic (saddle brown)
		0xFFFF00FF,  // 4: Apoptotic (magenta - highly distinct)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};

// A "flip" of 4th above
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFFB3D9,  // 1: Normal (light pink - lightened deep pink)
		0xFFFF1493,  // 2: Hypoxic (deep pink/magenta - your reference color)
		0xFFC70069,  // 3: Necrotic (dark magenta - darkened deep pink)
		0xFF20B2AA,  // 4: Apoptotic (light sea green/teal - your reference color)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark maroon)
	};
	
// A "flip" of 4th above with a wider spacing between norm-hypo-necro
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFFCCE5,  // 1: Normal (very light pink - more pastel)
		0xFFFF1493,  // 2: Hypoxic (deep pink/magenta)
		0xFF99004D,  // 3: Necrotic (very dark magenta)
		0xFF20B2AA,  // 4: Apoptotic (light sea green/teal)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark maroon)
	};
// *** Slightly darker teal - apoptotic  cells compared to the previous scheme ^
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFFCCE5,  // 1: Normal (very light pink)
		0xFFFF1493,  // 2: Hypoxic (deep pink/magenta)
		0xFF99004D,  // 3: Necrotic (very dark magenta)
		0xFF008B8B,  // 4: Apoptotic (dark cyan - deeper than light sea green)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark maroon)
	};

// *** 
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFFFCCE5,  // 1: Normal (very light pink)
		0xFFFF1493,  // 2: Hypoxic (deep pink/magenta)
		0xFF99004D,  // 3: Necrotic (very dark magenta)
		0xFFFFAC00,  // 4: Apoptotic (orange 0xFFFF8C00 - saddle brown 0xFF8B4513)
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark maroon)
	};
*/
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFD2B48C,  // 1: Normal tumor (light tan)
		0xFFB87333,  // 2: Hypoxic (copper)
		0xFF654321,  // 3: Necrotic (dark brown)
		0xFFB8F888,  // 4: Apoptotic ()
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};



	public static final double FONT_SCALE_FACTOR = 2.0;  // Adjust this to scale all fonts

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
			out.println("| Days | " + java.util.Arrays.toString(INJECTION_SCHEDULE) + " |");
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
			
//			System.out.println("Parameter report generated: " + filepath);
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
			out.println("injection_days," + java.util.Arrays.toString(INJECTION_SCHEDULE) + ",days,");
			out.println("dose_per_injection," + DOSE_PER_INJECTION * 1e9 + ",nmol,");
			out.println("hot_fraction," + HOT_FRACTION + ",,");
			
			// Grid
			out.println("grid_size," + GRID_SIZE + ",cells,");
			out.println("cell_length," + CELL_LENGTH * 1e6 + ",um,");
			out.println("domain_size," + DOMAIN_SIZE * 1e3 + ",mm,");
			out.println("initial_tumour_radius," + INITIAL_TUMOR_RADIUS * 1e3 + ",mm,");
			
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
			
//			System.out.println("Parameters exported to: " + filepath);
		}
	}


}
