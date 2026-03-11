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
    // Per-experiment settings: vary these via setupExperiment() in Main.java
    // or directly at the top of each sweep's main() method.
    // =================================================================

	public static String EXPERIMENT_NAME = "default";
	public static String EXPERIMENT_DESCRIPTION = "Default run";
	public static int[] INJECTION_SCHEDULE = {45};
	public static double INITIAL_TUMOR_RADIUS = 333e-6;

    // Vessel density configuration
    // Use scripts/GenerateVessels/GenerateUniformVessels.py to create new csv files.
    // R_REPEL values refer to the repel radius parameter used in GenerateUniformVessels.py
    // Density of vessels is "inversely" proportional to R_REPEL
    // I switched to denoting csv files by the density of capillaries reported as output from
    // GenerateUniformVessels.py. The csv files are in src/main/resources/vasculature/.
	// R_REPEL = 20 --> Density = 625 vessels/mm^2 
	// R_REPEL = 40 --> Density = 616 vessels/mm^2 
	// R_REPEL = 50 --> Density = 605 vessels/mm^2 
	// R_REPEL = 60 --> Density = 592 vessels/mm^2 
	// R_REPEL = 80 --> Density = 566 vessels/mm^2 
	//                  Density = 374 vessels/mm^2 
	//                  Density = 291 vessels/mm^2 
	public static String VESSEL_DENSITY_CONFIG = "605";  // Default configuration

	// Setter method:
	public static void setExperiment(String name, String desc, int[] injections,
									  double radius, int hypoxiaDevDays) {
		EXPERIMENT_NAME = name;
		EXPERIMENT_DESCRIPTION = desc;
		INJECTION_SCHEDULE = injections;
		INITIAL_TUMOR_RADIUS = radius;
		HYPOXIA_DEV_DAYS = hypoxiaDevDays;
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
	public static boolean FREEZE_TUMOR = false;          // Set true to test PK without cell dynamics    

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
	// Originally used K_ON=0.046 * ... and K_OFF=0.386 but updated based on literature search.
	// The change should not make a noticeable difference because these are still in the QSS regime.
    public static final double K_ON = 1.5e-3 * 1e6 / 60.0;     // m^3/(mol s) (binding) (1.5e-3 is in lit/nmol/min)
    public static final double K_OFF = 1.2e-2 / 60.0;          // 1/s (unbinding) (1.2e-2 is in 1/min)
    public static final double K_INT = 0.001 / 60.0;          // 1/s (internalization)
    public static final double K_REL = 2e-4 / 60.0;           // 1/s (release from cells)
    
    // =================================================================
    // PHARMACOKINETICS - Volumes
    // =================================================================
    public static final double V_CENTRAL = 0.5e-3;          // m^3 (central/arterial, 0.458 L)
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
    public static final double RECEPTORS_PER_CELL = 3e5;      // receptors/cell
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
	
	// Number of days to run the pre-simulation to allow hypoxia to develop in the initial tumour state.
	public static int HYPOXIA_DEV_DAYS = 40;  // Default value

    // =================================================================
    // RADIOBIOLOGY
    // =================================================================

	public static double ALPHA_NORMAL = 0.15;           // Gy^(-1)
    public static double BETA_NORMAL = 0.05;           // Gy^(-2)
    public static final double ALPHA_HYPOXIC = 0.1;          // Gy^(-1) (0.15/2.5, OER=1.5 - consistent with clinic context)
    public static final double BETA_HYPOXIC = 0.02;          // Gy^(-2) (0.048/2.5, OER=1.5)


    public static final double REPAIR_RATE = 0.7 / 3600.0;   // 1/s (0.7 per hour)    

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
	
//	public static final int maxLookupAge = 200;  // days
    
    // =================================================================
    // INJECTION PROTOCOL
    // =================================================================
    // **FOR SINGLE RUNS:** Change these parameters to run individual experiments
    // **FOR PARAMETER SWEEPS:** Use IntervalSkewSweep.java which programmatically 
    //                          varies inter-injection interval and dose distribution
    //                           or DoseReceptorSweep.java (injected amount and receptor density)
    // Main.java uses these values for single-run mode
    
    public static final double DOSE_PER_INJECTION = 100e-9;   // mol (100 nmol per injection is baseline)
    public static double HOT_FRACTION = 0.1;            // Fraction that is radioactive (0.1 = 10%)
/**	Moved to Main
    public static final int[] INJECTION_SCHEDULE = {5};  // Days to inject (e.g., {5, 35, 65, 95})
*/    
    // Examples for different treatment schedules:
    // Single dose:      {35}
    // 2 doses:          {5, 35}
    // 4 doses:          {5, 35, 65, 95}
    // Weekly x 4:       {7, 14, 21, 28}
    
    // Simulation length (days to simulate after last injection)
    public static final int DAYS_AFTER_LAST_INJECTION = 40;
    
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

	// A colour-blind friendly upgrade:
	public static final int[] COLORLIST = {
		0x88888888,  // 0: Healthy (gray)
		0xFFD2B48C,  // 1: Normal tumor (light tan)
		0xFFB87333,  // 2: Hypoxic (copper)
		0xFF654321,  // 3: Necrotic (dark brown)
		0xFFB8F888,  // 4: Apoptotic ()
		0xFFFF0000,  // 5: Open vessel (bright red)
		0xFF8B0000   // 6: Occluded vessel (dark red/maroon)
	};

	// =================================================================
	// ZOOMED IMAGE EXPORT PARAMETERS
	// Edit these to tune the zoomed figure images.
	// ZOOM_HALF_WIDTH: half-side of the crop window in cells (centred on tumour)
	//   e.g. 50 => 100x100 cell crop = 1mm x 1mm at 10um/cell
	// ZOOM_PIXEL_SCALE: pixels per cell in the output image (upscaling factor)
	//   e.g. 4 => each cell becomes a 4x4 pixel block
	// =================================================================
	public static final int ZOOM_HALF_WIDTH   = 50;   // cells (50 --> 500 um half-width)
	public static final int ZOOM_PIXEL_SCALE  = 4;    // pixels per cell
	public static final int ZOOM_FONT_SIZE       = 40;   // pt — day label and scale bar text
	public static final int ZOOM_SCALE_BAR_HEIGHT = 8;   // pixels — height of scale bar rectangle
	

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
 * Export all simulation parameters to CSV
 * 
 * @param filepath Output file path
 * @param injectionTimes Days when injections occur (sweep-specific)
 * @param injectionDoses Amounts injected in mol (sweep-specific)
 * @param hotFraction Fraction that is radioactive (could vary)
 * @param receptorDensity Receptors per cell in mol (sweep-specific for DoseReceptorSweep)
 */
public static void exportParametersToCSV(String filepath, 
                                         int[] injectionTimes, 
                                         double[] injectionDoses,
                                         double hotFraction,
                                         double receptorDensity) throws IOException {
    try (PrintWriter out = new PrintWriter(new FileWriter(filepath))) {
        out.println("parameter,value,units,description");
        
        // ===== EXPERIMENT METADATA =====
        out.println("experiment_name," + EXPERIMENT_NAME + ",,Experiment identifier");
        out.println("experiment_description," + EXPERIMENT_DESCRIPTION.replace(",", ";") + ",,");
        out.println("experiment_date," + java.time.LocalDateTime.now() + ",,Run timestamp");

        // ===== HYPOXIA DEVELOPMENT =====
        out.println("hypoxia_dev_days," + HYPOXIA_DEV_DAYS + ",days,Days to develop hypoxia before treatment");
        
        // ===== INJECTION PROTOCOL (SWEEP-SPECIFIC) =====
        out.println("num_injections," + injectionTimes.length + ",,");
        
        // Format arrays as quoted strings
        String injTimesStr = arrayToString(injectionTimes);
        String injDosesStr = arrayToString(injectionDoses);
        
        out.println("injection_days,\"" + injTimesStr + "\",days,");
        out.println("injection_doses,\"" + injDosesStr + "\",nmol,");
        
        // Total dose
        double totalDose = 0;
        for (double d : injectionDoses) totalDose += d;
        out.println("total_injected_amount," + (totalDose * 1e9) + ",nmol,");

        out.println("hot_fraction," + hotFraction + ",,Fraction that is radioactive");
        
        // ===== TUMOR INITIAL CONDITIONS =====
        out.println("initial_tumor_radius," + (INITIAL_TUMOR_RADIUS * 1e6) + ",um,");
        out.println("vessel_density_config," + VESSEL_DENSITY_CONFIG + ",,Vessel configuration file ID");
        
        // ===== GRID/DOMAIN =====
        out.println("grid_size," + GRID_SIZE + ",cells,");
        out.println("cell_length," + (CELL_LENGTH * 1e6) + ",um,");
        out.println("domain_size," + (DOMAIN_SIZE * 1e3) + ",mm,");
        out.println("time_step," + (TIME_STEP / 3600) + ",hours,");
        out.println("cell_cycle," + (CELL_CYCLE / 3600) + ",hours,");
        
        // ===== PHARMACOKINETICS =====
        out.println("lambda_bio," + (LAMBDA_BIO * 3600) + ",1/hr,Biological clearance");
        out.println("lambda_decay," + (LAMBDA_DECAY * 3600) + ",1/hr,Lu-177 radioactive decay");
        out.println("k_on," + (K_ON / 1e6 * 60) + ",L/(nmol*min),Association rate");
        out.println("k_off," + (K_OFF * 60) + ",1/min,Dissociation rate");
        out.println("k_int," + (K_INT * 60) + ",1/min,Internalization rate");
        out.println("k_rel," + (K_REL * 60) + ",1/min,Release rate");
        out.println("v_central," + (V_CENTRAL * 1e3) + ",L,Central compartment volume");
        out.println("interstitial_fraction," + INTERSTITIAL_FRACTION + ",,Fraction of tumor that is extracellular");
        
        // ===== RADIOBIOLOGY =====
        out.println("alpha_normoxic," + ALPHA_NORMAL + ",Gy^-1,Linear term normoxic");
        out.println("beta_normoxic," + BETA_NORMAL + ",Gy^-2,Quadratic term normoxic");
        out.println("alpha_hypoxic," + ALPHA_HYPOXIC + ",Gy^-1,Linear term hypoxic");
        out.println("beta_hypoxic," + BETA_HYPOXIC + ",Gy^-2,Quadratic term hypoxic");
        out.println("repair_rate," + (REPAIR_RATE * 3600) + ",1/hr,DNA repair rate (mu)");
        out.println("e_beta_lu177," + E_BETA_LU177 + ",J,Average beta particle energy");
        
        // ===== OXYGEN =====
        out.println("p_o2_vessel," + (P_O2_VESSEL / MMHG_TO_PA) + ",mmHg,Capillary oxygen pressure");
        out.println("p_o2_hypoxic," + (P_O2_HYPOXIC / MMHG_TO_PA) + ",mmHg,Hypoxia threshold");
        out.println("p_o2_necrotic," + (P_O2_NECROTIC / MMHG_TO_PA) + ",mmHg,Necrosis threshold");
        out.println("d_o2," + (D_O2 * 1e4) + ",cm^2/s,Oxygen diffusion coefficient");
        
        out.println("consumption_healthy," + CONSUMPTION_HEALTHY + ",Pa/s,");
        out.println("consumption_normal," + CONSUMPTION_NORMAL + ",Pa/s,");
        out.println("consumption_hypoxic," + CONSUMPTION_HYPOXIC + ",Pa/s,");
        out.println("consumption_necrotic," + CONSUMPTION_NECROTIC + ",Pa/s,");
        
        // ===== CELL BIOLOGY (SWEEP-SPECIFIC) =====
        // Use passed receptor density (varies in DoseReceptorSweep)
        double receptorsPerCell = receptorDensity * AVOGADRO;
        out.println("receptors_per_cell," + String.format("%.0f", receptorsPerCell) + ",,Number per cell");
        out.println("receptors_per_cell_mol," + receptorDensity + ",mol/cell,ACTUAL USED IN THIS RUN");
        
        // Also export baseline for reference
        out.println("receptors_baseline_mol," + RECEPTORS_PER_CELL_MOL + ",mol/cell,Default baseline value");
        
        out.println("vessel_influence_radius," + VESSEL_INFLUENCE_RADIUS + ",cells,For vessel occlusion");
        
        // ===== COMPUTATIONAL SETTINGS =====
        out.println("export_tumour_ox_images," + EXPORT_TUMOUR_OX_IMAGES + ",,");
        out.println("export_ox_images," + EXPORT_OX_IMAGES + ",,");
        out.println("export_sf_images," + EXPORT_SF_IMAGES + ",,");
        out.println("freeze_tumor," + FREEZE_TUMOR + ",,For PK testing");
        out.println("enable_pbpk_logging," + ENABLE_PBPK_LOGGING + ",,");
    }
}

/**
 * Helper: Convert int array to string "[1, 2, 3]"
 */
private static String arrayToString(int[] array) {
    StringBuilder sb = new StringBuilder("[");
    for (int i = 0; i < array.length; i++) {
        sb.append(array[i]);
        if (i < array.length - 1) sb.append(", ");
    }
    sb.append("]");
    return sb.toString();
}

/**
 * Helper: Convert double array to string (in nmol) "[45.0, 55.0]"
 */
private static String arrayToString(double[] array) {
    StringBuilder sb = new StringBuilder("[");
    for (int i = 0; i < array.length; i++) {
        sb.append(String.format("%.1f", array[i] * 1e9));  // Convert mol to nmol
        if (i < array.length - 1) sb.append(", ");
    }
    sb.append("]");
    return sb.toString();
}

/**
 * Backward-compatible wrapper that uses SimParams defaults
 * Use this only for single runs where SimParams values are actually used
 */
public static void exportParametersToCSV(String filepath) throws IOException {
    // Create doses array from DOSE_PER_INJECTION
    double[] doses = new double[INJECTION_SCHEDULE.length];
    for (int i = 0; i < doses.length; i++) {
        doses[i] = DOSE_PER_INJECTION;
    }
    
    exportParametersToCSV(filepath, INJECTION_SCHEDULE, doses, HOT_FRACTION,
                         RECEPTORS_PER_CELL_MOL);
}


}
