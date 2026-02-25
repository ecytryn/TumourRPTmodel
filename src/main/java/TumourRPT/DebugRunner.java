package TumorRPT;

import HAL.Rand;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Debug Runner - Investigate specific parameter combinations from sweeps
 * 
 * Use this to understand WHY a specific point in your parameter sweep heatmap
 * succeeded or failed. Runs a single simulation with full visualization.
 * 
 * TWO WAYS TO USE:
 * 
 * 1. COMMAND LINE (Quick):
 *    ./gradlew runDebug --args="interval=30 skew=20"
 *    ./gradlew runDebug --args="dose=150 receptors=0.8"
 *    ./gradlew runDebug --args="interval=30 skew=20 seed=43"
 * 
 * 2. EDIT METHOD (Detailed):
 *    Edit setDebugParameters() method below
 *    ./gradlew runDebug
 */

public class DebugRunner {
    
	// Class variables
	private static double baselineReceptors = 5.0e-19;  // Default, can be loaded from file

    // Storage for doses when we need variable doses per injection
    private static double[] customDoses = null;
    private static int totalDays = 120;

	private static Integer customSeed = null;    
    private static boolean saveSFdata = true;

    // ===================================================================
    // METHOD 2: EDIT THESE PARAMETERS DIRECTLY
    // ===================================================================
    
    private static void setDebugParameters_IntervalSkew() {
        // INTERVAL-SKEW SWEEP DEBUG
        int interval = 30;          // Days between injections
        double skew = 20e-9;        // Dose skew in mol (20 nmol)
        
        // Calculate injection schedule
        int firstInjectionDay = 5;
        int secondInjectionDay = firstInjectionDay + interval;
        
        // Calculate doses
        double baseDose = 50e-9;   // mol
        double firstDose = baseDose + skew;
        double secondDose = baseDose - skew;
        
        // Use setExperiment to configure simulation
        String expName = String.format("Debug_I%d_S%.0f", interval, skew * 1e9);
        String expDesc = String.format("interval=%d, skew=%.0fnmol", interval, skew * 1e9);
        SimParams.setExperiment(expName, expDesc,
                               new int[]{firstInjectionDay, secondInjectionDay},
                               100e-6, 0.05, 0.02, SimParams.HYPOXIA_DEV_DAYS);
        
        // Store custom doses for injection
        customDoses = new double[]{firstDose, secondDose};
        totalDays = secondInjectionDay + 30;
        
        // Enable visualization
        SimParams.EXPORT_TUMOUR_OX_IMAGES = true;
        SimParams.EXPORT_OX_IMAGES = true;
        
        System.out.println("=== INTERVAL-SKEW DEBUG ===");
        System.out.println("Interval: " + interval + " days");
        System.out.println("Skew: " + (skew * 1e9) + " nmol");
        System.out.println("Injections: Day " + firstInjectionDay + " (" + (firstDose*1e9) + " nmol), " +
                          "Day " + secondInjectionDay + " (" + (secondDose*1e9) + " nmol)");
        System.out.println("Simulation: " + totalDays + " days");
        System.out.println("===========================\n");
    }
    
    private static void setDebugParameters_DoseReceptor() {
        // DOSE-RECEPTOR SWEEP DEBUG
        double totalDose = 150;     // nmol
        double receptorMultiplier = 0.8;
        
        double baselineReceptors = 6.64e-19;  // mol/cell
        double receptorDensity = receptorMultiplier * baselineReceptors;
        
        // Configure simulation
        int firstInjectionDay = 5;
        String expName = String.format("Debug_D%.0f_R%.2f", totalDose, receptorMultiplier);
        String expDesc = String.format("dose=%.0fnmol, receptors=%.0f%%", totalDose, receptorMultiplier*100);
        SimParams.setExperiment(expName, expDesc,
                               new int[]{firstInjectionDay},
                               250e-6,   // Match sweep: 250 μm radius
                               0.06, 0.019, SimParams.HYPOXIA_DEV_DAYS);
        
        // Override receptor density
        SimParams.RECEPTORS_PER_CELL_MOL = receptorDensity;
        
        // Store dose (convert nmol to mol)
        customDoses = new double[]{totalDose * 1e-9};
        totalDays = firstInjectionDay + 90;
        
        // Enable visualization  
        SimParams.EXPORT_TUMOUR_OX_IMAGES = true;
        SimParams.EXPORT_OX_IMAGES = true;
        
        System.out.println("=== DOSE-RECEPTOR DEBUG ===");
        System.out.println("Dose: " + totalDose + " nmol");
        System.out.println("Receptor density: " + (receptorMultiplier * 100) + "% of baseline");
        System.out.println("Simulation: " + totalDays + " days");
        System.out.println("===========================\n");
    }
    
    // ===================================================================
    // COMMAND LINE PARSING
    // ===================================================================
    
    private static void parseCommandLine(String[] args) {
        if (args.length == 0) {
            System.out.println("No command line args - using setDebugParameters() method");
            return;
        }
        
        String sweepType = null;
		String paramFile = null;
		Integer interval = null;
        Double skew = null;
        Double dose = null;
        Double receptors = null;
        Integer xIndex = null;
        Integer yIndex = null;
        
        // Parse key=value pairs
        for (String arg : args) {
            String[] parts = arg.split("=",2);
            if (parts.length != 2) continue;
            
            String key = parts[0].toLowerCase();
            String value = parts[1];
            
            switch (key) {
                case "sweep":
                    sweepType = value.toLowerCase();
                    break;
                case "paramfile":
                    paramFile = value;
                    break;
                case "interval":
                    interval = Integer.parseInt(value);
                    break;
                case "skew":
                    skew = Double.parseDouble(value) * 1e-9;  // Input in nmol, convert to mol
                    break;
                case "dose":
                    dose = Double.parseDouble(value);  // nmol
                    break;
                case "receptors":
                    receptors = Double.parseDouble(value);  // Multiplier of baseline
                    break;
				case "seed":
					customSeed = Integer.parseInt(value);
					break;
				case "x":
                    xIndex = Integer.parseInt(value);
                    break;
                case "y":
                    yIndex = Integer.parseInt(value);
                    break;
            }
        }

		// Load parameters if file specified
		if (paramFile != null && !paramFile.isEmpty()) {
			System.out.println("Loading configuration from: " + paramFile);
			Map<String, String> params = loadParametersFromCSV(paramFile);
			applyParametersFromFile(params);
		}
		
        // Handle sweep type with indices
        if (sweepType != null && xIndex != null && yIndex != null) {
            if (sweepType.contains("interval") || sweepType.contains("skew")) {
                // Copy arrays from IntervalSkewSweep - MUST MATCH YOUR ACTUAL SWEEP!
                int[] INTERVALS = {20, 28, 36, 44, 52, 60, 68, 76};
                double[] SKEWS = {-60e-9, -40e-9, -20e-9, 0, 20e-9, 40e-9, 60e-9};
                
                if (xIndex < INTERVALS.length && yIndex < SKEWS.length) {
                    interval = INTERVALS[xIndex];
                    skew = SKEWS[yIndex];
                    System.out.println("Using sweep indices: x=" + xIndex + ", y=" + yIndex);
                } else {
                    System.err.println("ERROR: Indices out of bounds!");
                    System.err.println("Valid x range: 0-" + (INTERVALS.length - 1));
                    System.err.println("Valid y range: 0-" + (SKEWS.length - 1));
                    System.exit(1);
                }
            } else if (sweepType.contains("dose") || sweepType.contains("receptor")) {
                // Copy arrays from DoseReceptorSweep - MUST MATCH YOUR ACTUAL SWEEP!
                double[] DOSES = {50, 75, 100, 125, 150, 175, 200, 225, 250};
                double BASELINE = 6.64e-19;
                double[] RECEPTOR_MULTS = {0.76, 0.84, 0.92, 1.0, 1.08, 1.16, 1.24, 1.32, 1.4, 1.48};
                
                if (xIndex < DOSES.length && yIndex < RECEPTOR_MULTS.length) {
                    dose = DOSES[xIndex];
                    receptors = RECEPTOR_MULTS[yIndex];
                    System.out.println("Using sweep indices: x=" + xIndex + ", y=" + yIndex);
                } else {
                    System.err.println("ERROR: Indices out of bounds!");
                    System.err.println("Valid x range: 0-" + (DOSES.length - 1));
                    System.err.println("Valid y range: 0-" + (RECEPTOR_MULTS.length - 1));
                    System.exit(1);
                }
            }
        }
        
        // Configure based on parsed parameters
        if (interval != null && skew != null) {
            configureIntervalSkew(interval, skew);
        } else if (dose != null && receptors != null) {
            configureDoseReceptor(dose, receptors);
        } else {
            System.err.println("ERROR: Incomplete parameters!");
            printUsage();
            System.exit(1);
        }
    }
    
    private static void configureIntervalSkew(int interval, double skew) {
        int firstInjectionDay = 5;
        int secondInjectionDay = firstInjectionDay + interval;
        double baseDose = 100e-9;
        double firstDose = baseDose + skew;
        double secondDose = baseDose - skew;
        
        String expName = String.format("Debug_I%d_S%.0f", interval, skew * 1e9);
        String expDesc = String.format("interval=%d, skew=%.0fnmol", interval, skew * 1e9);

		SimParams.EXPERIMENT_NAME = expName;
		SimParams.EXPERIMENT_DESCRIPTION = expDesc;
		SimParams.INJECTION_SCHEDULE = new int[]{firstInjectionDay, secondInjectionDay};
        
        customDoses = new double[]{firstDose, secondDose};
        totalDays = secondInjectionDay + 60;
        
        SimParams.EXPORT_TUMOUR_OX_IMAGES = true;
        SimParams.EXPORT_OX_IMAGES = true;
        
        System.out.println("=== INTERVAL-SKEW DEBUG (from command line) ===");
        System.out.println("Interval: " + interval + " days");
        System.out.println("Skew: " + (skew * 1e9) + " nmol");
        System.out.println("Injections: Day " + firstInjectionDay + " (" + (firstDose*1e9) + " nmol), " +
                          "Day " + secondInjectionDay + " (" + (secondDose*1e9) + " nmol)");
        System.out.println("===========================\n");
    }
    
    private static void configureDoseReceptor(double dose, double receptorMult) {
//        double baselineReceptors = 6.64e-19;  // <-- this is now defined by reading parameters.csv
        double receptorDensity = receptorMult * baselineReceptors;
        int firstInjectionDay = 5;
        
        String expName = String.format("Debug_D%.0f_R%.2f", dose, receptorMult);
        String expDesc = String.format("dose=%.0fnmol, receptors=%.0f%%", dose, receptorMult*100);

		SimParams.EXPERIMENT_NAME = expName;
		SimParams.EXPERIMENT_DESCRIPTION = expDesc;
		SimParams.INJECTION_SCHEDULE = new int[]{firstInjectionDay};
            
        SimParams.RECEPTORS_PER_CELL_MOL = receptorDensity;
        customDoses = new double[]{dose * 1e-9};
        totalDays = firstInjectionDay + 90;
        
        SimParams.EXPORT_TUMOUR_OX_IMAGES = true;
        SimParams.EXPORT_OX_IMAGES = true;
        
        System.out.println("=== DOSE-RECEPTOR DEBUG (from command line) ===");
        System.out.println("Dose: " + dose + " nmol");
        System.out.println("Receptor density: " + (receptorMult * 100) + "% of baseline");
        System.out.println("===========================\n");
    }
    
	/**
	 * Load parameters from a parameters.csv file
	 * Returns map of parameter_name -> value
	 */
	private static Map<String, String> loadParametersFromCSV(String filepath) {
		Map<String, String> params = new HashMap<>();
		
		try (BufferedReader br = new BufferedReader(new FileReader(filepath))) {
			String line;
			boolean firstLine = true;
			
			while ((line = br.readLine()) != null) {
				if (firstLine) {
					firstLine = false;  // Skip header
					continue;
				}
				
				// Handle quoted values (arrays like "[1, 2]")
				// Split on commas NOT inside quotes
				String[] parts = line.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
				
				if (parts.length >= 2) {
					String param = parts[0].trim();
					String value = parts[1].trim();
					
					// Remove quotes if present
					if (value.startsWith("\"") && value.endsWith("\"")) {
						value = value.substring(1, value.length() - 1);
					}
					
					params.put(param, value);
				}
			}
		} catch (IOException e) {
			System.err.println("ERROR loading parameters: " + e.getMessage());
		}
		
		return params;
	}
	
	/**
	 * Apply loaded parameters to SimParams
	 * Only applies FIXED parameters (not sweep-specific ones)
	 */
	private static void applyParametersFromFile(Map<String, String> params) {
		System.out.println("\nLoading fixed parameters from file:");
		
		if (params.containsKey("hypoxia_dev_days")) {
			SimParams.HYPOXIA_DEV_DAYS = Integer.parseInt(params.get("hypoxia_dev_days"));
			System.out.println("  hypoxia_dev_days = " + SimParams.HYPOXIA_DEV_DAYS);
		}
		
		if (params.containsKey("initial_tumor_radius")) {
			double radius_um = Double.parseDouble(params.get("initial_tumor_radius"));
			SimParams.INITIAL_TUMOR_RADIUS = radius_um * 1e-6;  // um -> m
			System.out.println("  initial_tumor_radius = " + radius_um + " um");
		}
		
		if (params.containsKey("vessel_density_config")) {
			SimParams.VESSEL_DENSITY_CONFIG = params.get("vessel_density_config");
			System.out.println("  vessel_density_config = " + SimParams.VESSEL_DENSITY_CONFIG);
		}
		
		if (params.containsKey("alpha_normoxic")) {
			SimParams.ALPHA_NORMAL = Double.parseDouble(params.get("alpha_normoxic"));
			System.out.println("  alpha_normoxic = " + SimParams.ALPHA_NORMAL);
		}
		
		if (params.containsKey("beta_normoxic")) {
			SimParams.BETA_NORMAL = Double.parseDouble(params.get("beta_normoxic"));
			System.out.println("  beta_normoxic = " + SimParams.BETA_NORMAL);
		}
		
		if (params.containsKey("alpha_hypoxic")) {
			SimParams.ALPHA_HYPOXIC = Double.parseDouble(params.get("alpha_hypoxic"));
			System.out.println("  alpha_hypoxic = " + SimParams.ALPHA_HYPOXIC);
		}
		
		if (params.containsKey("beta_hypoxic")) {
			SimParams.BETA_HYPOXIC = Double.parseDouble(params.get("beta_hypoxic"));
			System.out.println("  beta_hypoxic = " + SimParams.BETA_HYPOXIC);
		}
		
		if (params.containsKey("hot_fraction")) {
			SimParams.HOT_FRACTION = Double.parseDouble(params.get("hot_fraction"));
			System.out.println("  hot_fraction = " + SimParams.HOT_FRACTION);
		}
		
		// Get baseline receptors for calculating multipliers
		if (params.containsKey("receptors_baseline_mol")) {
			baselineReceptors = Double.parseDouble(params.get("receptors_baseline_mol"));
			System.out.println("  receptors_baseline = " + baselineReceptors + " mol/cell");
		}
		
		System.out.println();
	}    
    
    private static void printUsage() {
        System.out.println("\n=== DEBUG RUNNER USAGE ===\n");
        System.out.println("Method 1 - Direct parameters:");
        System.out.println("  ./gradlew runDebug --args=\"interval=30 skew=20\"");
        System.out.println("  ./gradlew runDebug --args=\"dose=150 receptors=0.8\"");
        System.out.println();
        System.out.println("Method 2 - From heatmap indices:");
        System.out.println("  ./gradlew runDebug --args=\"sweep=interval x=4 y=6\"");
        System.out.println("  ./gradlew runDebug --args=\"sweep=dose x=3 y=2\"");
        System.out.println();
        System.out.println("Method 3 - Edit setDebugParameters() method and run:");
        System.out.println("  ./gradlew runDebug");
        System.out.println();
        System.out.println("Parameters:");
        System.out.println("  interval   - Days between injections (e.g., 30)");
        System.out.println("  skew       - Dose skew in nmol (e.g., 20 for +20/-20)");
        System.out.println("  dose       - Total dose in nmol (e.g., 150)");
        System.out.println("  receptors  - Receptor multiplier (e.g., 0.8 for 80% of baseline)");
        System.out.println("  sweep      - Sweep type: 'interval' or 'dose'");
        System.out.println("  x, y       - Array indices from your sweep parameters");
        System.out.println();
    }
    
    // ===================================================================
    // MAIN - Run the debug simulation
    // ===================================================================
    
    public static void main(String[] args) throws IOException {
        System.out.println("\n" + "=".repeat(60));
        System.out.println("        TUMOR RPT DEBUG RUNNER");
        System.out.println("=".repeat(60) + "\n");
                
        // Parse command line OR use edited method
        if (args.length > 0) {
            parseCommandLine(args);
        } else {
            System.out.println("Using parameters from setDebugParameters() method\n");
            
            // CHOOSE WHICH DEBUG MODE:
            setDebugParameters_IntervalSkew();  // For interval-skew sweeps
            // setDebugParameters_DoseReceptor();   // For dose-receptor sweeps
        }

        // Initialize random seed
		int seed = (customSeed != null) ? customSeed : 42;
		Rand rng = new Rand(seed);
		System.out.println("Using random seed: " + seed);
        
        // Create output directory
        String timestamp = java.time.LocalDateTime.now()
                .format(java.time.format.DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"));
        String outputDir = "results/debug_runs/" + SimParams.EXPERIMENT_NAME + "_" + timestamp;
        new File(outputDir).mkdirs();

		saveGitInfo(outputDir);
        
        // Set up output directories (matches Main.java pattern)
        SimParams.OUTPUT_DIR_BASE = outputDir;
        SimParams.OUTPUT_DIR_TUMOUR_IMAGES = outputDir + "/tumour_images";
        SimParams.OUTPUT_DIR_OXYGEN_IMAGES = outputDir + "/oxygen_images";
        SimParams.OUTPUT_DIR_SF_IMAGES = outputDir + "/sf_images";
        
        new File(SimParams.OUTPUT_DIR_TUMOUR_IMAGES).mkdirs();
        new File(SimParams.OUTPUT_DIR_OXYGEN_IMAGES).mkdirs();
        new File(SimParams.OUTPUT_DIR_SF_IMAGES).mkdirs();
        
        System.out.println("Output directory: " + outputDir + "\n");
        
        // Run the simulation
        runDebugSimulation(outputDir, rng);
        
        System.out.println("\n" + "=".repeat(60));
        System.out.println("DEBUG SIMULATION COMPLETE");
        System.out.println("=".repeat(60));
        System.out.println("\nCheck outputs in: " + outputDir);
        System.out.println("Review tumour_images/ to understand the dynamics\n");
    }
    
    private static void runDebugSimulation(String outputDir, Rand rng) throws IOException {
        // Initialize grid (matches Main.java pattern)
        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);
        DataLogger logger = new DataLogger();
        DaVinci drawer = new DaVinci(model);
        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};
        
        // Initialize simulation
        int dayCount = -1;
        model.Init(dayCount, drawer, logger);

		// Develop hypoxia without growth
		if (SimParams.HYPOXIA_DEV_DAYS > 0) {
			model.developHypoxiaWithoutGrowth(SimParams.HYPOXIA_DEV_DAYS, false);
		}
        
        int initialPop = model.countTumorCells();
        System.out.println("Initial tumor cells: " + initialPop);
        System.out.println("Starting simulation for " + totalDays + " days...\n");
        
        // SF logging: track SF for day-0 cohort and previous-day cohort each day
        // Columns: day, SF_cohort0_norm, SF_cohort0_hypo, SF_prevDay_norm, SF_prevDay_hypo
        ArrayList<double[]> sfData = saveSFdata ? new ArrayList<>() : null;
        
        // Main simulation loop
        for (dayCount = 0; dayCount <= totalDays; dayCount++) {
            int day = dayCount;
            
            // Check for injections
            for (int i = 0; i < SimParams.INJECTION_SCHEDULE.length; i++) {
                if (SimParams.INJECTION_SCHEDULE[i] == dayCount) {
                    // Use custom dose if specified, otherwise use default
                    double injectionDose = (customDoses != null && i < customDoses.length) ? 
                                          customDoses[i] : SimParams.DOSE_PER_INJECTION;
                    
                    double hotDose = injectionDose * SimParams.HOT_FRACTION;
                    double coldDose = injectionDose * (1.0 - SimParams.HOT_FRACTION);
                    
                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    currentPK[0] += hotDose;   // N_cen_hot
                    currentPK[1] += coldDose;  // N_cen_cold
                    model.PKStateVariables.set(model.PKStateVariables.size() - 1, currentPK);
                    
                    System.out.printf("Day %d: Injected %.1f nmol (%.1f%% hot)%n",
                                    dayCount, injectionDose * 1e9, SimParams.HOT_FRACTION * 100);
                }
            }
            
            // Export images every 5 days
            if (SimParams.EXPORT_TUMOUR_OX_IMAGES && day % 5 == 0) {
                drawer.gridDraw(visualizationMaskList);
                double[] valueList = MyUtils.lastElementOfDoubleArrayList(model.DoseRateList);
                drawer.plot(dayCount, valueList[0]);
                
                String imageFile = String.format("%s/day_%03d.png", 
                                                SimParams.OUTPUT_DIR_TUMOUR_IMAGES, dayCount);
                logger.saveFigureTotal(imageFile, drawer, dayCount, false, false);
            }
            
            // Step simulation
            model.Step(dayCount);
            
            // Log SF data for key cohorts
            // Columns: day, SF_c0_norm, SF_c0_hypo, SF_prev_norm, SF_prev_hypo,
            //          D_c0, A_c0, Gnum_c0, D_prev, A_prev, Gnum_prev
            if (saveSFdata) {
                double sf_c0_norm = model.radioBio.calculateSF(0, SimParams.NORMAL);
                double sf_c0_hypo = model.radioBio.calculateSF(0, SimParams.HYPOXIC);
                // Previous-day cohort: born 24 hours ago. Guard against negative on day 0.
                int prevDayCohort = Math.max(0, (dayCount - 1) * 24);
                double sf_prev_norm = model.radioBio.calculateSF(prevDayCohort, SimParams.NORMAL);
                double sf_prev_hypo = model.radioBio.calculateSF(prevDayCohort, SimParams.HYPOXIC);

                double[] state_c0   = model.radioBio.getCohortState(0);
                double[] state_prev = model.radioBio.getCohortState(prevDayCohort);

                sfData.add(new double[] {
                    dayCount,
                    sf_c0_norm, sf_c0_hypo, sf_prev_norm, sf_prev_hypo,
                    state_c0[0],   state_c0[1],   state_c0[2],    // D, A, G_num for cohort 0
                    state_prev[0], state_prev[1], state_prev[2]   // D, A, G_num for prev-day cohort
                });
            }
            
            // Report progress every 10 days
            if (dayCount % 10 == 0 || dayCount == totalDays) {
                int currentPop = model.countTumorCells();
                System.out.printf("Day %d: %d cells%n", dayCount, currentPop);
            }
        }
        
        // Final report
        System.out.println("\n=== FINAL RESULTS ===");
        int[] typeCounts = new int[SimParams.NUM_CELL_TYPES];
        for (Cell cell : model) {
            if (cell != null && cell.type != SimParams.VESSEL) {
                typeCounts[cell.type]++;
            }
        }
        
        System.out.println("Final cell counts:");
        System.out.println("  Normal: " + typeCounts[SimParams.NORMAL]);
        System.out.println("  Hypoxic: " + typeCounts[SimParams.HYPOXIC]);
        System.out.println("  Necrotic: " + typeCounts[SimParams.NECROTIC]);
        System.out.println("  Apoptotic: " + typeCounts[SimParams.APOPTOTIC]);
        System.out.println("  TOTAL: " + model.countTumorCells());
        
        boolean cured = model.countTumorCells() < 10;
        System.out.println("\nOutcome: " + (cured ? "CURE" : "FAILURE"));
        
        // Save data files
        String popcsvFileName = outputDir + "/populations.csv";
        String dosecsvFileName = outputDir + "/dose.csv";
        String PKvarscsvFileName = outputDir + "/pkStateVariables.csv";
        
        logger.log(model.PopsOverTime, popcsvFileName);
        logger.log(model.DoseRateList, dosecsvFileName);
        logger.log(model.PKStateVariables, PKvarscsvFileName);
        
        // Save SF data if enabled
        if (saveSFdata && sfData != null) {
            String sfCsvFileName = outputDir + "/sfData.csv";
            logger.log(sfData, sfCsvFileName);
            System.out.println("SF data saved to: " + sfCsvFileName);
        }
        
        // Generate parameter report
        SimParams.generateParameterReport(outputDir + "/parameters.md");
		SimParams.exportParametersToCSV(outputDir + "/parameters.csv",
                                SimParams.INJECTION_SCHEDULE,
                                customDoses != null ? customDoses : getDefaultDoses(),
                                SimParams.HOT_FRACTION, SimParams.RECEPTORS_PER_CELL_MOL);
    }

	private static double[] getDefaultDoses() {
		double[] doses = new double[SimParams.INJECTION_SCHEDULE.length];
		for (int i = 0; i < doses.length; i++) {
			doses[i] = SimParams.DOSE_PER_INJECTION;
		}
		return doses;
	}

	private static void saveGitInfo(String outputDir) {
		try {
			Process process = Runtime.getRuntime().exec("git rev-parse HEAD");
			BufferedReader reader = new BufferedReader(
				new InputStreamReader(process.getInputStream()));
			String commit = reader.readLine();
			
			try (PrintWriter out = new PrintWriter(
					new FileWriter(outputDir + "/git_info.txt"))) {
				out.println("Git commit: " + commit);
				out.println("Date: " + java.time.LocalDateTime.now());
			}
		} catch (Exception e) {
			// Git not available or not a git repo - silently skip
		}
	}


}
