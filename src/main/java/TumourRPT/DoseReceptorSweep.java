package TumorRPT;

import HAL.Rand;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.stream.IntStream;

/**
 * Parameter sweep over dose and receptor density
 * 
 * Explores 2D parameter space:
 * - Axis 1: Total dose (nmol) - split equally between two injections
 * - Axis 2: Receptor density (receptors per cell)
 * 
 * Research question: How does required dose to eliminate tumor vary with 
 * patient-specific receptor expression?
 * 
 * Stage 1: Coarse grid to map out dose-response landscape
 */

public class DoseReceptorSweep {
    
    // =======================================================================
    // SWEEP PARAMETERS - Edit these to define the parameter space
    // =======================================================================
    
	private static final double[] DOSES = 
		IntStream.range(0, 13)  // 13 points
				 .mapToDouble(i -> (50 + i * 12.5))  // 50 --> 200 nmol
				 .toArray();
    
    // RECEPTOR_DENSITIES: Receptors per cell (in moles)
    // Expressed as multiples of SimParams.RECEPTORS_PER_CELL_MOL (the canonical baseline).

	private static final double[] RECEPTOR_DENSITIES = 
		IntStream.range(3, 6)  // 3 points
				 .mapToDouble(i -> (1.28 + i * 0.04) * SimParams.RECEPTORS_PER_CELL_MOL)
				 .toArray();
//		IntStream.range(0, 16)  // 16 points
//				 .mapToDouble(i -> (0.68 + i * 0.04) * SimParams.RECEPTORS_PER_CELL_MOL)
//				 .toArray();
                 
    // Output suffix - change this when refining to avoid overwriting
    private static final String OUTPUT_SUFFIX = "";
    
    // Number of replicates per parameter combination
    private static final int NUM_REPLICATES = 20;
    
    // =======================================================================
    // FIXED PARAMETERS - Constant across all sweep runs
    // =======================================================================
    
    private static final int NUM_INJECTIONS = 1;
    private static final double HOT_FRACTION = 0.1;
    private static final int FIRST_INJECTION_DAY = 5;
//    private static final int INJECTION_INTERVAL = 30;  // Fixed 30-day interval
    private static final int MIN_DAYS_AFTER_LAST_INJECTION = 60;
    
    // Output configuration
    private static final String BASE_OUTPUT_DIR = "results/DoseReceptorSweep/DoseReceptorSweep" + OUTPUT_SUFFIX;
    
    // Image export control (set to false for faster sweeps)
    private static final boolean EXPORT_IMAGES = false;
    
    public static void main(String[] args) throws IOException {

        // ===================================================================
        // RUN CONFIGURATION - edit here to switch between high/low cap density
        // All other parameters come from SimParams canonical values.
        // ===================================================================
        SimParams.INITIAL_TUMOR_RADIUS = 21e-6;    // 10 um: small tumour a bit below the cure threshold
        SimParams.HYPOXIA_DEV_DAYS = 40;            // Pre-sim steps to establish hypoxia
        SimParams.VESSEL_DENSITY_CONFIG = "605";    // "605" = high density; "374" = low density
        // ===================================================================

        String timestamp = LocalDateTime.now().format(
            DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss"));
        
        String sweepDir = BASE_OUTPUT_DIR + "_" + timestamp;
        new File(sweepDir).mkdirs();
        
        String summaryPath = sweepDir + "/sweep_summary.csv";
        
        // Initialize CSV output
        PrintWriter csvWriter;
        try {
            csvWriter = new PrintWriter(new FileWriter(summaryPath));
        } catch (IOException e) {
            throw new RuntimeException("Could not open summary CSV at: " + summaryPath, e);
        }

        // CSV header
        csvWriter.println("dose_nmol,receptors_per_cell_mol,dose1_nmol,replicate,finalTumorCount,outcome");
        
        int totalCombinations = DOSES.length * RECEPTOR_DENSITIES.length;
        int totalRuns = totalCombinations * NUM_REPLICATES;
        int currentRun = 0;
        
        System.out.println("=== DOSE-RECEPTOR SWEEP START ===");
        System.out.println("Sweep timestamp: " + timestamp);
        System.out.println("Parameter combinations: " + totalCombinations);
        System.out.println("Replicates per combination: " + NUM_REPLICATES);
        System.out.println("Total runs: " + totalRuns);
        System.out.println("Output directory: " + sweepDir);
        System.out.println();
        
        // Loop over parameter space
        for (double dose : DOSES) {
            for (double receptorDensity : RECEPTOR_DENSITIES) {
                
                // Generate injection schedule
                int[] injectionTimes = generateInjectionTimes();
                double[] injectionDoses = generateInjectionDoses(dose); // in mol
                
                // Calculate simulation days
                int lastInjectionDay = injectionTimes[NUM_INJECTIONS - 1];
                int simulationDays = lastInjectionDay + MIN_DAYS_AFTER_LAST_INJECTION;
                
                // Run replicates for this parameter combination
                for (int rep = 0; rep < NUM_REPLICATES; rep++) {
                    currentRun++;
                    
                    // Create output directory for this run
                    String runDir = String.format("%s/dose_%.0f_recep_%.2e_rep_%d", 
                                                 sweepDir, dose, receptorDensity, rep + 1);
                    new File(runDir).mkdirs();

					saveGitInfo(runDir);
                    
                    System.out.printf("[%d/%d] Running: dose=%.0f nmol, receptors=%.2e mol/cell, replicate=%d/%d%n", 
                                     currentRun, totalRuns, dose, receptorDensity, rep + 1, NUM_REPLICATES);
//                    System.out.printf("  Injection times: %d, %d days%n",
//                                     injectionTimes[0], injectionTimes[1]);
//                    System.out.printf("  Doses: %.1f, %.1f nmol%n",
//                                     injectionDoses[0] * 1e9, injectionDoses[1] * 1e9);
                    System.out.printf("  Simulation length: %d days%n", simulationDays);

					// Use replicate-based seed for reproducibility
//					int seed = 42 + rep;  // Seeds: 42, 43, 44, 45, 46 for reps 1-5
					int seed = rep * 999983 + currentRun;
					Rand rng = new Rand(seed);
                    
                    // Run simulation with this receptor density
                    SimResult result = runSingleSimulation(rng, injectionTimes, injectionDoses, 
                                                          receptorDensity, runDir, simulationDays);

                    // Classify outcome
                    String outcome = (result.finalTumorCount == 0) ? "CURE" : "FAILURE";                

                    System.out.printf("  Result: %d tumor cells remaining -> %s%n%n", 
                                     result.finalTumorCount, outcome);

                    // Write to CSV
                    csvWriter.printf("%.1f,%.3e,%.1f,%d,%d,%s%n",
                           dose, receptorDensity,
                           injectionDoses[0] * 1e9,
                           rep + 1, result.finalTumorCount, outcome);
                    csvWriter.flush(); // Flush after each run
                }
            }
        }
        
        csvWriter.close();
        
        System.out.println("=== DOSE-RECEPTOR SWEEP COMPLETE ===");
        System.out.println("Results saved to: " + summaryPath);
        System.out.println("Individual run outputs in: " + sweepDir);
    }
    
    /**
     * Generates injection times for two injections at fixed interval
     */
	private static int[] generateInjectionTimes() {
		int[] times = new int[NUM_INJECTIONS];
		times[0] = FIRST_INJECTION_DAY;
		return times;
	}
    
    /**
     * Generates injection doses - splits total dose equally
     * 
     * @param totalDose Total dose across both injections in nmol
     * @return Array of injection doses in mol (SI units)
     */
	private static double[] generateInjectionDoses(double totalDose) {
		double totalDoseMol = totalDose * 1e-9; // Convert nmol to mol
		double dosePerInjection = totalDoseMol / NUM_INJECTIONS;
		
		double[] doses = new double[NUM_INJECTIONS];
		for (int i = 0; i < NUM_INJECTIONS; i++) {
			doses[i] = dosePerInjection;
		}
		return doses;
	}
            
    /**
     * Runs a single simulation with specified receptor density
     * 
     * @param injectionTimes Array of injection days
     * @param injectionDoses Array of injection amounts in mol
     * @param receptorDensity Receptors per cell in mol
     * @param outputDir Directory for this run's output files
     * @param simulationDays Number of days to simulate
     * @return SimResult with final tumor count
     */
    private static SimResult runSingleSimulation(Rand rng, int[] injectionTimes, 
    									  double[] injectionDoses,
                                          double receptorDensity, String outputDir, 
                                          int simulationDays) throws IOException {
        
        // CRITICAL: Override receptor density BEFORE creating Grid
        // This must happen before any PK calculations
        SimParams.RECEPTORS_PER_CELL_MOL = receptorDensity;
        
        // Initialize simulation
		Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);

        DataLogger logger = new DataLogger();
        DaVinci drawer = new DaVinci(model);
        
        // Set up output directory paths (always, even if not exporting)
        SimParams.OUTPUT_DIR_BASE = outputDir;
        SimParams.OUTPUT_DIR_TUMOUR_IMAGES = outputDir + "/tumour_images";
        SimParams.OUTPUT_DIR_OXYGEN_IMAGES = outputDir + "/oxygen_images";
        SimParams.OUTPUT_DIR_SF_IMAGES = outputDir + "/sf_images";
        
        // Only create directories if exporting images
        if (EXPORT_IMAGES) {
            new File(SimParams.OUTPUT_DIR_TUMOUR_IMAGES).mkdirs();
            new File(SimParams.OUTPUT_DIR_OXYGEN_IMAGES).mkdirs();
            new File(SimParams.OUTPUT_DIR_SF_IMAGES).mkdirs();
        }
        
        // Suppress oxygen field image exports for faster sweeps
        if (!EXPORT_IMAGES) {
            SimParams.EXPORT_OX_IMAGES = false;
            SimParams.EXPORT_TUMOUR_OX_IMAGES = false;
            SimParams.EXPORT_SF_IMAGES = false;
        }
        
        int dayCount = -1;
        model.Init(dayCount, drawer, logger);

		// Develop hypoxia without growth
		if (SimParams.HYPOXIA_DEV_DAYS > 0) {
			model.developHypoxiaWithoutGrowth(SimParams.HYPOXIA_DEV_DAYS, false);
		}

        // Print initial state
        int initialCells = model.countTumorCells();
        int numVessels = model.countVessels();
        double receptorMoles = SimParams.computeReceptorMoles(initialCells, numVessels);
        System.out.printf("  Initial tumor cells: %d (R_T = %.3e nmol, %.1f receptors/cell)%n", 
                         initialCells, receptorMoles * 1e9, receptorDensity * 1e15);

        // Visualization settings
        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};

        // Run simulation
        for (int day = 0; day < simulationDays; day++) {
            dayCount++;
            
            // Check if an injection is scheduled today
            for (int j = 0; j < injectionTimes.length; j++) {
                if (injectionTimes[j] == dayCount) {
                    // Inject into PK model
                    double hotDose = injectionDoses[j] * HOT_FRACTION;
                    double coldDose = injectionDoses[j] * (1.0 - HOT_FRACTION);
                    
                    // Get current PK state and inject
                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    
                    // Inject into central compartment
                    currentPK[0] += hotDose;   // N_cen_hot (mol)
                    currentPK[1] += coldDose;  // N_cen_cold (mol)
                    
                    model.PKStateVariables.set(model.PKStateVariables.size() - 1, currentPK);
                    
                    System.out.printf("  Day %d: Injected %.1f nmol (%.1f hot, %.1f cold)%n",
                                    dayCount, injectionDoses[j] * 1e9, hotDose * 1e9, coldDose * 1e9);
                }
            }
            
            // Step simulation
            model.Step(dayCount);

            // EARLY STOP: tumor eliminated
            double liveNormoxic = model.CurrentCellsPops[SimParams.NORMAL];
            double liveHypoxic  = model.CurrentCellsPops[SimParams.HYPOXIC];
            
            int viableTumor = (int)(liveNormoxic + liveHypoxic);
            
            if (viableTumor == 0) {
                System.out.printf("  Early stop at day %d - tumor eliminated%n", dayCount);
                break;
            }

            // Draw and save visualizations (if enabled)
            if (EXPORT_IMAGES && SimParams.EXPORT_TUMOUR_OX_IMAGES) {
                if (day % 10 == 0) {
                    drawer.gridDraw(visualizationMaskList);
                    
                    double[] valueList = MyUtils.lastElementOfDoubleArrayList(model.DoseRateList);
                    drawer.plot(dayCount, valueList[0]);
                    
                    String imageFile = String.format("%s/day_%03d.png", 
                                                    SimParams.OUTPUT_DIR_TUMOUR_IMAGES, dayCount);
                    logger.saveFigureTotal(imageFile, drawer, dayCount, false, false);
                }
            }
        }
        
        // Save final data files
        String popsFile = outputDir + "/populations.csv";
        String doseFile = outputDir + "/dose.csv";
        String pkFile = outputDir + "/pkStateVariables.csv";
        
        logger.log(model.PopsOverTime, popsFile);
        logger.log(model.DoseRateList, doseFile);
        logger.log(model.PKStateVariables, pkFile);
        
        // Save parameter info for this run
        SimParams.generateParameterReport(outputDir + "/parameters.md");
        SimParams.exportParametersToCSV(outputDir + "/parameters.csv", 
                                injectionTimes, 
                                injectionDoses, 
                                HOT_FRACTION,
                                receptorDensity);
                                
        // Report final viable tumor count
        int finalCount = (int)(model.CurrentCellsPops[SimParams.NORMAL] +
                               model.CurrentCellsPops[SimParams.HYPOXIC]);
        
        return new SimResult(finalCount, NUM_INJECTIONS);
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
