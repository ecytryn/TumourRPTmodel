package TumorRPT;

import HAL.Rand;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.stream.IntStream;

/**
 * Parameter sweep over TWO-INJECTION radiopharmaceutical therapy schedules
 * 
 * Explores 2D parameter space:
 * - Axis 1: Inter-injection interval (days between doses)
 * - Axis 2: Dose distribution skew (front-loading vs back-loading)
 * 
 * Each simulation runs dynamically based on injection schedule and classifies 
 * outcome as CURE or FAILURE based on final tumor cell count.
 * 
 * === HOW TO REFINE THE SWEEP ===
 * 1. Run initial coarse sweep
 * 2. Look at heatmap to identify transition regions (CURE <-> FAILURE boundaries)
 * 3. Edit INTERVALS and SKEWS arrays below to add finer resolution in those regions
 * 4. Change OUTPUT_SUFFIX to avoid overwriting previous results (e.g., "_v2")
 * 5. Re-run sweep
 * 6. Repeat as needed
 */

public class IntervalSkewSweep {
    
    // =======================================================================
    // SWEEP PARAMETERS - Edit these to define the parameter space
    // =======================================================================
    // These parameters are SWEPT (varied systematically) to explore treatment space
    // 
    // INTERVALS: Days between injections (e.g., 30 → injections at 5, 35, 65, 95)
    // SKEWS: Dose distribution offsets in nmol (e.g., 10 → doses of 110, 90 nmol)
    //        Note: Average dose is always 100 nmol regardless of skew
    //        Skew pattern: [100+s, 100-s]
    //
    // For a 2D sweep, each (interval, skew) combination will be run NUM_REPLICATES times
        
    // Current sweep - edit as needed:
    private static final int[] INTERVALS = 
    		IntStream.range(0, 8)  // 9 points
				 .map(i -> (20 + i * 2))  // 20 --> 36 days
				 .toArray();
    private static final double[] SKEWS = 
    		IntStream.range(0, 10)  // 11 points
				 .mapToDouble(i -> (-25e-9 + i * 5e-9))  // -25 --> 25 nmol
				 .toArray();

    // Output suffix - change this when refining to avoid overwriting (e.g., "_v2", "_fine")
    private static final String OUTPUT_SUFFIX = "";
    
    // Number of replicates per parameter combination
    private static final int NUM_REPLICATES = 20;
    
    // =======================================================================
    // FIXED PARAMETERS - Constant across all sweep runs
    // =======================================================================
    // Note: These override SimParams values for sweep runs
    //       Single-run mode (Main.java) uses SimParams values
    
    private static final int NUM_INJECTIONS = 2;
    private static final double BASE_DOSE = 50.0e-9; // mol (nmol converted to mol)
    private static final double HOT_FRACTION = 0.1;
	private static final int FIRST_INJECTION_DAY = 5;
    private static final int MIN_DAYS_AFTER_LAST_INJECTION = 50; // Days to observe after final injection
    
    // Initial tumor size and radiosensitivity (from SimParams defaults)
    // These can be overridden if needed, but default to SimParams values
    
    // Output configuration
    private static final String BASE_OUTPUT_DIR = "results/IntervalSkewSweep/IntervalSkewSweep" + OUTPUT_SUFFIX;
    private static final String SUMMARY_CSV = BASE_OUTPUT_DIR + "/sweep_summary" + OUTPUT_SUFFIX + ".csv";
    
    // Image export control (set to false for faster sweeps)
    private static final boolean EXPORT_IMAGES = false;
    
    public static void main(String[] args) throws IOException {

		SimParams.INITIAL_TUMOR_RADIUS = 100e-6;  // Small vulnerable tumor just below treatment threshold

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

        // CSV header (doses shown in nmol for readability)
        csvWriter.println("interval,skew,dose1_nmol,dose2_nmol,replicate,finalTumorCount,outcome,injectionsUsed");
        
        int totalCombinations = INTERVALS.length * SKEWS.length;
        int totalRuns = totalCombinations * NUM_REPLICATES;
        int currentRun = 0;
        
        System.out.println("=== PARAMETER SWEEP START ===");
        System.out.println("Sweep timestamp: " + timestamp);
        System.out.println("Parameter combinations: " + totalCombinations);
        System.out.println("Replicates per combination: " + NUM_REPLICATES);
        System.out.println("Total runs: " + totalRuns);
        System.out.println("Output directory: " + sweepDir);
        System.out.println();
        
        // Loop over parameter space
        for (int interval : INTERVALS) {
            for (double skew : SKEWS) {
                
                // Generate injection schedule
                int[] injectionTimes = generateInjectionTimes(interval);
                double[] injectionDoses = generateInjectionDoses(skew); // in mol
                
                // Check for positive dose constraint
                boolean allPositive = true;
                for (double dose : injectionDoses) {
                    if (dose < 0) {
                        allPositive = false;
                        break;
                    }
                }
                
                if (!allPositive) {
                    for (int rep = 0; rep < NUM_REPLICATES; rep++) {
                        currentRun++;
                        System.out.printf("[%d/%d] SKIPPED: interval=%d, skew=%.0f nmol, rep=%d (negative dose)%n", 
                                         currentRun, totalRuns, interval, skew * 1e9, rep + 1);
                        csvWriter.printf("%d,%.1f,%.1f,%.1f,%d,%d,%s,%d%n",
                                       interval, skew * 1e9,  // Convert to nmol for display
                                       injectionDoses[0] * 1e9, injectionDoses[1] * 1e9, 
                                       rep + 1, -1, "INVALID", 0);
                        csvWriter.flush();
                    }
                    continue;
                }
                
                // Calculate simulation days dynamically based on injection schedule
                int lastInjectionDay = injectionTimes[NUM_INJECTIONS - 1];
                int simulationDays = lastInjectionDay + MIN_DAYS_AFTER_LAST_INJECTION;
                
                // Run replicates for this parameter combination
                for (int rep = 0; rep < NUM_REPLICATES; rep++) {
                    currentRun++;
                    
                    // Create output directory for this run (includes replicate number)
                    String runDir = String.format("%s/interval_%d_skew_%.0f_rep_%d", 
                                                 sweepDir, interval, skew * 1e9, rep + 1);
                    new File(runDir).mkdirs();
					
					saveGitInfo(runDir);                      
                    
                    System.out.printf("[%d/%d] Running: interval=%d days, skew=%.0f nmol, replicate=%d/%d%n", 
                                     currentRun, totalRuns, interval, skew * 1e9, rep + 1, NUM_REPLICATES);
                    System.out.printf("  Injection times: %d, %d days%n",
                                     injectionTimes[0], injectionTimes[1]);
                    System.out.printf("  Doses: %.1f, %.1f nmol%n",
                                     injectionDoses[0] * 1e9, injectionDoses[1] * 1e9);
                    System.out.printf("  Simulation length: %d days%n", simulationDays);

					// Use replicate-based seed for reproducibility across parameter combinations
					// This means replicate 1 across all parameters uses seed 43, rep 2 uses 44, etc.
					// This creates a "controlled experiment" where only parameters change, not stochasticity
					int seed = 42 + rep;  // rep is 1-indexed, so seeds are 43, 44, 45, 46, 47
					Rand rng = new Rand(seed);
                    
                    // Run simulation
                    SimResult result = runSingleSimulation(rng, injectionTimes, injectionDoses, runDir, simulationDays);

                    // Classify outcome
                    String outcome = (result.finalTumorCount == 0) ? "CURE" : "FAILURE";                

                    System.out.printf("  Result: %d tumor cells remaining -> %s%n%n", 
                                     result.finalTumorCount, outcome);

                    // Write to CSV (doses in nmol for readability)
                    csvWriter.printf("%d,%.1f,%.1f,%.1f,%d,%d,%s,%d%n",
                           interval, skew * 1e9,
                           injectionDoses[0] * 1e9, injectionDoses[1] * 1e9,
                           rep + 1, result.finalTumorCount, outcome,
                           result.injectionsUsed);				
                    csvWriter.flush(); // Flush after each run to preserve data if interrupted
                }
            }
        }
        
        csvWriter.close();
        
        System.out.println("=== PARAMETER SWEEP COMPLETE ===");
        System.out.println("Results saved to: " + summaryPath);
        System.out.println("Individual run outputs in: " + sweepDir);
    }
    
    /**
     * Generates injection times for two injections
     * @param interval Days between injections
     * @return Array of injection times (days)
     */
    private static int[] generateInjectionTimes(int interval) {
        int[] times = new int[NUM_INJECTIONS];
        times[0] = FIRST_INJECTION_DAY;
        times[1] = FIRST_INJECTION_DAY + interval;
        return times;
    }
    
    /**
     * Generates injection doses with linear skew pattern
     * Pattern: [100+s, 100-s] where s is the skew parameter
     * 
     * @param skew Dose skew parameter in mol (input as nmol in main, converted internally)
     * @return Array of injection doses in mol
     */
    private static double[] generateInjectionDoses(double skew) {
        double[] doses = new double[NUM_INJECTIONS];
        doses[0] = BASE_DOSE + skew;
        doses[1] = BASE_DOSE - skew;
        return doses;
    }
    
    /**
     * Runs a single simulation with specified injection schedule
     * 
     * @param injectionTimes Array of injection days
     * @param injectionDoses Array of injection amounts in mol
     * @param outputDir Directory for this run's output files
     * @param simulationDays Number of days to simulate
     * @return SimResult with final tumor count and injections used
     */
    private static SimResult runSingleSimulation(Rand rng, int[] injectionTimes, double[] injectionDoses, 
                                          String outputDir, int simulationDays) throws IOException {
        
        // Initialize simulation (similar to Main.java structure)

        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);
        
        DataLogger logger = new DataLogger();
        DaVinci drawer = new DaVinci(model);
        
        // Set up output subdirectories if exporting images
        
		SimParams.OUTPUT_DIR_BASE = outputDir;
		SimParams.OUTPUT_DIR_TUMOUR_IMAGES = outputDir + "/tumour_images";
		SimParams.OUTPUT_DIR_OXYGEN_IMAGES = outputDir + "/oxygen_images";
		SimParams.OUTPUT_DIR_SF_IMAGES = outputDir + "/sf_images";

        if (EXPORT_IMAGES) {            
            new File(SimParams.OUTPUT_DIR_TUMOUR_IMAGES).mkdirs();
            new File(SimParams.OUTPUT_DIR_OXYGEN_IMAGES).mkdirs();
            new File(SimParams.OUTPUT_DIR_SF_IMAGES).mkdirs();
        }

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

        // Print initial receptor count for verification
        int initialCells = model.countTumorCells();
        int numVessels = model.countVessels();
        double receptorMoles = SimParams.computeReceptorMoles(initialCells, numVessels);
        System.out.printf("  Initial tumor cells: %d (R_T = %.3e nmol)%n", 
                         initialCells, receptorMoles * 1e9);

        // Visualization settings
        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};

        int injectionsUsed = 0;        

        // Run simulation
        for (int day = 0; day < simulationDays; day++) {
            dayCount += 1;
            
            // Check if an injection is scheduled today
            for (int j = 0; j < injectionTimes.length; j++) {
                if (injectionTimes[j] == dayCount) {
                    injectionsUsed++;
                    // Inject into PK model
                    double hotDose = injectionDoses[j] * HOT_FRACTION;
                    double coldDose = injectionDoses[j] * (1.0 - HOT_FRACTION);
                    
                    // Get current PK state and inject
                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    
                    // Inject into central compartment
                    // Index 0 = N_cen_hot, Index 1 = N_cen_cold
                    currentPK[0] += hotDose;   // mol
                    currentPK[1] += coldDose;
                    
                    // Update QSS compartments (indices 2-9)
                    // These will be recalculated in next PK step
                    
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
                System.out.printf("  Early stop at day %d - tumor eliminated after %d injections%n",
                                dayCount, injectionsUsed);
                break;
            }

            // Draw and save visualizations (if enabled)
            if (EXPORT_IMAGES && SimParams.EXPORT_TUMOUR_OX_IMAGES) {
                if (day % 10 == 0) { // Every 10 days for sweep (less frequent than Main)
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
                                SimParams.RECEPTORS_PER_CELL_MOL);
                                        
        // Report final viable tumor count (exclude vessels, necrotic, and apoptotic)
        int finalCount =
            (int)(model.CurrentCellsPops[SimParams.NORMAL] +
                  model.CurrentCellsPops[SimParams.HYPOXIC]);
        
        return new SimResult(finalCount, injectionsUsed);
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

/**
 * Helper class to store simulation results
 */
class SimResult {
    int finalTumorCount;
    int injectionsUsed;
    
    SimResult(int finalTumorCount, int injectionsUsed) {
        this.finalTumorCount = finalTumorCount;
        this.injectionsUsed = injectionsUsed;
    }
}
