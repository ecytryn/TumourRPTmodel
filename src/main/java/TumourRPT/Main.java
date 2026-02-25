package TumorRPT;

import HAL.Rand;

// File operations
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

// Timestamps
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

/**
 * Minimal test harness for frozen tumor PK validation
 * 
 * Tests:
 * - Vessel loading
 * - Tumor initialization
 * - PK model with frozen tumor (no cell dynamics)
 * - Injection mechanism
 * - Dose rate calculations
 * 
 * Run with: ./gradlew run
 */
public class Main {
    
    public static void main(String[] args) throws IOException {
        
		// Parse command line arguments
		String experimentName = "default";
		if (args.length > 0) {
			experimentName = args[0];
		}
		
		// Set parameters based on experiment
		setupExperiment(experimentName);

		// Initialize grid
        Rand rng = new Rand(42);  // Fixed seed for reproducibility
        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);

        DataLogger logger = new DataLogger();

        // Visualization settings
		DaVinci drawer = new DaVinci(model);
	    
        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};

		// Define output dir
		String timestamp = LocalDateTime.now().format(
		   DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss"));
		String outputDir = String.format("results/single_runs/%s_%s",
									   SimParams.EXPERIMENT_NAME, timestamp);
		new File(outputDir).mkdirs();

		// Store in SimParams for global access
		SimParams.OUTPUT_DIR_BASE = outputDir;
		SimParams.OUTPUT_DIR_TUMOUR_IMAGES = outputDir + "/tumour_images";
		SimParams.OUTPUT_DIR_OXYGEN_IMAGES = outputDir + "/oxygen_images";
		SimParams.OUTPUT_DIR_SF_IMAGES = outputDir + "/sf_images";		

		new File(SimParams.OUTPUT_DIR_TUMOUR_IMAGES).mkdirs();
		new File(SimParams.OUTPUT_DIR_OXYGEN_IMAGES).mkdirs();
		new File(SimParams.OUTPUT_DIR_SF_IMAGES).mkdirs();

		// Generate parameter report
		SimParams.generateParameterReport(outputDir + "/parameters.md");
		SimParams.exportParametersToCSV(outputDir + "/parameters.csv",
                                SimParams.INJECTION_SCHEDULE,
                                // Need to create dose array from DOSE_PER_INJECTION:
                                createDoseArray(SimParams.INJECTION_SCHEDULE.length),
                                SimParams.HOT_FRACTION,
                                SimParams.RECEPTORS_PER_CELL_MOL);
		saveGitInfo(outputDir);
		
        // Initialize simulation (creates vessels, seeds tumor, sets up PK)
        int dayCount = -1;
        model.Init(dayCount, drawer, logger);
//        model.Init(dayCount, null, null);

		// Develop hypoxia without growth
		if (SimParams.HYPOXIA_DEV_DAYS > 0) {
			model.developHypoxiaWithoutGrowth(SimParams.HYPOXIA_DEV_DAYS, false);
		}

        // Report initial state
        int initialCells = model.countTumorCells();
        int numVessels = model.countVessels();
        int vesselsNearTumor = model.countVesselsNearTumor();

		System.out.println("\n=== INITIAL STATE ===");
		System.out.println("Tumor cells (2D): " + initialCells);
//		System.out.println("Vessels (near tumor): " + vesselsNearTumor);
//		System.out.println("Vessels (total): " + numVessels);

        double receptorMoles = SimParams.computeReceptorMoles(initialCells, numVessels);
//        System.out.printf("Total receptors: %.3e mol (%.1f nmol)%n", receptorMoles, receptorMoles * 1e9);

//		System.out.println("\n=== GEOMETRY OUTPUT ===");		
//		System.out.printf("V_ec = %.6e m³\n", model.PBPK.getV_ec());
//		System.out.printf("V_v = %.6e m³\n", model.PBPK.getV_v());
//		System.out.printf("R_total = %.6e mol\n", model.PBPK.getR_total());
		System.out.printf("Tumor volume = %.6e m³ (cylindrical extension)\n", model.PBPK.getTumorVolume());
		
		// Calculate derived quantities for comparison
		double R_T_tilde = model.PBPK.getR_total() / model.PBPK.getV_ec();
		double beta = (SimParams.K_OFF + SimParams.K_INT) / SimParams.K_ON;
//		System.out.printf("R_T_tilde = %.6e mol/m³\n", R_T_tilde);
//		System.out.printf("beta = %.6e mol/m³\n", beta);
//		System.out.println("==============================================\n");
        
        // Injection protocol - ALL parameters defined in SimParams.java
        // (No local variables - use SimParams directly for clarity)
        
//        System.out.println("\n--- Injection Protocol ---");
//        System.out.printf("Dose per injection: %.1f nmol%n", SimParams.DOSE_PER_INJECTION * 1e9);
//        System.out.printf("Hot fraction: %.1f%%%n", SimParams.HOT_FRACTION * 100);
        System.out.printf("Injection days: ");
        for (int day : SimParams.INJECTION_SCHEDULE) {
        	System.out.printf("%d  ",day);
		}
        System.out.printf("%n");

        
        // Determine simulation length (from last injection + follow-up period)
        int lastInjectionDay = SimParams.INJECTION_SCHEDULE[SimParams.INJECTION_SCHEDULE.length - 1];
        int totalDays = lastInjectionDay + SimParams.DAYS_AFTER_LAST_INJECTION;
        
        System.out.printf("Simulation length: %d days%n", totalDays);
        
        // Run simulation
        System.out.println("\n--- Running Simulation ---");
        for (int day = 0; day < totalDays; day++) {
            dayCount++;
            
            // Check for injection
            for (int injDay : SimParams.INJECTION_SCHEDULE) {
                if (injDay == dayCount) {
                    // Inject into PK model
                    double hotDose = SimParams.DOSE_PER_INJECTION * SimParams.HOT_FRACTION;
                    double coldDose = SimParams.DOSE_PER_INJECTION * (1.0 - SimParams.HOT_FRACTION);
                    
                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    currentPK[0] += hotDose;   // N_cen_hot
                    currentPK[1] += coldDose;  // N_cen_cold
                    model.PKStateVariables.set(model.PKStateVariables.size() - 1, currentPK);
                    
//                    System.out.printf("Day %d: Injected %.1f nmol (%.1f hot + %.1f cold)%n",
//                                    dayCount, SimParams.DOSE_PER_INJECTION * 1e9, hotDose * 1e9, coldDose * 1e9);
                }
            }

			if (SimParams.EXPORT_TUMOUR_OX_IMAGES) {
            // Draw and save visualizations every 5 days
				if (day % 5 == 0  || (day>=4 && day<=10) || (day>=44 && day<=50)) { // || day<=10
					drawer.gridDraw(visualizationMaskList);
					
					double[] valueList = MyUtils.lastElementOfDoubleArrayList(model.DoseRateList);
					drawer.plot(dayCount, valueList[0]);
					
					// Save tumour image with day number in filename
					String imageFile = String.format("%s/day_%03d.png", SimParams.OUTPUT_DIR_TUMOUR_IMAGES, dayCount);
					if (day == 0 && experimentName.equals("WatchGrow")) {
						logger.saveFigureTotal(imageFile, drawer, dayCount, false, true);
					} else {
						logger.saveFigureTotal(imageFile, drawer, dayCount, false, false);
					}

					// Save SF visualization
					if (SimParams.EXPORT_SF_IMAGES) {
						String sfImageFile = String.format("%s/sf_day_%03d.png", SimParams.OUTPUT_DIR_SF_IMAGES, dayCount);
						logger.saveSFVisualization(sfImageFile, model, dayCount);
					}					
				}
			}


            // Step simulation (tumor frozen, only PK evolves)
            model.Step(dayCount);
            
            // Report every 20 days AND at day 0 for early diagnostics
            if (dayCount % 20 == 0 || dayCount == 0 || dayCount == totalDays - 1) {
//                model.printDiagnostics(dayCount);

                System.out.printf("Day %d%n", dayCount); 

                // Extra detail for first few days
                if (dayCount <= 1) {
//                    System.out.printf("Detailed cell type breakdown:\n");
                    int[] typeCounts = new int[SimParams.NUM_CELL_TYPES];
                    for (Cell cell : model) {
                        if (cell != null && cell.type != SimParams.VESSEL) {
                            typeCounts[cell.type]++;
                        }
                    }
//                    System.out.printf("  Normal=%d, Hypoxic=%d, Necrotic=%d, Apoptotic=%d\n",
//                                     typeCounts[SimParams.NORMAL], typeCounts[SimParams.HYPOXIC],
//                                     typeCounts[SimParams.NECROTIC], typeCounts[SimParams.APOPTOTIC]);
//                    System.out.println();
                }
            }

        }
        
        // Final report
//        System.out.println("\n--- Final State ---");
//        System.out.printf("Simulation completed: %d days%n", totalDays);
//        System.out.printf("PK data points: %d%n", model.PKStateVariables.size());
//        System.out.printf("Dose rate data points: %d%n", model.DoseRateList.size());
        
        double[] finalPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
//        System.out.printf("Final C_cen_hot: %.3e mol/m³ (%.3e nmol/L)%n", 
//                         finalPK[0], finalPK[0] * 1e9);
//        System.out.printf("Final C_ic_hot: %.3e mol/m³ (%.3e nmol/L)%n", 
//                         finalPK[8], finalPK[8] * 1e9);
        
/**
        System.out.println("\n=== Test Complete ===");
        System.out.println("Check output for:");
        System.out.println("1. C_cen decays exponentially after injection");
        System.out.println("2. C_b saturates then decays (receptor binding)");
        System.out.println("3. C_ic accumulates slowly (internalization)");
        System.out.println("4. Dose rate follows activity curve");
        System.out.println("\nIf all looks reasonable, PK core is working!");
*/        
		// Final output of data:
//        DataLogger logger = new DataLogger();
		String popcsvFileName = String.format("%s/populations.csv", outputDir);
		String dosecsvFileName = String.format("%s/dose.csv", outputDir);
		String PKvarscsvFileName = String.format("%s/pkStateVariables.csv", outputDir);

        logger.log(model.PopsOverTime, popcsvFileName); // population over time
        logger.log(model.DoseRateList, dosecsvFileName); // dose rate over time
        logger.log(model.PKStateVariables, PKvarscsvFileName); // PK variables (radioligand in each compartment)

        System.out.println("All data logged!-Close the simulation");
        
    }

	private static double[] createDoseArray(int numInjections) {
		double[] doses = new double[numInjections];
		for (int i = 0; i < numInjections; i++) {
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
			// Git not available or not a git repo
		}
	}

	private static void setupExperiment(String name) {
		switch(name) {
			case "WatchGrow":
				SimParams.setExperiment("WatchATumourGrow", 
					"Watch a tumour grow to see how growth, vessel occlusion, cell type transitions occur.",
					new int[]{195}, 100e-6, 0.06, 0.02, 0);
				break;
			case "WatchGrowLarge":
				SimParams.setExperiment("WatchATumourGrowLarge", 
					"Watch a tumour grow to see how growth, vessel occlusion, cell type transitions occur.",
					new int[]{40}, 950e-6, 0.06, 0.02, 0);
				break;
			case "NormoxicSmall":
				SimParams.setExperiment("NormoxicSmallTumour",
					"Treatment without pre-simulation hypoxia development to demonstrate what happens with treatment of a normoxic tumour. Small initial tumour.",
					new int[]{5}, 100e-6, 0.05, 0.02, 0);
				break;
			case "NormoxicMedium":
				SimParams.setExperiment("NormoxicMediumTumour",
					"Treatment without pre-simulation hypoxia development to demonstrate what happens with treatment of a normoxic tumour. Medium initial tumour.",
					new int[]{5}, 333e-6, 0.05, 0.02, 0);
				break;
			case "NormoxicLarge":
				SimParams.setExperiment("NormoxicLargeTumour",
					"Treatment without pre-simulation hypoxia development to demonstrate what happens with treatment of a normoxic tumour. Large initial tumour.",
					new int[]{5}, 950e-6, 0.05, 0.02, 0);
				break;
			case "HypoxicSmall":
				SimParams.setExperiment("HypoxicSmallTumour",
					"Treatment with pre-simulation hypoxia development to demonstrate what happens with treatment of a hypoxic tumour. Small initial tumour.",
					new int[]{5}, 100e-6, 0.05, 0.02, 40);
				break;
			case "HypoxicMedium":
				SimParams.setExperiment("HypoxicMediumTumour",
					"Treatment with pre-simulation hypoxia development to demonstrate what happens with treatment of a hypoxic tumour. Medium initial tumour.",
					new int[]{5}, 333e-6, 0.05, 0.02, 40);
				break;
			case "HypoxicLarge":
				SimParams.setExperiment("HypoxicLargeTumour",
					"Treatment with pre-simulation hypoxia development to demonstrate what happens with treatment of a hypoxic tumour. Large initial tumour",
					new int[]{5}, 950e-6, 0.05, 0.02, 40);
				break;
			case "Reoxygenation":
				SimParams.setExperiment("ReoxygenationExperiment",
					"Set alpha_hypoxic = beta_hypoxic =0. This allows us to see how, once reoxygenated, hypoxic cells convert back to normoxic and die.",
					new int[]{5}, 333e-6, 0.0, 0.0, 40);  // Zero hypoxic sensitivity!
				break;
			case "CustomRun":
				SimParams.setExperiment("CustomRunToMakeAFig",
					"This case is for one-off runs to create specific plots for a figure.",
					new int[]{5}, 100e-6, 0.05, 0.02, 40);
				break;
			default:
				// Use defaults from SimParams
				break;
		}
	}
}
