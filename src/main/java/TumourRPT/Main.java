package TumorRPT;

import HAL.Rand;
import java.io.IOException;

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
        
        // Initialize grid
        Rand rng = new Rand(42);  // Fixed seed for reproducibility
        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);

        DataLogger logger = new DataLogger();

        // Visualization settings
        DaVinci drawer = new DaVinci(model);
        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};
		// Define output dir
		String outputDir = String.format("results/single_runs");
                
        // Initialize simulation (creates vessels, seeds tumor, sets up PK)
        int dayCount = -1;
        model.Init(dayCount, drawer, logger);
//        model.Init(dayCount, null, null);

        // Report initial state
        int initialCells = model.countTumorCells();
        int numVessels = model.countVessels();
        int vesselsNearTumor = model.countVesselsNearTumor();

		System.out.println("\n=== INITIAL STATE ===");
		System.out.println("Tumor cells (2D): " + initialCells);
		System.out.println("Vessels (near tumor): " + vesselsNearTumor);
		System.out.println("Vessels (total): " + numVessels);

        double receptorMoles = SimParams.computeReceptorMoles(initialCells, numVessels);
        System.out.printf("Total receptors: %.3e mol (%.1f nmol)%n", receptorMoles, receptorMoles * 1e9);

		System.out.println("\n=== GEOMETRY OUTPUT ===");		
		System.out.printf("V_ec = %.6e m³\n", model.PBPK.getV_ec());
		System.out.printf("V_v = %.6e m³\n", model.PBPK.getV_v());
		System.out.printf("R_total = %.6e mol\n", model.PBPK.getR_total());
		System.out.printf("Tumor volume = %.6e m³\n", model.PBPK.getTumorVolume());
		
		// Calculate derived quantities for comparison
		double R_T_tilde = model.PBPK.getR_total() / model.PBPK.getV_ec();
		double beta = (SimParams.K_OFF + SimParams.K_INT) / SimParams.K_ON;
		System.out.printf("R_T_tilde = %.6e mol/m³\n", R_T_tilde);
		System.out.printf("beta = %.6e mol/m³\n", beta);
		System.out.println("==============================================\n");
        
        // Injection protocol
        int[] injectionDays = {5, 35, 65, 95};
        double injectionDose = 100e-9;  // mol (100 nmol)
        double hotFraction = 0.1;
        
        System.out.println("\n--- Injection Protocol ---");
        System.out.printf("Dose per injection: %.1f nmol%n", injectionDose * 1e9);
        System.out.printf("Hot fraction: %.1f%%%n", hotFraction * 100);
        System.out.printf("Injection days: ");
        for (int day : injectionDays) {
        	System.out.printf("%d  ",day);
		}
        System.out.printf("%n");

        
        // Determine simulation length
        int lastInjectionDay = injectionDays[injectionDays.length - 1];
        int daysAfterLastInjection = 30;
        int totalDays = lastInjectionDay + daysAfterLastInjection;
        
        System.out.printf("Simulation length: %d days%n", totalDays);
        
        // Run simulation
        System.out.println("\n--- Running Simulation ---");
        for (int day = 0; day < totalDays; day++) {
            dayCount++;
            
            // Check for injection
            for (int injDay : injectionDays) {
                if (injDay == dayCount) {
                    // Inject into PK model
                    double hotDose = injectionDose * hotFraction;
                    double coldDose = injectionDose * (1.0 - hotFraction);
                    
                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    currentPK[0] += hotDose / SimParams.V_CENTRAL;   // C_cen_hot
                    currentPK[1] += coldDose / SimParams.V_CENTRAL;  // C_cen_cold
                    model.PKStateVariables.set(model.PKStateVariables.size() - 1, currentPK);
                    
                    System.out.printf("Day %d: Injected %.1f nmol (%.1f hot + %.1f cold)%n",
                                    dayCount, injectionDose * 1e9, hotDose * 1e9, coldDose * 1e9);
                }
            }

            // Step simulation (tumor frozen, only PK evolves)
            model.Step(dayCount);
            
            // Report every 20 days AND at day 0 for early diagnostics
            if (dayCount % 20 == 0 || dayCount == 0 || dayCount == totalDays - 1) {
//                model.printDiagnostics(dayCount);
                
                // Extra detail for first few days
                if (dayCount <= 1) {
                    System.out.printf("Detailed cell type breakdown:\n");
                    int[] typeCounts = new int[SimParams.NUM_CELL_TYPES];
                    for (Cell cell : model) {
                        if (cell != null && cell.type != SimParams.VESSEL) {
                            typeCounts[cell.type]++;
                        }
                    }
                    System.out.printf("  Normal=%d, Hypoxic=%d, Necrotic=%d, Apoptotic=%d\n",
                                     typeCounts[SimParams.NORMAL], typeCounts[SimParams.HYPOXIC],
                                     typeCounts[SimParams.NECROTIC], typeCounts[SimParams.APOPTOTIC]);
                    System.out.println();
                }
            }

            // Draw and save visualizations every 5 days
            if (day % 5 == 0 || day<=10) {
                drawer.gridDraw(visualizationMaskList);
//                drawer.gridDrawAge(visualizationMaskList);
                
                double[] valueList = MyUtils.lastElementOfDoubleArrayList(model.DoseRateList);
                drawer.plot(dayCount, valueList[0]);
                
                // Save with day number in filename
                String imageFile = String.format("%s/day_%03d.png", outputDir, dayCount);
                logger.saveFigureTotal(imageFile, drawer, dayCount, false);
            }

        }
        
        // Final report
        System.out.println("\n--- Final State ---");
        System.out.printf("Simulation completed: %d days%n", totalDays);
//        System.out.printf("PK data points: %d%n", model.PKStateVariables.size());
//        System.out.printf("Dose rate data points: %d%n", model.DoseRateList.size());
        
        // Print radiobiology validation report
        model.radioBio.printValidationReport();
        
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
        logger.log(model.PopsOverTime, "results/single_runs/populations.csv"); // population over time
        logger.log(model.DoseRateList, "results/single_runs/dose.csv"); // dose rate over time
        logger.log(model.PKStateVariables, "results/single_runs/pkStateVariables.csv"); // PK variables (radioligand in each compartment)

        System.out.println("All data logged!-Close the simulation");
        
    }
}
