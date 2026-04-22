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

/**
 * Sweep over initial tumour size to characterize the cure/failure threshold.
 *
 * Motivation: single-seed figure runs found threshold between 19 µm (7 cells,
 * failure) and 20 µm (11 cells, success) with seed 42. This sweep runs
 * multiple replicates per size to determine whether that threshold is
 * seed-dependent and to quantify the transition as a cure rate.
 *
 * Initial radii tested: 10, 15, 20, 25 µm → 3, 7, 11, 19 cells (from
 * CircleHood geometry with one vessel collision at centre).
 *
 * All other parameters match the figure runs:
 * - 605 cap/mm², 40 pre-sim days, TJ retention table
 * - 100 nmol single injection, 10% hot fraction, injection day 5
 * - 40-day follow-up after injection
 */
public class TumourSizeThresholdSweep {

    // =======================================================================
    // SWEEP PARAMETERS
    // =======================================================================

    // Initial radii in um to be converted below into metres. 11 and 21 are to avoid floating point issues at edge cases.
    private static final double[] INITIAL_RADII_UM = {11, 15, 21, 25, 30, 40, 50};
//    private static final double[] INITIAL_RADII_UM = {70, 90, 110, 130, 160, 200};

    private static final int NUM_REPLICATES = 150;

//    private static final double[] INITIAL_RADII_UM = {21, 25, 28, 30};
//    private static final int NUM_REPLICATES = 5;

    // =======================================================================
    // FIXED PARAMETERS - match figure runs exactly
    // =======================================================================

    private static final double DOSE_NMOL      = 50.0;       // nmol
//    private static final double DOSE_NMOL      = 100.0;       // nmol
    private static final double HOT_FRACTION   = 0.1;
    private static final int    INJECTION_DAY  = 5;
    private static final int    FOLLOW_UP_DAYS = 40;
    private static final int    SIMULATION_DAYS = INJECTION_DAY + FOLLOW_UP_DAYS;

    private static final String BASE_OUTPUT_DIR =
        "results/TumourSizeThresholdSweep/TumourSizeThresholdSweep";

    private static final boolean EXPORT_IMAGES = false;

    // =======================================================================

    public static void main(String[] args) throws IOException {

        // ===================================================================
        // RUN CONFIGURATION
        // ===================================================================
        SimParams.HYPOXIA_DEV_DAYS     = 40;
//        SimParams.VESSEL_DENSITY_CONFIG = "605";
        SimParams.VESSEL_DENSITY_CONFIG = "420";
        // ===================================================================

        String timestamp = LocalDateTime.now().format(
            DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss"));

        String sweepDir = BASE_OUTPUT_DIR + "_" + timestamp;
        new File(sweepDir).mkdirs();

        String summaryPath = sweepDir + "/sweep_summary.csv";
        PrintWriter csvWriter;
        try {
            csvWriter = new PrintWriter(new FileWriter(summaryPath));
        } catch (IOException e) {
            throw new RuntimeException("Could not open summary CSV: " + summaryPath, e);
        }

        csvWriter.println("radius_um,cell_count,replicate,finalTumorCount,outcome");

        int totalRuns = INITIAL_RADII_UM.length * NUM_REPLICATES;
        int currentRun = 0;

        System.out.println("=== TUMOUR SIZE THRESHOLD SWEEP ===");
        System.out.println("Timestamp: " + timestamp);
        System.out.println("Radii (um): 10, 15, 20, 25");
        System.out.println("Replicates per radius: " + NUM_REPLICATES);
        System.out.println("Total runs: " + totalRuns);
        System.out.println("Output: " + sweepDir);
        System.out.println();

        double[] injectionDoses = new double[]{DOSE_NMOL * 1e-9};  // mol
        int[]    injectionTimes = new int[]{INJECTION_DAY};

        for (double radius_um : INITIAL_RADII_UM) {

            SimParams.INITIAL_TUMOR_RADIUS = radius_um * 1e-6;  // m

            System.out.printf("--- Radius = %.0f µm ---\n", radius_um);

            int cellCountObserved = -1;  // will be set from first replicate

            for (int rep = 0; rep < NUM_REPLICATES; rep++) {
                currentRun++;

				int seed = rep * 999983 + currentRun;
                Rand rng = new Rand(seed);

                String runDir = String.format("%s/radius_%.0fum_rep_%d",
                                             sweepDir, radius_um, rep + 1);
                new File(runDir).mkdirs();
                saveGitInfo(runDir);

                System.out.printf("[%d/%d] radius=%.0f µm, rep=%d/%d\n",
                                  currentRun, totalRuns,
                                  radius_um, rep + 1, NUM_REPLICATES);

                SimResult result = runSingleSimulation(rng, injectionTimes,
                                                       injectionDoses, runDir);

                // Record observed cell count from first replicate for CSV
                if (rep == 0) {
                    cellCountObserved = result.initialCellCount;
                    System.out.printf("  Initial cell count: %d\n", cellCountObserved);
                }

                String outcome = (result.finalTumorCount == 0) ? "CURE" : "FAILURE";

                System.out.printf("  Result: %d cells remaining -> %s\n\n",
                                  result.finalTumorCount, outcome);

                csvWriter.printf("%.0f,%d,%d,%d,%s\n",
                                 radius_um, cellCountObserved,
                                 rep + 1, result.finalTumorCount, outcome);
                csvWriter.flush();
            }
        }

        csvWriter.close();
        System.out.println("=== SWEEP COMPLETE ===");
        System.out.println("Results: " + summaryPath);
    }

    private static SimResult runSingleSimulation(Rand rng, int[] injectionTimes,
                                                  double[] injectionDoses,
                                                  String outputDir) throws IOException {

        SimParams.EXPORT_OX_IMAGES      = false;
        SimParams.EXPORT_TUMOUR_OX_IMAGES = false;
        SimParams.EXPORT_SF_IMAGES       = false;

        SimParams.OUTPUT_DIR_BASE          = outputDir;
        SimParams.OUTPUT_DIR_TUMOUR_IMAGES = outputDir + "/tumour_images";
        SimParams.OUTPUT_DIR_OXYGEN_IMAGES = outputDir + "/oxygen_images";
        SimParams.OUTPUT_DIR_SF_IMAGES     = outputDir + "/sf_images";

        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);
        DataLogger logger = new DataLogger();
        DaVinci drawer = new DaVinci(model);

        int dayCount = -1;
        model.Init(dayCount, drawer, logger);
                  
        if (SimParams.HYPOXIA_DEV_DAYS > 0) {
            model.developHypoxiaWithoutGrowth(SimParams.HYPOXIA_DEV_DAYS, false);
        }

        int initialCells = model.countTumorCells();
        System.out.printf("  Seeded cells: %d\n", initialCells);

        for (int day = 0; day < SIMULATION_DAYS; day++) {
            dayCount++;

            for (int j = 0; j < injectionTimes.length; j++) {
                if (injectionTimes[j] == dayCount) {
                    double hotDose  = injectionDoses[j] * HOT_FRACTION;
                    double coldDose = injectionDoses[j] * (1.0 - HOT_FRACTION);
                    double[] currentPK = model.PKStateVariables.get(
                                             model.PKStateVariables.size() - 1);
                    currentPK[0] += hotDose;
                    currentPK[1] += coldDose;
                    model.PKStateVariables.set(
                        model.PKStateVariables.size() - 1, currentPK);
                    System.out.printf("  Day %d: injected %.1f nmol\n",
                                      dayCount, injectionDoses[j] * 1e9);
                }
            }

            model.Step(dayCount);

            // Early stop if tumour eliminated
            int viable = (int)(model.CurrentCellsPops[SimParams.NORMAL] +
                               model.CurrentCellsPops[SimParams.HYPOXIC]);
            if (viable == 0) {
                System.out.printf("  Early stop day %d - tumour eliminated\n", dayCount);
                break;
            }
        }

        // Save data
        logger.log(model.PopsOverTime,      outputDir + "/populations.csv");
        logger.log(model.DoseRateList,      outputDir + "/dose.csv");
        logger.log(model.PKStateVariables,  outputDir + "/pkStateVariables.csv");

        SimParams.generateParameterReport(outputDir + "/parameters.md");
        SimParams.exportParametersToCSV(outputDir + "/parameters.csv",
                                        injectionTimes, injectionDoses,
                                        HOT_FRACTION,
                                        SimParams.RECEPTORS_PER_CELL_MOL);

        int finalCount = (int)(model.CurrentCellsPops[SimParams.NORMAL] +
                               model.CurrentCellsPops[SimParams.HYPOXIC]);

        return new SimResult(finalCount, 1, initialCells);
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
            // Git not available - silently skip
        }
    }
}