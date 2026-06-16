package TumourRPT;

import HAL.Rand;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.stream.IntStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Parameter sweep over injected amount and per-cell receptor content.
 *
 * Usage:
 *   java DoseReceptorSweep
 *       -> New timestamped folder, reps 1..NUM_REPLICATES
 *
 *   java DoseReceptorSweep 2026-01-27_15-30-00
 *       -> Appends NUM_REPLICATES more reps to existing folder,
 *          starting from (max existing rep + 1) per parameter combination.
 *
 * Explores 2D parameter space:
 *   Axis 1: Injected amount (nmol) - single injection
 *   Axis 2: Per-cell receptor content (mol/cell)
 */
public class DoseReceptorSweep {

    // =======================================================================
    // SWEEP PARAMETERS - Edit these to define the parameter space
    // =======================================================================

    private static final double[] DOSES =
        IntStream.range(1, 17)
                 .mapToDouble(i -> (0 + i * 12.5))
                 .toArray();

    // Per-cell receptor content: expressed as multiples of SimParams.RECEPTORS_PER_CELL_MOL
    private static final double[] RECEPTOR_DENSITIES =
        IntStream.range(0, 25)
                 .mapToDouble(i -> (0.68 + i * 0.04) * SimParams.RECEPTORS_PER_CELL_MOL)
                 .toArray();

    // Number of replicates to add per run (new folder: reps 1..N; append: next N reps)
    private static final int NUM_REPLICATES = 1;

    // =======================================================================
    // FIXED PARAMETERS
    // =======================================================================

    private static final int NUM_INJECTIONS = 1;
    private static final double HOT_FRACTION = 0.1;
    private static final int FIRST_INJECTION_DAY = 5;
    private static final int MIN_DAYS_AFTER_LAST_INJECTION = 60;
    private static final String BASE_OUTPUT_DIR = "results/DoseReceptorSweep/DoseReceptorSweep";
    private static final boolean EXPORT_IMAGES = false;

    public static void main(String[] args) throws IOException {

        // ===================================================================
        // RUN CONFIGURATION
        // ===================================================================
        SimParams.INITIAL_TUMOR_RADIUS = 21e-6;
        SimParams.HYPOXIA_DEV_DAYS = 40;
        SimParams.VESSEL_DENSITY_CONFIG = "605";
        // ===================================================================

        boolean appendMode = (args.length >= 1);
        String sweepDir;
        String summaryPath;
        int startRepOffset; // global rep offset: new reps are numbered startRepOffset+1 .. startRepOffset+NUM_REPLICATES

        if (appendMode) {
            // --- APPEND MODE ---
            String timestamp = args[0];
            sweepDir = BASE_OUTPUT_DIR + "_" + timestamp;

            if (!new File(sweepDir).exists()) {
                System.err.println("ERROR: Directory not found: " + sweepDir);
                System.exit(1);
            }

            summaryPath = sweepDir + "/sweep_summary.csv";

            if (!new File(summaryPath).exists()) {
                System.err.println("ERROR: sweep_summary.csv not found in: " + sweepDir);
                System.exit(1);
            }

            // Determine the global max replicate number already in the file
            startRepOffset = findMaxReplicate(summaryPath);
            System.out.println("=== DOSE-RECEPTOR SWEEP APPEND MODE ===");
            System.out.println("Appending to: " + sweepDir);
            System.out.println("Existing replicates found: " + startRepOffset);
            System.out.println("Adding replicates: " + (startRepOffset + 1) + " to " + (startRepOffset + NUM_REPLICATES));

        } else {
            // --- NEW FOLDER MODE ---
            String timestamp = LocalDateTime.now().format(
                DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss"));
            sweepDir = BASE_OUTPUT_DIR + "_" + timestamp;
            new File(sweepDir).mkdirs();
            summaryPath = sweepDir + "/sweep_summary.csv";
            startRepOffset = 0;

            // Write CSV header (only for new files)
            try (PrintWriter hw = new PrintWriter(new FileWriter(summaryPath, false))) {
                hw.println("dose_nmol,receptors_per_cell_mol,dose1_nmol,replicate,finalTumorCount,outcome");
            }

            System.out.println("=== DOSE-RECEPTOR SWEEP START ===");
            System.out.println("Output directory: " + sweepDir);
        }

        int totalCombinations = DOSES.length * RECEPTOR_DENSITIES.length;
        int totalRuns = totalCombinations * NUM_REPLICATES;
        int currentRun = 0;

        System.out.println("Parameter combinations: " + totalCombinations);
        System.out.println("Replicates per combination this run: " + NUM_REPLICATES);
        System.out.println("Total simulations this run: " + totalRuns);
        System.out.println();

        // Open CSV in append mode
        PrintWriter csvWriter = new PrintWriter(new FileWriter(summaryPath, true));

        try {
            for (double dose : DOSES) {
                for (double receptorDensity : RECEPTOR_DENSITIES) {

                    int[] injectionTimes = generateInjectionTimes();
                    double[] injectionDoses = generateInjectionDoses(dose);

                    int lastInjectionDay = injectionTimes[NUM_INJECTIONS - 1];
                    int simulationDays = lastInjectionDay + MIN_DAYS_AFTER_LAST_INJECTION;

                    for (int repIdx = 0; repIdx < NUM_REPLICATES; repIdx++) {
                        currentRun++;
                        int globalRep = startRepOffset + repIdx + 1; // 1-based, globally unique

                        String runDir = String.format("%s/dose_%.0f_recep_%.2e_rep_%d",
                                                     sweepDir, dose, receptorDensity, globalRep);
                        new File(runDir).mkdirs();
                        saveGitInfo(runDir);

                        System.out.printf("[%d/%d] dose=%.0f nmol, receptors=%.2e mol/cell, rep=%d%n",
                                         currentRun, totalRuns, dose, receptorDensity, globalRep);
                        System.out.printf("  Simulation length: %d days%n", simulationDays);

                        // Seed is globally unique per replicate number (same scheme as before)
						int seed = globalRep * 999983 + currentRun;
                        Rand rng = new Rand(seed);

                        SimResult result = runSingleSimulation(rng, injectionTimes, injectionDoses,
                                                              receptorDensity, runDir, simulationDays);

                        String outcome = (result.finalTumorCount == 0) ? "CURE" : "FAILURE";

                        System.out.printf("  Result: %d tumor cells remaining -> %s%n%n",
                                         result.finalTumorCount, outcome);

                        csvWriter.printf("%.1f,%.3e,%.1f,%d,%d,%s%n",
                               dose, receptorDensity,
                               injectionDoses[0] * 1e9,
                               globalRep, result.finalTumorCount, outcome);
                        csvWriter.flush();
                    }
                }
            }
        } finally {
            csvWriter.close();
        }

        System.out.println("=== DOSE-RECEPTOR SWEEP COMPLETE ===");
        System.out.println("Results saved to: " + summaryPath);
    }

    /**
     * Reads an existing sweep_summary.csv and returns the maximum replicate
     * number found across all parameter combinations.
     */
    private static int findMaxReplicate(String csvPath) throws IOException {
        int maxRep = 0;
        try (BufferedReader br = new BufferedReader(new java.io.FileReader(csvPath))) {
            String line = br.readLine(); // skip header
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] parts = line.split(",");
                if (parts.length >= 4) {
                    try {
                        int rep = Integer.parseInt(parts[3].trim());
                        if (rep > maxRep) maxRep = rep;
                    } catch (NumberFormatException e) {
                        // skip malformed rows
                    }
                }
            }
        }
        return maxRep;
    }

    private static int[] generateInjectionTimes() {
        int[] times = new int[NUM_INJECTIONS];
        times[0] = FIRST_INJECTION_DAY;
        return times;
    }

    private static double[] generateInjectionDoses(double totalDose) {
        double totalDoseMol = totalDose * 1e-9;
        double dosePerInjection = totalDoseMol / NUM_INJECTIONS;
        double[] doses = new double[NUM_INJECTIONS];
        for (int i = 0; i < NUM_INJECTIONS; i++) {
            doses[i] = dosePerInjection;
        }
        return doses;
    }

    private static SimResult runSingleSimulation(Rand rng, int[] injectionTimes,
                                          double[] injectionDoses,
                                          double receptorDensity, String outputDir,
                                          int simulationDays) throws IOException {

        SimParams.RECEPTORS_PER_CELL_MOL = receptorDensity;

        Grid model = new Grid(SimParams.GRID_SIZE, SimParams.GRID_SIZE, rng, null);
        DataLogger logger = new DataLogger();
        DaVinci drawer = new DaVinci(model);

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

        if (SimParams.HYPOXIA_DEV_DAYS > 0) {
            model.developHypoxiaWithoutGrowth(SimParams.HYPOXIA_DEV_DAYS, false);
        }

        int initialCells = model.countTumorCells();
        int numVessels = model.countVessels();
        double receptorMoles = SimParams.computeReceptorMoles(initialCells, numVessels);
        System.out.printf("  Initial tumor cells: %d (R_T = %.3e nmol, %.1f receptors/cell)%n",
                         initialCells, receptorMoles * 1e9, receptorDensity * 1e15);

        boolean[] visualizationMaskList = {true, true, true, true, true, true, true};

        for (int day = 0; day < simulationDays; day++) {
            dayCount++;

            for (int j = 0; j < injectionTimes.length; j++) {
                if (injectionTimes[j] == dayCount) {
                    double hotDose = injectionDoses[j] * HOT_FRACTION;
                    double coldDose = injectionDoses[j] * (1.0 - HOT_FRACTION);

                    double[] currentPK = model.PKStateVariables.get(model.PKStateVariables.size() - 1);
                    currentPK[0] += hotDose;
                    currentPK[1] += coldDose;
                    model.PKStateVariables.set(model.PKStateVariables.size() - 1, currentPK);

                    System.out.printf("  Day %d: Injected %.1f nmol (%.1f hot, %.1f cold)%n",
                                    dayCount, injectionDoses[j] * 1e9, hotDose * 1e9, coldDose * 1e9);
                }
            }

            model.Step(dayCount);

            double liveNormoxic = model.CurrentCellsPops[SimParams.NORMAL];
            double liveHypoxic  = model.CurrentCellsPops[SimParams.HYPOXIC];
            int viableTumor = (int)(liveNormoxic + liveHypoxic);

            if (viableTumor == 0) {
                System.out.printf("  Early stop at day %d - tumor eliminated%n", dayCount);
                break;
            }

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

        String popsFile = outputDir + "/populations.csv";
        String doseFile = outputDir + "/dose.csv";
        String pkFile = outputDir + "/pkStateVariables.csv";

        logger.log(model.PopsOverTime, popsFile);
        logger.log(model.DoseRateList, doseFile);
        logger.log(model.PKStateVariables, pkFile);

        SimParams.generateParameterReport(outputDir + "/parameters.md");
        SimParams.exportParametersToCSV(outputDir + "/parameters.csv",
                                injectionTimes,
                                injectionDoses,
                                HOT_FRACTION,
                                receptorDensity);

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
            // Git not available - silently skip
        }
    }
}