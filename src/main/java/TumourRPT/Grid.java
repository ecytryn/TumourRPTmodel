package TumorRPT;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.Grid2Ddouble;
import HAL.Rand;

import java.io.InputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.HashSet;
import java.util.Set;

import static HAL.Util.CircleHood;

/**
 * Main spatial grid for tumor simulation
 * 
 * Manages:
 * - Cell agents (tumor, healthy, vessels)
 * - Oxygen field (Grid2Ddouble)
 * - PK model integration
 * - Radiobiology calculations
 * - Output data collection
 */
public class Grid extends AgentGrid2D<Cell> {
	        
	public Rand rng;
    public Main main;
    public Oxygen oxygen;
    public PK PBPK;
    
    public Grid2Ddouble oxygenGrid;
    public RadioBio radioBio;
        
    // Data collection
    public ArrayList<double[]> DoseRateList = new ArrayList<>();
    public ArrayList<double[]> PKStateVariables = new ArrayList<>();
    public ArrayList<double[]> PopsOverTime = new ArrayList<>();
    public double[] CurrentCellsPops = new double[SimParams.NUM_CELL_TYPES];
    
    // Vessel tracking
    public ArrayList<Integer> vesselsIndex = new ArrayList<>();
    public int currentVesselCounts;
    
    // Tumor geometry (used by PK for receptor calculation)
    public double tempAverageRad;
    
    // Oxygen visualization
//    public OxygenCrossSectionPlotter oxygenPlotter;

    /**
     * Initialize grid with dimensions and random number generator
     */
    public Grid(int xDim, int yDim, Rand rng, Main main) {
        super(xDim, yDim, Cell.class);
        this.rng = rng;
        this.main = main;
        this.oxygen = new Oxygen(this);
        this.PBPK = new PK(this);
        this.radioBio = new RadioBio(this);
        this.oxygenGrid = new Grid2Ddouble(xDim, yDim);
    }

    /**
     * Initialize simulation: create vessels, seed tumor, set up PK
     */
    void Init(int currentDay, DaVinci drawer, DataLogger logger) throws IOException {
        // Create vessels
        this.currentVesselCounts = this.GenVessels();

        // Count blocked vessels
        int blocked = 0;
        for (int idx : vesselsIndex) {
            if (GetAgent(idx).blockedVessel) blocked++;
        }
        System.out.printf("Vessels: %d total, %d blocked (%.1f%%)%n", 
                         vesselsIndex.size(), blocked, 100.0*blocked/vesselsIndex.size());

        // CALIBRATE oxygen boundary condition from healthy tissue
        // Must happen after vessels but before tumor seeding
		if (SimParams.USE_CALIBRATED_BC) {
			this.oxygen.calibrateBoundaryCondition();
		}       
        // Seed tumor
        this.SeedTumor();

        // Initial oxygen solve
        double r_average = updatePKGeometry();
        this.oxygen.UpdateSteadyStateOxygen(currentDay, r_average);

        // Initialize PK state vector
        // [N_cen_hot, N_cen_cold, N_v_hot, N_v_cold, N_ec_hot, N_ec_cold, 
        //  N_b_hot, N_b_cold, N_ic_hot, N_ic_cold, A_blob]
        PKStateVariables.add(new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                         0.0, 0.0, 0.0, 0.0, 1.0}); 
        DoseRateList.add(new double[]{0.0});
        PopsOverTime.add(new double[SimParams.NUM_CELL_TYPES]);
    }

    /**
     * Main simulation step - advances one day

    public void Step(int currentDay) {
        int counter = 0;
        int hourCount = -1;
        
        // Update per hour, not all cells at once
        int perHourCellUpdateBudget = Math.floorDiv((int) this.Pop(), 23);  // 24 hours - 1
        
        ArrayList<double[]> survivalProbLookupTable = new ArrayList<>();
        
        for (Cell cell : this) {
            if (counter % perHourCellUpdateBudget == 0) {
                // Hour boundary - update fields and PK
                hourCount++;
				SimParams.updateGlobalTime(currentDay, hourCount);                
                // Update PK geometry based on current tumor size
                double r_average = updatePKGeometry();
                
                // Solve oxygen steady state
                this.oxygen.UpdateSteadyStateOxygen(currentDay, r_average);
                
                // Step PK forward (time now in seconds!)
                double t_0 = hourCount * SimParams.TIME_STEP;
                double t_f = (hourCount + 1) * SimParams.TIME_STEP;
                this.PBPK.DoseRateCalc(t_0, t_f, currentDay, hourCount);
                
                // Compute radiation survival probabilities
                survivalProbLookupTable = this.radioBio.SurvivalProbLookupTableCalc(currentDay, hourCount);

                // Log population
                PopsOverTime.add(copyArray(CurrentCellsPops));
            }
            
            // Update cell (unless tumor is frozen for testing)
            if (!SimParams.FREEZE_TUMOR) {
                cell.Step(hourCount, currentDay, survivalProbLookupTable);
            }
            counter++;
        }
        
        CleanAgents();
        ShuffleAgents(this.rng);
    }
	*/
    public void Step(int currentDay) {
        long stepStart = System.nanoTime();
        
        int counter = 0;
        int hourCount = -1;
        
        int perHourCellUpdateBudget = Math.floorDiv((int) this.Pop(), 23);
        ArrayList<double[]> survivalProbLookupTable = new ArrayList<>();
        
        long oxygenTime = 0, pkTime = 0, radioBioTime = 0, cellIterTime = 0;
        int oxygenCalls = 0, pkCalls = 0;
        
        for (Cell cell : this) {
            if (counter % perHourCellUpdateBudget == 0) {
                hourCount++;
                SimParams.updateGlobalTime(currentDay, hourCount);
                
                // Time oxygen solver
                long t1 = System.nanoTime();
                double r_average = updatePKGeometry();
                this.oxygen.UpdateSteadyStateOxygen(currentDay, r_average);
                oxygenTime += System.nanoTime() - t1;
                oxygenCalls++;
                
                // Time PK solver
                long t2 = System.nanoTime();
                double t_0 = hourCount * SimParams.TIME_STEP;
                double t_f = (hourCount + 1) * SimParams.TIME_STEP;
                this.PBPK.DoseRateCalc(t_0, t_f, currentDay, hourCount);
                pkTime += System.nanoTime() - t2;
                pkCalls++;
                
                // Time radiobiology
                long t3 = System.nanoTime();
// ----->          Temporary work-around for the growth in compute time - only fixes FREEZE_TUMOR=true case.
//				if (!SimParams.FREEZE_TUMOR) {
					survivalProbLookupTable = this.radioBio.SurvivalProbLookupTableCalc(currentDay, hourCount);
//				} else {
//					survivalProbLookupTable = new ArrayList<>();  // Empty - not used
//				}
                radioBioTime += System.nanoTime() - t3;
                
                PopsOverTime.add(copyArray(CurrentCellsPops));
            }
            
            // Time cell iteration
            long t4 = System.nanoTime();
            if (!SimParams.FREEZE_TUMOR) {
                cell.Step(hourCount, currentDay, survivalProbLookupTable);
            }
            cellIterTime += System.nanoTime() - t4;
            
            counter++;
        }
        
        long cleanStart = System.nanoTime();
        CleanAgents();
        ShuffleAgents(this.rng);
        long cleanTime = System.nanoTime() - cleanStart;
        
        long totalStepTime = System.nanoTime() - stepStart;
        
        // Print timing breakdown every 5 days
        if (currentDay % 5 == 6) {
            System.out.printf("\nDay %d timing breakdown:\n", currentDay);
            System.out.printf("  Oxygen:    %.3f s (%d calls, avg %.0f ms)\n", 
                             oxygenTime/1e9, oxygenCalls, oxygenTime/1e6/oxygenCalls);
            System.out.printf("  PK:        %.3f s (%d calls, avg %.0f ms)\n", 
                             pkTime/1e9, pkCalls, pkTime/1e6/pkCalls);
            System.out.printf("  RadioBio:  %.3f s\n", radioBioTime/1e9);
            System.out.printf("  CellIter:  %.3f s\n", cellIterTime/1e9);
            System.out.printf("  Clean:     %.3f s\n", cleanTime/1e9);
            System.out.printf("  TOTAL:     %.3f s\n", totalStepTime/1e9);
            
            // Check for growing data structures
            System.out.printf("  PK states: %d, DoseRates: %d, Pops: %d\n",
                             PKStateVariables.size(), DoseRateList.size(), PopsOverTime.size());
        }
    }
	
	
    /**
     * Update PK model with current tumor geometry
     * Called every hour to keep PK synchronized with tumor growth
     * 
     * Computes:
     * - Average tumor radius (for oxygen solver)
     * - Total receptor count (for PK binding calculations)
     * - Updates PK model with new geometry
     * 
     * @return Average tumor radius in meters
     */
    private double updatePKGeometry() {
        // Count tumor cells (excluding vessels)
        int totalCellCount = (int)this.Pop() - currentVesselCounts;
        
        // Estimate average tumor radius from 2D cell count
        // Assuming cells are in a circle: N = π·r²
        // r_average is in units of cell lengths (grid spacing)
        double r_average_cells = Math.sqrt(totalCellCount / Math.PI);
        
        // Convert to meters
        double r_average_m = r_average_cells * SimParams.CELL_LENGTH;  // m
        
        // Store for access by other components (in cell lengths for compatibility)
        this.tempAverageRad = r_average_cells;

        
        // Update PK model with current tumor geometry
        // This updates receptor count and volumes internally
		double h_tumor_m = 2.0 * r_average_m;  // tumor height for cylindrical extrusion
        double tumorVolume = SimParams.computeTumorVolume(totalCellCount);
        double receptorMoles = SimParams.computeReceptorMoles(totalCellCount, currentVesselCounts);        
        int vesselsNearTumor = countVesselsNearTumor();
		this.PBPK.updateGeometry(tumorVolume, receptorMoles, vesselsNearTumor, h_tumor_m);

        return r_average_m;
    }

    /**
     * Seed initial tumor at center of domain
     */
    public void SeedTumor() {
        int[] hood = CircleHood(false, SimParams.INITIAL_TUMOR_RADIUS_CELLS);
        int seedCellCount = MapHood(hood, xDim / 2, yDim / 2);
        
        for (int i = 0; i < seedCellCount; i++) {
            Cell c = GetAgent(hood[i]);
            if (c == null) {
                NewAgentSQ(hood[i]).Init(SimParams.NORMAL, 0);
            }
        }
    }

    /**
     * Generate vessels from configuration file
     * Uses vessel density setting from SimParams to select appropriate CSV
     */
    public int GenVessels() throws IOException {
        // Get vessel configuration file path based on density setting
        String vesselConfigPath = getVesselConfigPath();
        
        // Generate configuration if it doesn't exist
        this.GenerateVesselConf(vesselConfigPath);
        
        // Load vessel positions from CSV
        VesselConfigConvertor convertor = new VesselConfigConvertor(vesselConfigPath);
        double[][] grid = convertor.grid;
        int vesselCount = 0;

        for (int i = 0; i < this.length; i++) {
            int x = ItoX(i);
            int y = ItoY(i);
            
            if (grid[x][y] == 1) {
                NewAgentSQ(i).Init(SimParams.VESSEL, 0);
                this.vesselsIndex.add(i);
                ++vesselCount;
            }
        }

        return vesselCount;
    }

    /**
     * Get vessel configuration file path based on density setting
     * 
     * @return Path to vessel CSV file (resource path format)
     */
    private String getVesselConfigPath() {
        String density = SimParams.VESSEL_DENSITY_CONFIG;
        return String.format("vasculature/%s.csv", density);
    }

    /**
     * Generate vessel configuration file if needed
     * 
     * Currently checks if file exists in resources.
     * Future: Call Python script to generate custom vessel patterns
     * 
     * @param configPath Path to check/generate
     */
    private void GenerateVesselConf(String configPath) {
        // Check if resource exists
        InputStream testStream = getClass().getClassLoader().getResourceAsStream(configPath);
        
        if (testStream != null) {
            // File exists, close and return
            try {
                testStream.close();
            } catch (IOException e) {
                // Ignore
            }
            return;
        }
        
        // TODO: File doesn't exist - call Python script to generate it
        System.err.println("WARNING: Vessel configuration not found: " + configPath);
        System.err.println("Expected location: src/main/resources/" + configPath);
        System.err.println("Please ensure the file exists or implement GenerateVesselConf() to create it.");
    }

    /**
     * Count tumor cells (excluding vessels and healthy tissue)
     * Used by PK model for receptor calculation
     */
    public int countTumorCells() {
        int count = 0;
        for (Cell cell : this) {
            if (cell != null && 
                (cell.type == SimParams.NORMAL || 
                 cell.type == SimParams.HYPOXIC ||
                 cell.type == SimParams.NECROTIC ||
                 cell.type == SimParams.APOPTOTIC)) {
                count++;
            }
        }
        return count;
    }

    /**
     * Count vessels in domain
     * Used by PK model for flow/PS calculation
     */
    public int countVessels() {
        return vesselsIndex.size();
    }
	
	/**
	 * Count vessels near tumor cells (within influence radius)
	 * Handles fragmented tumors by checking neighborhood of each tumor cell
	 * 
	 * @return Number of unique vessels adjacent to tumor
	 */
	public int countVesselsNearTumor() {
		Set<Integer> vesselsNearTumor = new HashSet<>();
		int[] neighborhood = CircleHood(true, SimParams.VESSEL_INFLUENCE_RADIUS);
		
		for (Cell cell : this) {
			if (cell != null && isTumorCell(cell)) {
				int hoodSize = MapHood(neighborhood, cell.Isq());
				for (int i = 0; i < hoodSize; i++) {
					Cell neighbor = GetAgent(neighborhood[i]);
					if (neighbor != null && neighbor.type == SimParams.VESSEL) {
						vesselsNearTumor.add(neighborhood[i]);
					}
				}
			}
		}
		return vesselsNearTumor.size();
	}
	
	/**
	 * Check if cell is a living tumor cell
	 */
	private boolean isTumorCell(Cell c) {
		return c.type == SimParams.NORMAL || c.type == SimParams.HYPOXIC;
	}

    /**
     * Utility: copy array
     */
    private double[] copyArray(double[] original) {
        return Arrays.copyOf(original, original.length);
    }

    // ===================================================================
    // DIAGNOSTICS AND REPORTING
    // ===================================================================
/**
 * Print diagnostic information about current simulation state
 */
	public void printDiagnostics(int currentDay) {
		System.out.println("\n=== DAY " + currentDay + " STATUS ===");
		
		// Vessel status
		int blocked = 0;
		for (int idx : vesselsIndex) {
			if (GetAgent(idx).blockedVessel) blocked++;
		}
		System.out.printf("Vessels: %d total, %d blocked (%.1f%%)%n", 
						 vesselsIndex.size(), blocked, 
						 100.0*blocked/vesselsIndex.size());
		
		// Cell populations
		System.out.printf("Populations: Normal=%d, Hypoxic=%d, Necrotic=%d, Apoptotic=%d%n",
						 (int)CurrentCellsPops[SimParams.NORMAL],
						 (int)CurrentCellsPops[SimParams.HYPOXIC],
						 (int)CurrentCellsPops[SimParams.NECROTIC],
						 (int)CurrentCellsPops[SimParams.APOPTOTIC]);
		
		// PK state (if relevant)
		if (PKStateVariables.size() > 0) {
			double[] pk = PKStateVariables.get(PKStateVariables.size() - 1);
			double N_cen_total = (pk[0] + pk[1]) * 1e9;  // nmol
			double N_ic_total = (pk[8] + pk[9]) * 1e9;
			System.out.printf("PK: N_cen=%.2e nmol, N_ic=%.2e nmol%n",
							 N_cen_total, N_ic_total);
		}
		
		// Oxygen statistics
		double oxygenMin = Double.MAX_VALUE;
		double oxygenMax = 0;
		double oxygenSum = 0;
		int oxygenCount = 0;
		
		for (Cell cell : this) {
			if (cell != null && cell.type != SimParams.VESSEL) {
				double o2 = oxygenGrid.Get(cell.Isq());
				oxygenMin = Math.min(oxygenMin, o2);
				oxygenMax = Math.max(oxygenMax, o2);
				oxygenSum += o2;
				oxygenCount++;
			}
		}
		
		if (oxygenCount > 0) {
			double oxygenAvg = oxygenSum / oxygenCount;
			System.out.printf("Oxygen: min=%.2e, avg=%.2e, max=%.2e (Pa)%n",
							 oxygenMin, oxygenAvg, oxygenMax);
		}
	}


}
