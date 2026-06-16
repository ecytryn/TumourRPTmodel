package TumourRPT;

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
//        System.out.printf("Vessels: %d total, %d blocked (%.1f%%)%n", 
//                         vesselsIndex.size(), blocked, 100.0*blocked/vesselsIndex.size());

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
		
		long oxygenTime = 0, pkTime = 0, radioBioTime = 0, cellIterTime = 0;
		int oxygenCalls = 0, pkCalls = 0;
		
		// ===== GUARANTEED 24 HOURLY UPDATES =====
		for (int hourCount = 0; hourCount < 24; hourCount++) {
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
			
			// Time radiobiology - update all hour states
			long t3 = System.nanoTime();
			this.radioBio.updateAllStates();
			radioBioTime += System.nanoTime() - t3;
			
			PopsOverTime.add(copyArray(CurrentCellsPops));
		
			// ===== CELL UPDATES =====
			for (Cell cell : this) {
				long t4 = System.nanoTime();
				if (!SimParams.FREEZE_TUMOR) {
					cell.Step(hourCount, currentDay);  
				}
				cellIterTime += System.nanoTime() - t4;
			}
			
			CleanAgents();
			ShuffleAgents(this.rng);
		}

		
		// Timing report (if enabled)
		if (SimParams.VERBOSE_ON && currentDay % 10 == 0) {
			long stepTime = System.nanoTime() - stepStart;
			System.out.printf("Day %d timing:\n", currentDay);
			System.out.printf("  Oxygen:    %.3f s (%d calls)\n", oxygenTime/1e9, oxygenCalls);
			System.out.printf("  PK:        %.3f s (%d calls)\n", pkTime/1e9, pkCalls);
			System.out.printf("  RadioBio:  %.3f s\n", radioBioTime/1e9);
			System.out.printf("  Cells:     %.3f s\n", cellIterTime/1e9);
			System.out.printf("  Total:     %.3f s\n", stepTime/1e9);
			System.out.printf("  PKStateVariables.size=%d, DoseRateList.size=%d\n",
							 PKStateVariables.size(), DoseRateList.size());
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
//		if (SimParams.globalTime == 0) {
//			System.out.printf("Vessels near tumor (for PK): %d%n", vesselsNearTumor);
//		}
		this.PBPK.updateGeometry(tumorVolume, receptorMoles, vesselsNearTumor, h_tumor_m);

        return r_average_m;
    }

    /**
     * Seed initial tumor at center of domain
     */
	public void SeedTumor() {
		double r_m = SimParams.INITIAL_TUMOR_RADIUS;
		double r_cells = r_m / SimParams.CELL_LENGTH;
		
		int cx = xDim / 2;
		int cy = yDim / 2;

		// For small radii, converting "a ball of that radius" into a set of cells seems 
		// to be sensitive procedure for HAL, generating inconsistent results for unclear 
		// reasons. Possible large radii too but harder to see that in the results. So 
		// for radii below 2*sqrt(2), we do it manually. Recall that the cells designated 
		// here might not end up in the actual initial configuration if there's already a 
		// vessel in that node.
		
		// Thresholds in cell-length units
		double R1      = 1.0;
		double R_SQRT2 = Math.sqrt(2.0);
		double R2      = 2.0;
		double R_SQRT5 = Math.sqrt(5.0);
		double R_2SQRT2 = 2.0 * Math.sqrt(2.0);
		
		if (r_cells < R_2SQRT2) {
			// Hard-coded discrete shapes, centred on (cx, cy)
			// Offsets from centre to include, in (dx, dy) pairs
			int[][] offsets;
			
			if (r_cells < R1) {
				// 1 cell: centre only
				offsets = new int[][]{{0, 0}};
			} else if (r_cells < R_SQRT2) {
				// 5 cells: centre + 4 cardinal neighbours
				offsets = new int[][]{
					{0, 0},
					{1, 0}, {-1, 0}, {0, 1}, {0, -1}
				};
			} else if (r_cells < R2) {
				// 9 cells: 3x3 square
				offsets = new int[][]{
					{-1,-1}, {0,-1}, {1,-1},
					{-1, 0}, {0, 0}, {1, 0},
					{-1, 1}, {0, 1}, {1, 1}
				};
			} else if (r_cells < R_SQRT5) {
				// 13 cells: 3x3 + 4 cardinal extensions
				offsets = new int[][]{
					{-1,-1}, {0,-1}, {1,-1},
					{-1, 0}, {0, 0}, {1, 0},
					{-1, 1}, {0, 1}, {1, 1},
					{0, -2}, {0, 2}, {-2, 0}, {2, 0}
				};
			} else {
				// 21 cells: 5x5 minus corners
				offsets = new int[][]{
							  {0,-2}, {1,-2}, {-1,-2},
					{-2,-1}, {-1,-1}, {0,-1}, {1,-1}, {2,-1},
					{-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0},
					{-2, 1}, {-1, 1}, {0, 1}, {1, 1}, {2, 1},
							  {0, 2}, {1, 2}, {-1, 2}
				};
			}
			
			for (int[] offset : offsets) {
				int idx = I(cx + offset[0], cy + offset[1]);
				Cell c = GetAgent(idx);
				if (c == null) {
					NewAgentSQ(idx).Init(SimParams.NORMAL, 0);
				}
			}
						
		} else {
			// Large tumour: use CircleHood with int radius
			int r_int = (int) Math.round(r_cells);
			int[] hood = CircleHood(true, r_int);
			int seedCellCount = MapHood(hood, cx, cy);
			for (int i = 0; i < seedCellCount; i++) {
				Cell c = GetAgent(hood[i]);
				if (c == null) {
					NewAgentSQ(hood[i]).Init(SimParams.NORMAL, 0);
				}
			}
		}
	}

	/**
	 * Develop hypoxia in tumor WITHOUT growth
	 * Strategy: Freeze cell division/death, but run oxygen solver and cell type transitions
	 * This allows vessels to occlude and hypoxia to develop based on existing tumor size
	 * 
	 * @param burnInDays Number of days to develop hypoxia (without growth)
	 * @param verbose Print progress updates
	 */
	public void developHypoxiaWithoutGrowth(int burnInDays, boolean verbose) {
		if (verbose) {
			System.out.println("\n=== HYPOXIA DEVELOPMENT PHASE (NO GROWTH) ===");
			System.out.println("Developing hypoxia for " + burnInDays + " days...");
			System.out.println("Tumor size is FROZEN - only oxygen/hypoxia dynamics active");
			
			// Initial state
			int normoxic = (int) CurrentCellsPops[SimParams.NORMAL];
			int hypoxic = (int) CurrentCellsPops[SimParams.HYPOXIC];
			int total = normoxic + hypoxic;
			System.out.printf("Day 0: %d cells (%.1f%% hypoxic)%n", 
							 total, total > 0 ? 100.0 * hypoxic / total : 0);
		}
		
		// Freeze tumor dynamics
		boolean originalFreezeSetting = SimParams.FREEZE_TUMOR;
		SimParams.FREEZE_TUMOR = true;
		
		// Run simulation for burn-in period
		// Even though FREEZE_TUMOR=true, we still need to:
		// 1. Update oxygen field (vessels can occlude, oxygen changes)
		// 2. Update cell types (normoxic -> hypoxic transitions)
		// 3. NOT update PK (no drug yet)
		
		for (int day = 0; day < burnInDays; day++) {
			// Manual hour loop (similar to Step() but without PK)
			for (int hourCount = 0; hourCount < 24; hourCount++) {
				SimParams.updateGlobalTime(day, hourCount);
				
				// Update oxygen field
				double r_average = updatePKGeometry();
				this.oxygen.UpdateSteadyStateOxygen(day, r_average);
				
				// Update cell oxygen levels and allow type transitions
				for (Cell cell : this) {
					if (cell != null && cell.type != SimParams.VESSEL) {
						// Update oxygen level
						cell.oxygen = this.oxygenGrid.Get(cell.Isq());
						
						// Allow normoxic <-> hypoxic transitions based on oxygen
						if (cell.type == SimParams.NORMAL && 
							cell.oxygen < SimParams.P_O2_HYPOXIC) {
							// Become hypoxic
							cell.ChangeType(SimParams.HYPOXIC);
						} else if (cell.type == SimParams.HYPOXIC && 
								   cell.oxygen >= SimParams.P_O2_HYPOXIC) {
							// Reoxygenate
							cell.ChangeType(SimParams.NORMAL);
						} else if (cell.oxygen < SimParams.P_O2_NECROTIC &&
								   cell.type != SimParams.NECROTIC &&
								   cell.type != SimParams.APOPTOTIC) {
							// Die from extreme hypoxia
							cell.ChangeType(SimParams.NECROTIC);
						}
					}
				}
				
				// Log populations (for tracking)
//				PopsOverTime.add(copyArray(CurrentCellsPops));
			}
			
			// Progress report
			if (verbose && ((day + 1) % 5 == 0 || day == burnInDays - 1)) {
				int normoxic = (int) CurrentCellsPops[SimParams.NORMAL];
				int hypoxic = (int) CurrentCellsPops[SimParams.HYPOXIC];
				int necrotic = (int) CurrentCellsPops[SimParams.NECROTIC];
				int total = normoxic + hypoxic;
				double hypoxicFrac = total > 0 ? (double)hypoxic / total : 0;
				
				// Count occluded vessels
				int occluded = 0;
				for (int idx : vesselsIndex) {
					if (GetAgent(idx).blockedVessel) occluded++;
				}
				
				System.out.printf("Day %d: %d viable cells (%.1f%% hypoxic, %d necrotic), %d/%d vessels occluded%n", 
								 day + 1, total, hypoxicFrac * 100, necrotic,
								 occluded, vesselsIndex.size());
			}
		}
		
		// Restore original freeze setting
		SimParams.FREEZE_TUMOR = originalFreezeSetting;
		
		if (verbose) {
			// Final summary
			int normoxic = (int) CurrentCellsPops[SimParams.NORMAL];
			int hypoxic = (int) CurrentCellsPops[SimParams.HYPOXIC];
			int necrotic = (int) CurrentCellsPops[SimParams.NECROTIC];
			int total = normoxic + hypoxic;
			double hypoxicFrac = total > 0 ? (double)hypoxic / total : 0;
			
			int occluded = 0;
			for (int idx : vesselsIndex) {
				if (GetAgent(idx).blockedVessel) occluded++;
			}
			
			System.out.println("\n=== HYPOXIA DEVELOPMENT COMPLETE ===");
			System.out.printf("Final state after %d days (NO GROWTH):%n", burnInDays);
			System.out.printf("  Viable tumor cells: %d (UNCHANGED)%n", total);
			System.out.printf("  Normoxic: %d (%.1f%%)%n", normoxic, 100.0 * normoxic / total);
			System.out.printf("  Hypoxic: %d (%.1f%%)%n", hypoxic, hypoxicFrac * 100);
			System.out.printf("  Necrotic: %d%n", necrotic);
			System.out.printf("  Occluded vessels: %d/%d (%.1f%%)%n", 
							 occluded, vesselsIndex.size(), 
							 100.0 * occluded / vesselsIndex.size());
			System.out.println("Starting treatment phase...\n");
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
        return String.format("vasculature/Capillaries_Density%s.csv", density);
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
        
        // Requested vessel file doesn't exist. Prompt the user to create it.
        System.err.println("WARNING: Vessel configuration not found: " + configPath);
        System.err.println("Expected location: src/main/resources/" + configPath);
        System.err.println("Please use GenerateUniformVessels.py to create it.");
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
