package TumourRPT;

import HAL.GridsAndAgents.Grid2Ddouble;

import static HAL.Util.CircleHood;

public class Oxygen {
    public Grid grid;
    private OxygenDiffusionSolver diffusionSolver;
    
    // Old method parameters (kept for backward compatibility)
//    int r_max;
    
    // Vessel occlusion parameters (tunable)
    private static final int OCCLUSION_HOOD_RADIUS = 10;  // cells (100 um)
    private static final double BLOCKAGE_PROBABILITY_PER_HOUR = 0.002;  // At full tumor density
    private static final double DENSITY_THRESHOLD = 0.5;  // Min tumor fraction for occlusion
	private static final double REOPENING_THRESHOLD = 0.3;  // Reopen when <30% tumor
	private static final double REOPENING_PROBABILITY_PER_HOUR = 0.005;  // ~50% reopen in 1 week at 0% tumor
    
    public Oxygen(Grid grid) {
        this.grid = grid;
//        this.r_max = SimParams.r_max;
        
        // Initialize new diffusion solver
        this.diffusionSolver = new OxygenDiffusionSolver(grid);
    }
    
    /**
     * Calibrates boundary condition from healthy tissue before tumor seeding
     * Call this once during initialization, after vessels are created but before tumor
     */
    public void calibrateBoundaryCondition() {
        diffusionSolver.calibrateBoundaryCondition(this.grid.oxygenGrid);
    }
    
    /**
     * Main interface - updates oxygen field using selected method
     */
    void UpdateSteadyStateOxygen(int currentDay, double r_average) {
        // Apply pressure-based vessel occlusion first
        applyPressureBasedOcclusion();
		updateWithDiffusionSolver(currentDay);
    }
    
    /**
     * Solve diffusion-consumption PDE
     */
    private void updateWithDiffusionSolver(int currentDay) {
        // Solver modifies oxygenGrid in place
        int iterations = diffusionSolver.solve(this.grid.oxygenGrid, currentDay);
        
        // Optional: warn if solver didn't converge
        if (iterations >= 5000) {
            System.out.println("WARNING: Oxygen solver reached max iterations on day " + currentDay);
        }
    }
    
    /**
     * Applies vessel occlusion based on local tumor burden
     * 
     * Physical model: Tumor cells exert mechanical pressure on nearby vessels.
     * When tumor density is high in a vessel's neighborhood, the vessel has 
     * a probability of being occluded each hour.
     */
    private void applyPressureBasedOcclusion() {
        int[] localHood = CircleHood(true, OCCLUSION_HOOD_RADIUS);
        
        for (int vesselIndex : this.grid.vesselsIndex) {
            
            int hoodSize = this.grid.MapHood(localHood, vesselIndex);
            
            // Count tumor cells in neighborhood
            int tumorCellCount = 0;
            for (int i = 0; i < hoodSize; i++) {
                Cell c = this.grid.GetAgent(localHood[i]);
                if (c != null && isTumorCell(c)) {
                    tumorCellCount++;
                }
            }
            
            double densityFactor = tumorCellCount / (double) hoodSize;

			Cell vessel = this.grid.GetAgent(vesselIndex);

			if (vessel.blockedVessel) {
				// VESSEL IS BLOCKED - check for reopening
				if (densityFactor < REOPENING_THRESHOLD) {
					// Low tumor burden - vessel can reopen
					double reopeningProbability = (REOPENING_THRESHOLD - densityFactor) / REOPENING_THRESHOLD 
												 * REOPENING_PROBABILITY_PER_HOUR;
					
					if (this.grid.rng.Double() < reopeningProbability) {
						vessel.blockedVessel = false;
					}
				}
			} else {
				// VESSEL IS OPEN - check for occlusion
				if (densityFactor > DENSITY_THRESHOLD) {
					double excessDensity = (densityFactor - DENSITY_THRESHOLD) / (1.0 - DENSITY_THRESHOLD);
					double occlusionProbability = excessDensity * BLOCKAGE_PROBABILITY_PER_HOUR;
					
					if (this.grid.rng.Double() < occlusionProbability) {
						vessel.blockedVessel = true;
					}
				}
			}
        }
    }
    
    /**
     * Check if cell is a growing tumor cell (contributes to pressure)
     */
    private boolean isTumorCell(Cell c) {
        return c.type == SimParams.NORMAL || c.type == SimParams.HYPOXIC;
    }
    
    /**
     * Diagnostic: Count and print vessel occlusion statistics
     */
    public void printVesselStats() {
        int totalVessels = this.grid.vesselsIndex.size();
        int blockedCount = 0;
        
        for (int vesselIndex : this.grid.vesselsIndex) {
            if (this.grid.GetAgent(vesselIndex).blockedVessel) {
                blockedCount++;
            }
        }
        
        double blockedPercent = 100.0 * blockedCount / totalVessels;
//        System.out.printf("Vessels: %d total, %d blocked (%.1f%%)\n", 
//                         totalVessels, blockedCount, blockedPercent);
    }
}
