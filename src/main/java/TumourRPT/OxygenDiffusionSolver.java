package TumorRPT;

import HAL.GridsAndAgents.Grid2Ddouble;
import java.awt.FontMetrics;

/**
 * Solves the steady-state oxygen diffusion-consumption equation:
 * D∇²u - C(x,y)u + aV(x,y)(u_B - u) = 0
 * 
 * Uses Red-Black Successive Over-Relaxation (SOR) for efficient iterative solution.
 * Suitable for large grids (400x400) with hourly updates.
 */
public class OxygenDiffusionSolver {
    
    private Grid grid;
    private int xDim, yDim;
    private double h; // grid spacing 
    
    // PDE parameters
    private double D; // diffusion coefficient (converted to grid units²/time)
    private double a_vessel2tissue; // vessel-tissue transfer rate (1/time)
    private double u_B; // blood oxygen concentration
    
    // Solver parameters
    private double omega = 1.8; // SOR relaxation parameter (1.0 = Gauss-Seidel, ~1.8-1.95 optimal)
    private int maxIterations = 5000;
    private double tolerance = 1e-6;
    
    // Whether to enforce vessels as hard Dirichlet BC or let them be determined by the PDE
    private boolean vesselsAsDirichlet = false; // Set to false for physically correct behavior
    
    // Calibrated boundary condition (set by calibrateBoundaryCondition method)
    private double calibratedBoundary = -1; // -1 means not yet calibrated

	// THIS IS NOT A CONTROL ON USING A CALIBRATED BC FOR THE PDE. 
	// IT IS A DEFAULT VALUE - DO NOT CHANGE IT. THE BOOLEAN TO CHANGE IS IN SymParams.
    private boolean useCalibrated = false;
    
    // Diagnostic flags
    private int consumptionBuildCount = 0;
    
    // Day on which to add a scale bar to the oxygen image
    private int ScaleBarDay = 55;
    
	// Boundary condition type: true = Neumann (zero flux), false = Dirichlet (fixed value)
	// NOTE: Neumann BC is NOT IMPLEMENTED - always uses Dirichlet regardless of this flag
	private boolean useNeumannBC = false; // Currently only Dirichlet is implemented
    
    public OxygenDiffusionSolver(Grid grid) {
        this.grid = grid;
        this.xDim = grid.xDim;
        this.yDim = grid.yDim;
		this.h = SimParams.CELL_LENGTH;
        
        // Initialize PDE parameters from Params
        this.D = SimParams.D_O2; // Already in correct units from Params
        this.u_B = SimParams.P_O2_VESSEL;
        
        // Estimate vessel-tissue transfer rate 'a_vessel2tissue'
        // The CELLS_CONSUMPTION_RATE_LIST values are already in simulation units (1/T = 1/hour)
        // From observed decay length ~30 um: k = sqrt(C/D) ≈ 1/30 um^-1
        // This gives C ≈ 7500 hr^-1 for normal cells, matching the code values
        // 
        // To approximate the old model's behavior (vessels as perfect sources),
        // set a >> C so vessel term dominates near vessels
        // Start with a = 10*C_typical
        double C_typical = SimParams.CELLS_CONSUMPTION_RATE_LIST[SimParams.NORMAL];
        this.a_vessel2tissue = 4.0 * C_typical; // Should give a_vessel2tissue ~ 75000 hr^-1

//        System.out.printf("Oxygen solver init: a=%.3e 1/s, u_B=%.3e Pa (%.1f mmHg)\n",
//                         this.a_vessel2tissue, u_B, u_B/SimParams.MMHG_TO_PA);
        
        // Verify D is in correct units (should be um^2/hour)
//        double D_check = this.D; // Already converted by Params
        
//        System.out.printf("Oxygen solver initialized: D=%.3e um^2/hr, a_vessel2tissue=%.3e 1/hr, u_B=%.3e nmol/um^3\n", D, a, u_B);
//        System.out.printf("Expected decay length: lambda ~ %.1f um (for C=%.0f 1/hr)\n", 
//                         Math.sqrt(D/C_typical), C_typical);
    }
    
    /**
     * Solves the steady-state oxygen equation using Red-Black SOR
     * 
     * @param oxygenGrid Grid to store solution (modified in place)
     * @param currentDay Current simulation day
     * @return Number of iterations to convergence
     */
    public int solve(Grid2Ddouble oxygenGrid, int currentDay) {
        
        // Diagnostic: Check oxygen BEFORE solving
        if (currentDay == 0) {
            double sum = 0, min = Double.MAX_VALUE, max = 0;
            for (int i = 0; i < oxygenGrid.length; i++) {
                double val = oxygenGrid.Get(i);
                sum += val;
                min = Math.min(min, val);
                max = Math.max(max, val);
            }
			if (SimParams.VERBOSE_ON) {
				System.out.printf("BEFORE solve: min=%.3e, avg=%.3e, max=%.3e Pa\n", 
								 min, sum/oxygenGrid.length, max);
			}
        }
        
        // Build consumption rate field C(x,y)
        Grid2Ddouble consumptionField = buildConsumptionField();
        
        // Build vessel indicator field V(x,y)
        Grid2Ddouble vesselField = buildVesselField();
        
        // Use previous solution as initial guess (warm start)
        // oxygenGrid already contains previous solution
        
        // Boundary conditions: Set boundary to healthy tissue oxygen level
		// TODO: Implement Neumann BC option - currently always uses Dirichlet
		double u_boundary;
		if (useCalibrated && calibratedBoundary > 0) {
			// use average from a tumour-free simulation carried out during initialization 
			// to get healthy oxygen level for current parameter values
			u_boundary = calibratedBoundary;  // Use calibrated value
		} else {
            // Use intermediate value between arterial and venous (healthy tissue average)
			u_boundary = 0.5 * (SimParams.P_O2_VESSEL + SimParams.P_O2_VEIN);
		}
		applyBoundaryConditions(oxygenGrid, u_boundary);
        
        // Red-Black SOR iterations
        int iter;
        double maxResidual = Double.MAX_VALUE;
        
        for (iter = 0; iter < maxIterations && maxResidual > tolerance; iter++) {
            // Red sweep (checkerboard pattern: i+j even)
            maxResidual = sorSweep(oxygenGrid, consumptionField, vesselField, true);
            
            // Black sweep (checkerboard pattern: i+j odd)
            maxResidual = Math.max(maxResidual, 
                                   sorSweep(oxygenGrid, consumptionField, vesselField, false));
            
            // Check convergence every 50 iterations to save time
            if (iter % 50 == 0 && iter > 0) {
                if (maxResidual < tolerance) break;
            }
        }
        
        // Optional: Print convergence info occasionally
        if (iter >= maxIterations || (currentDay < 5 && iter % 100 == 0)) {
//            System.out.printf("Oxygen solver: day %d, iter %d, residual %.3e\n", 
//                            currentDay, iter, maxResidual);
        }
        
        // Diagnostic: Check oxygen values on day 0
        if (currentDay == 0 && SimParams.VERBOSE_ON) {
            System.out.printf("\n=== OXYGEN DEBUG (Day 0) ===\n");
            System.out.printf("u_B (vessel source) = %.3e Pa (%.1f mmHg)\n", 
                             u_B, u_B/SimParams.MMHG_TO_PA);
            System.out.printf("Parameter a (vessel transfer rate) = %.3e 1/s\n", a_vessel2tissue);
            System.out.printf("Converged in %d iterations\n", iter);
            
            // Sample a few oxygen values
            double sum = 0, min = Double.MAX_VALUE, max = 0;
            int count = 0;
            for (int i = 0; i < oxygenGrid.length; i++) {
                double val = oxygenGrid.Get(i);
                sum += val;
                min = Math.min(min, val);
                max = Math.max(max, val);
                count++;
            }
            System.out.printf("O2 after solve: min=%.3e, avg=%.3e, max=%.3e Pa\n", min, sum/count, max);
            System.out.printf("In mmHg: min=%.1f, avg=%.1f, max=%.1f\n", 
                             min/SimParams.MMHG_TO_PA, (sum/count)/SimParams.MMHG_TO_PA, max/SimParams.MMHG_TO_PA);
            System.out.printf("============================\n\n");
        }
        
        // Save oxygen field visualization every N days
		int currentHour = SimParams.globalTime % 24;
        if ((currentDay % 5 == 0 || (currentDay >= 44 && currentDay <= 50)) && currentHour == 0) {
            saveOxygenFieldImage(oxygenGrid, currentDay);
        }
        
        return iter;
    }
    
    /**
     * Performs one Red-Black SOR sweep
     * 
     * @param u Current solution
     * @param C Consumption field
     * @param V Vessel indicator field
     * @param isRed true for red sweep (i+j even), false for black (i+j odd)
     * @return Maximum residual in this sweep
     */
    private double sorSweep(Grid2Ddouble u, Grid2Ddouble C, Grid2Ddouble V, boolean isRed) {
        double maxResidual = 0.0;
        double h2 = h * h;
        
        // Loop over interior points with checkerboard pattern
        for (int x = 1; x < xDim - 1; x++) {
            for (int y = 1; y < yDim - 1; y++) {
                
                // Skip if wrong color in checkerboard
                if (((x + y) % 2 == 0) != isRed) continue;
                
                int idx = grid.I(x, y);
                
                // If treating vessels as hard Dirichlet BC (optional)
                if (vesselsAsDirichlet && V.Get(idx) > 0.5) {
                    u.Set(idx, u_B);
                    continue;
                }
                
                // Get neighboring values
                double u_xp = u.Get(grid.I(x+1, y));
                double u_xm = u.Get(grid.I(x-1, y));
                double u_yp = u.Get(grid.I(x, y+1));
                double u_ym = u.Get(grid.I(x, y-1));
                double u_old = u.Get(idx);
                
                double C_val = C.Get(idx);
                double V_val = V.Get(idx);
                
                // Discretized equation (5-point stencil):
                // D(u_xp + u_xm + u_yp + u_ym - 4u)/h² - C*u + aV(u_B - u) = 0
                //
                // Rearranging for u:
                // [4D/h² + C + aV]*u = D(u_xp + u_xm + u_yp + u_ym)/h² + aV*u_B
                // u = [D(u_xp + u_xm + u_yp + u_ym)/h² + aV*u_B] / [4D/h² + C + aV]
                
                double numerator = D * (u_xp + u_xm + u_yp + u_ym) / h2 + a_vessel2tissue * V_val * u_B;
                double denominator = 4.0 * D / h2 + C_val + a_vessel2tissue * V_val;
                
                double u_new = numerator / denominator;
                
                // SOR update: u = u_old + omega * (u_new - u_old)
                u_new = u_old + omega * (u_new - u_old);
                
                // Ensure non-negative oxygen
                u_new = Math.max(0.0, u_new);
                
                u.Set(idx, u_new);
                
                // Track maximum residual for convergence check
                double residual = Math.abs(u_new - u_old);
                maxResidual = Math.max(maxResidual, residual);
            }
        }
        
        return maxResidual;
    }
    
    /**
     * Builds the consumption rate field C(x,y) based on cell types
     */
    private Grid2Ddouble buildConsumptionField() {
        Grid2Ddouble C = new Grid2Ddouble(xDim, yDim);
        
        // Track cell type counts for diagnostics
        int[] cellTypeCounts = new int[SimParams.NUM_CELL_TYPES];
        
        for (int i = 0; i < grid.length; i++) {
            Cell cell = grid.GetAgent(i);
            
            int cellType;
            if (cell == null) {
                cellType = SimParams.HEALTHY;
            } else if (cell.type == SimParams.VESSEL) {
                // Vessels don't consume oxygen (or very little)
                cellType = SimParams.HEALTHY;
            } else {
                cellType = cell.type;
                cellTypeCounts[cell.type]++;
            }
            
            // Get consumption rate for this cell type
            double consumptionRate = SimParams.CELLS_CONSUMPTION_RATE_LIST[cellType];
            C.Set(i, consumptionRate);
        }
        
        // Diagnostic: Print first few builds and save visualization
        consumptionBuildCount++;
        if (consumptionBuildCount <= 3) {
	        if (SimParams.VERBOSE_ON) {
				System.out.printf("\n=== CONSUMPTION FIELD BUILD #%d ===\n", consumptionBuildCount);
				System.out.printf("Cell type counts:\n");
				System.out.printf("  Healthy: %d (C=%.3e)\n", 
								 xDim*yDim - cellTypeCounts[SimParams.NORMAL] - cellTypeCounts[SimParams.HYPOXIC] 
								 - cellTypeCounts[SimParams.NECROTIC] - cellTypeCounts[SimParams.APOPTOTIC] 
								 - cellTypeCounts[SimParams.VESSEL],
								 SimParams.CONSUMPTION_HEALTHY);
				System.out.printf("  Normal:  %d (C=%.3e)\n", cellTypeCounts[SimParams.NORMAL], SimParams.CONSUMPTION_NORMAL);
				System.out.printf("  Hypoxic: %d (C=%.3e)\n", cellTypeCounts[SimParams.HYPOXIC], SimParams.CONSUMPTION_HYPOXIC);
				System.out.printf("  Necrotic:%d (C=%.3e)\n", cellTypeCounts[SimParams.NECROTIC], SimParams.CONSUMPTION_NECROTIC);
				System.out.printf("  Apoptotic:%d (C=%.3e)\n", cellTypeCounts[SimParams.APOPTOTIC], SimParams.CONSUMPTION_APOPTOTIC);
				System.out.printf("  Vessel:  %d (C=%.3e)\n", cellTypeCounts[SimParams.VESSEL], SimParams.CONSUMPTION_VESSEL);
			}				
            double sum = 0, min = Double.MAX_VALUE, max = 0;
            for (int i = 0; i < C.length; i++) {
                double val = C.Get(i);
                sum += val;
                min = Math.min(min, val);
                max = Math.max(max, val);
            }
			if (SimParams.VERBOSE_ON) {
				System.out.printf("Consumption field stats: min=%.3e, avg=%.3e, max=%.3e 1/s\n",
								 min, sum/C.length, max);
			}            
            // Save visualization
            saveConsumptionFieldImage(C, consumptionBuildCount);
            
			if (SimParams.VERBOSE_ON) {
				System.out.printf("=====================================\n\n");
            }
        }
        
        return C;
    }
    
    /**
     * Save consumption field as a PNG heatmap with fixed color scale
     * Uses fixed scale: 0 to 2.5 1/s (to capture healthy=0.42, normal=2.08)
     */
    private void saveConsumptionFieldImage(Grid2Ddouble C, int buildNumber) {
        if (SimParams.EXPORT_CONSUMP_IMAGES) {
            // Only save if visualization is enabled
            try {
                String filename = String.format("results/single_runs/consumption_field_%d.png", buildNumber);
                
                // Create buffered image with extra space for color bar
                int imageWidth = xDim + 100;  // Add 100 pixels on right for color bar
                int imageHeight = yDim;
                java.awt.image.BufferedImage img = new java.awt.image.BufferedImage(
                    imageWidth, imageHeight, java.awt.image.BufferedImage.TYPE_INT_RGB);
                
                // Fixed color scale: 0 to 2.5 1/s
                double minScale = 0.0;
                double maxScale = 2.5;  // Covers healthy (0.42) to normal (2.08)
                
                // Draw main heatmap
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < yDim; y++) {
                        int idx = grid.I(x, y);
                        double val = C.Get(idx);
                        
                        // Normalize to [0, 1]
                        double normalized = (val - minScale) / (maxScale - minScale);
                        normalized = Math.max(0.0, Math.min(1.0, normalized));
                        
                        // Use blue (low) to red (high) colormap
                        int color = HAL.Util.HeatMapRGB(normalized, 0, 1);
                        img.setRGB(x, yDim - 1 - y, color);  // Flip y for image coordinates
                    }
                }
                
                // Draw color bar
                int barX = xDim + 20;  // Start of color bar
                int barWidth = 30;
                int barY = 50;
                int barHeight = yDim - 100;
                
                for (int i = 0; i < barHeight; i++) {
                    double normalized = 1.0 - (double)i / barHeight;  // Top = high, bottom = low
                    int color = HAL.Util.HeatMapRGB(normalized, 0, 1);
                    for (int j = 0; j < barWidth; j++) {
                        img.setRGB(barX + j, barY + i, color);
                    }
                }
                
                // Add text annotations
                java.awt.Graphics2D g2d = img.createGraphics();
                g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, 
                                    java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 12*3/2));
                
                // Title
                g2d.setColor(java.awt.Color.WHITE);
                String title = String.format("Consumption Field #%d", buildNumber);
                g2d.drawString(title, 10, 20);
                g2d.setColor(java.awt.Color.BLACK);
                g2d.drawString(title, 9, 19);
                g2d.drawString(title, 11, 19);
                
                // Color bar labels
                g2d.setColor(java.awt.Color.BLACK);
                g2d.drawRect(barX, barY, barWidth, barHeight);  // Border around color bar
                
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 11));
                
                // Max label (top)
                String maxLabel = String.format("%.2f", maxScale);
                g2d.drawString(maxLabel, barX + barWidth + 5, barY + 5);
                
                // Mid label
                String midLabel = String.format("%.2f", (minScale + maxScale) / 2);
                g2d.drawString(midLabel, barX + barWidth + 5, barY + barHeight / 2);
                
                // Min label (bottom)
                String minLabel = String.format("%.2f", minScale);
                g2d.drawString(minLabel, barX + barWidth + 5, barY + barHeight - 5);
                
                // Units label
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 10));
                g2d.drawString("1/s", barX + barWidth + 5, barY + barHeight + 15);
                
                // Add key consumption rates as reference
                g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 10));
                int refY = barY + barHeight + 40;
                g2d.drawString("Healthy: 0.42", barX - 10, refY);
                g2d.drawString("Normal: 2.08", barX - 10, refY + 15);
                g2d.drawString("Necrotic: 0.0003", barX - 20, refY + 30);
                
                g2d.dispose();
                
                // Save image
                javax.imageio.ImageIO.write(img, "png", new java.io.File(filename));
                System.out.printf("Saved consumption field visualization: %s\n", filename);
                
            } catch (Exception e) {
                System.err.printf("Warning: Could not save consumption field image: %s\n", 
                                e.getMessage());
            }
        }
    }
    

	/**
	 * Builds the vessel indicator field V(x,y)
	 * V = vesselAmp at functioning vessel locations
	 * V = 0 at blocked vessels or non-vessel locations
	 */
	private Grid2Ddouble buildVesselField() {
		Grid2Ddouble V = new Grid2Ddouble(xDim, yDim);
		V.SetAll(0.0);
		
		for (int vesselIdx : grid.vesselsIndex) {
			Cell vessel = grid.GetAgent(vesselIdx);
			
			if (vessel != null && !vessel.blockedVessel) {
				// Functioning vessel: V = vesselAmp (accounts for variability)
				V.Set(vesselIdx, vessel.vesselAmp);
			}
			// Blocked vessels stay at V = 0 (no source term)
		}
		
		return V;
	}
    
    /**
     * Applies boundary conditions (Dirichlet with healthy tissue oxygen level)
     */
    private void applyBoundaryConditions(Grid2Ddouble u, double u_boundary) {
        // Set all boundary points to u_boundary
        for (int x = 0; x < xDim; x++) {
            u.Set(grid.I(x, 0), u_boundary);
            u.Set(grid.I(x, yDim-1), u_boundary);
        }
        for (int y = 0; y < yDim; y++) {
            u.Set(grid.I(0, y), u_boundary);
            u.Set(grid.I(xDim-1, y), u_boundary);
        }
    }
    
    // Setters for tuning solver parameters
    public void setOmega(double omega) { this.omega = omega; }
    public void setMaxIterations(int max) { this.maxIterations = max; }
    public void setTolerance(double tol) { this.tolerance = tol; }
    public void setVesselTransferRate(double a_vessel2tissue) { this.a_vessel2tissue = a_vessel2tissue; }
    public void setVesselsAsDirichlet(boolean flag) { this.vesselsAsDirichlet = flag; }
    public void setUseNeumannBC(boolean flag) { this.useNeumannBC = flag; }
    
    /**
     * Calibrates boundary condition by simulating healthy tissue (vessels + healthy cells only)
     * Should be called once at initialization before tumor is seeded
     * 
     * @param oxygenGrid Grid to use for calibration (will be modified)
     * @return Calibrated boundary oxygen pressure
     */
    public double calibrateBoundaryCondition(Grid2Ddouble oxygenGrid) {
//        System.out.println("Calibrating boundary condition from healthy tissue...");
        
        // Use initial guess for BC (won't affect final result much)
        double initialBC = 0.5 * (SimParams.P_O2_VESSEL + SimParams.P_O2_VEIN);
        
        // Initialize oxygen field
        oxygenGrid.SetAll(initialBC);
        applyBoundaryConditions(oxygenGrid, initialBC);
        
        // Build fields (should be all healthy cells + vessels at this point)
        Grid2Ddouble consumptionField = buildConsumptionField();
        Grid2Ddouble vesselField = buildVesselField();
        
        // Solve to steady state with generous iteration budget
        int iter;
        double maxResidual = Double.MAX_VALUE;
        int calibrationMaxIter = 10000; // More iterations for cold start
        
        for (iter = 0; iter < calibrationMaxIter && maxResidual > tolerance; iter++) {
            maxResidual = sorSweep(oxygenGrid, consumptionField, vesselField, true);
            maxResidual = Math.max(maxResidual, 
                                   sorSweep(oxygenGrid, consumptionField, vesselField, false));
            
            if (iter % 100 == 0 && iter > 0) {
                if (maxResidual < tolerance) break;
            }
        }
        
        // Compute statistics on the converged field
        double sum = 0, min = Double.MAX_VALUE, max = 0;
        int count = 0;
        
        // Exclude boundary points from statistics
        for (int x = 1; x < xDim - 1; x++) {
            for (int y = 1; y < yDim - 1; y++) {
                int idx = grid.I(x, y);
                double val = oxygenGrid.Get(idx);
                sum += val;
                min = Math.min(min, val);
                max = Math.max(max, val);
                count++;
            }
        }
        
        double avg = sum / count;
        double variability = (max - min) / avg;
        
//        System.out.printf("Calibration converged in %d iterations (residual %.3e)\n", iter, maxResidual);
//        System.out.printf("Healthy tissue oxygen: avg=%.3e Pa (%.1f mmHg), min=%.3e Pa, max=%.3e Pa\n", 
//                         avg, avg/SimParams.MMHG_TO_PA, min, max);
//        System.out.printf("Relative variability: %.2f (< 0.5 is good)\n", variability);
        
        if (variability > 1.0) {
            System.out.println("WARNING: High oxygen variability - consider increasing parameter 'a'");
        }
        
        // Store calibrated value
        this.calibratedBoundary = avg;
        this.useCalibrated = true;
        
        // Save calibration visualization
        saveOxygenFieldImage(oxygenGrid, -1);  // -1 indicates calibration
        
        return avg;
    }
    
    /**
     * Save oxygen field as a PNG heatmap with fixed color scale
     * Uses fixed scale: 0 to 15000 Pa (0 to ~112 mmHg) to show full range
     */
    private void saveOxygenFieldImage(Grid2Ddouble oxygenGrid, int day) {
        if (SimParams.EXPORT_OX_IMAGES) {
            try {
                String filename;
                if (day == -1) {
                    filename = String.format("%s/oxygen_field_calibration.png",SimParams.OUTPUT_DIR_OXYGEN_IMAGES);
                } else {
                    filename = String.format("%s/oxygen_field_day%d.png", SimParams.OUTPUT_DIR_OXYGEN_IMAGES, day);
                }

				// Export image data to csv to check oxygen bug:
//				DataLogger csvLogger = new DataLogger();
//				String csvFilename = filename.replace(".png", ".csv");
//				csvLogger.saveOxygenFieldCSV(csvFilename, grid);
            
				int imageWidth = xDim;                
				if (day == ScaleBarDay) {
					// Create buffered image with extra space for color bar on day 50 only
					imageWidth = xDim + 150;  // Add 150 pixels on right for color bar
				}
				int imageHeight = yDim;
				java.awt.image.BufferedImage img = new java.awt.image.BufferedImage(
					imageWidth, imageHeight, java.awt.image.BufferedImage.TYPE_INT_RGB);

                // Fill entire image with white background first
                java.awt.Graphics2D g2d_bg = img.createGraphics();
                g2d_bg.setColor(java.awt.Color.WHITE);
                g2d_bg.fillRect(0, 0, imageWidth, imageHeight);
                g2d_bg.dispose();
                                
                // Fixed color scale: 0 to 15000 Pa (0 to ~112 mmHg)
                double minScale = 0.0;
                double maxScale = 15000.0;  // Above vessel pressure
                
                // Draw main heatmap
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < yDim; y++) {
                        int idx = grid.I(x, y);
                        double val = oxygenGrid.Get(idx);
                        
                        // Normalize to [0, 1]
                        double normalized = (val - minScale) / (maxScale - minScale);
                        normalized = Math.max(0.0, Math.min(1.0, normalized));
                        
                        // Use blue (low) to red (high) colormap
                        int color = HAL.Util.HeatMapBGR(normalized, 0, 1);
                        img.setRGB(y, x, color);  // Flip y for image coordinates
                    }
                }
                
				// Draw color bar
				int barX = xDim + 20;  // Start of color bar
				int barWidth = 30;
				int barY = 50;
				int barHeight = yDim - 100;
                
				if (day == ScaleBarDay) {
					for (int i = 0; i < barHeight; i++) {
						double normalized = 1.0 - (double)i / barHeight;  // Top = high, bottom = low
						int color = HAL.Util.HeatMapBGR(normalized, 0, 1);
						for (int j = 0; j < barWidth; j++) {
							img.setRGB(barX + j, barY + i, color);
						}
					}
				}                

                // Add text annotations
                java.awt.Graphics2D g2d = img.createGraphics();
                g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, 
                                    java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
                g2d.setRenderingHint(java.awt.RenderingHints.KEY_TEXT_ANTIALIASING,
                                    java.awt.RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                                                    
                // Title
                String title;
                if (day == -1) {
                    title = "Oxygen Field (Calibration)";
                } else {
                    title = String.format("Day %d", day);
					if ( MyUtils.isElementPresent(SimParams.INJECTION_SCHEDULE, day) ) {
						title = title + " (injection)";
					}
                }
                
                // Use larger, bold font
				int titleSize = (int)(20 * SimParams.FONT_SCALE_FACTOR);
				java.awt.Font titleFont = new java.awt.Font("Arial", java.awt.Font.BOLD, titleSize);
                g2d.setFont(titleFont);

				// Get text dimensions for positioning
				FontMetrics fm = g2d.getFontMetrics();
				int textWidth = fm.stringWidth(title);
				int textHeight = fm.getHeight();
				
				// Position in top-left corner with padding
				int padding = 10;
				int titleX = padding;
				int titleY = padding + fm.getAscent();
                                
                // Draw black outline first for contrast
                g2d.setColor(java.awt.Color.BLACK);
                for (int dx = -2; dx <= 2; dx++) {
                    for (int dy = -2; dy <= 2; dy++) {
                        if (dx != 0 || dy != 0) {
                            g2d.drawString(title, titleX + dx, titleY + dy);
                        }
                    }
                }
                
                // Draw white text on top
                g2d.setColor(java.awt.Color.WHITE);
                g2d.drawString(title, titleX, titleY);
                // ======================================
                
                if (day == ScaleBarDay) {
					// Color bar labels
					g2d.setColor(java.awt.Color.BLACK);
					g2d.drawRect(barX, barY, barWidth, barHeight);  // Border around color bar
	
					int labelSize = (int)(14 * SimParams.FONT_SCALE_FACTOR);                
					g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, labelSize));
					
					// Max label (top) - kPa only
					double maxKPa = maxScale / 1000.0;  // Pa to kPa
					String maxLabel = String.format("%.0f kPa", maxKPa);
					g2d.drawString(maxLabel, barX + barWidth + 5, barY + 10);                
	
					// Mid label
					double midKPa = (minScale + maxScale) / 2000.0;
					String midLabel = String.format("%.1f", midKPa);
					g2d.drawString(midLabel, barX + barWidth + 5, barY + barHeight / 2);
					
					// Min label (bottom)
					String minLabel = "0 kPa";
					g2d.drawString(minLabel, barX + barWidth + 5, barY + barHeight + 3);
					
					// Add threshold markers directly on color bar
					// Calculate y-positions for thresholds on the color bar
					double vesselPos = (SimParams.P_O2_VESSEL - minScale) / (maxScale - minScale);
					double hypoxicPos = (SimParams.P_O2_HYPOXIC - minScale) / (maxScale - minScale);
					double necroticPos = (SimParams.P_O2_NECROTIC - minScale) / (maxScale - minScale);
					
					int vesselY = barY + (int)((1.0 - vesselPos) * barHeight);
					int hypoxicY = barY + (int)((1.0 - hypoxicPos) * barHeight);
					int necroticY = barY + (int)((1.0 - necroticPos) * barHeight);
					
					// Draw threshold lines across color bar (only if within scale)
					g2d.setStroke(new java.awt.BasicStroke(2.0f));
	
					// Threshold labels - simplified, larger font
					int thresholdLabelSize = (int)(11 * SimParams.FONT_SCALE_FACTOR);
					g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, thresholdLabelSize));
					g2d.setColor(java.awt.Color.BLACK);
					
					// Vessel threshold
					if (vesselPos >= 0 && vesselPos <= 1.0) {
						g2d.setColor(java.awt.Color.BLACK);
						g2d.setStroke(new java.awt.BasicStroke(5.0f));
						g2d.drawLine(barX, vesselY, barX + barWidth, vesselY);
						g2d.setColor(new java.awt.Color(1.0f, 0.3f, 0.3f));
//						g2d.setColor(java.awt.Color.GREEN);
						g2d.setStroke(new java.awt.BasicStroke(3.0f));
						g2d.drawLine(barX, vesselY, barX + barWidth, vesselY);
						g2d.setColor(java.awt.Color.BLACK);
						//g2d.drawString("← ", barX + barWidth + 5, vesselY + 4);
					}
					
					// Hypoxic threshold  
					if (hypoxicPos >= 0 && hypoxicPos <= 1.0) {
						g2d.setColor(java.awt.Color.BLACK);
						g2d.setStroke(new java.awt.BasicStroke(5.0f));
						g2d.drawLine(barX, hypoxicY, barX + barWidth, hypoxicY);
						g2d.setColor(java.awt.Color.YELLOW);
						g2d.setStroke(new java.awt.BasicStroke(3.0f));
						g2d.drawLine(barX, hypoxicY, barX + barWidth, hypoxicY);
						g2d.setColor(java.awt.Color.BLACK);
						//g2d.drawString("← ", barX + barWidth + 5, hypoxicY + 4);
					}
					
					// Necrotic threshold
					if (necroticPos >= 0 && necroticPos <= 1.0) {
						g2d.setColor(java.awt.Color.BLACK);
						g2d.setStroke(new java.awt.BasicStroke(5.0f));
						g2d.drawLine(barX, necroticY, barX + barWidth, necroticY);
						g2d.setColor(new java.awt.Color(0.0f, 0.7f, 0.0f));
//						g2d.setColor(java.awt.Color.RED);
						g2d.setStroke(new java.awt.BasicStroke(3.0f));
						g2d.drawLine(barX, necroticY, barX + barWidth, necroticY);
						g2d.setColor(java.awt.Color.BLACK);
						//g2d.drawString("← ", barX + barWidth + 5, necroticY + 4);
					}
									
					// Reset stroke for border redraw
					g2d.setStroke(new java.awt.BasicStroke(1.0f));
					g2d.setColor(java.awt.Color.BLACK);
					g2d.drawRect(barX, barY, barWidth, barHeight);  // Redraw border on top
				}                
                g2d.dispose();
                
                // Save image
                javax.imageio.ImageIO.write(img, "png", new java.io.File(filename));
//                System.out.printf("Saved oxygen field visualization: %s\n", filename);
                
            } catch (Exception e) {
                System.err.printf("Warning: Could not save oxygen field image: %s\n", 
                                e.getMessage());
            }
        }
    }
}
