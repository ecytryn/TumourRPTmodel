package TumorRPT;

import HAL.Gui.GridWindow;
import HAL.Gui.PlotLine;
import HAL.Gui.PlotWindow;

import java.util.ArrayList;

import static HAL.Util.*;

/**
 * Visualization manager for tumor simulation
 * 
 * Handles both live interactive display and image export
 * - Live mode (PLOT_LIVE_IMAGES=true): Creates windows for real-time viewing
 * - Headless mode (PLOT_LIVE_IMAGES=false): Uses offscreen buffer for export only
 */
class DaVinci {
    private Grid grid;
    
    // Pixel storage (always present)
    public int xDim;
    public int yDim;
    private int[] pixelBuffer;  // Offscreen buffer for headless mode
    
    // Live display components (null in headless mode)
    public GridWindow gridWin;

    public PlotWindow plotWin;

    public PlotLine testLine;

    public DaVinci(Grid grid){
        this.grid = grid;
        
        // Initialize dimensions
        this.xDim = SimParams.GRID_SIZE;
        this.yDim = SimParams.GRID_SIZE;
        
        if (SimParams.PLOT_LIVE_IMAGES) {
            // LIVE MODE: Create windows for interactive display
            this.gridWin = new GridWindow("Tumour and oxygen", xDim, yDim, 2);
            this.plotWin = new PlotWindow("Radioligand versus time", 250, 250, 4, 0, 0, 1, 0.000001);
            this.testLine = new PlotLine(this.plotWin, GREEN);
            this.pixelBuffer = null;  // Not needed - use gridWin
        } else {
            // HEADLESS MODE: No windows, use offscreen buffer
            this.gridWin = null;
            this.plotWin = null;
            this.testLine = null;
            this.pixelBuffer = new int[xDim * yDim];  // Offscreen buffer
        }
    }
    
    /* ------------------------------------------------------------
       Warm colour map for oxygen
       ------------------------------------------------------------ */
	private int WarmOxygenGradient(double normalized) {
		// Low oxygen (0.0) = dark brown/amber
		// High oxygen (1.0) = light yellow/gold
		int r = (int)(120 + 135 * normalized);  // 120 -> 255
		int g = (int)(80 + 175 * normalized);   // 80 -> 255
		int b = (int)(40 + 100 * normalized);   // 40 -> 140 (stays warm)
		return RGB(r/255.0, g/255.0, b/255.0);
	}
    /* ------------------------------------------------------------
       PLOTTING
       ------------------------------------------------------------ */
    public void plot(double t, double value){
        if (!SimParams.PLOT_LIVE_IMAGES) return;
        if (testLine == null) return;
        testLine.AddSegment(t, value);
    }

    /* ------------------------------------------------------------
       GRID DRAW
       ------------------------------------------------------------ */
    public void gridDraw(boolean[] maskList){
        // Draw to either live window or offscreen buffer
        
        if (SimParams.PLOT_LIVE_IMAGES) {
            // LIVE MODE: Draw to window
            gridWin.Clear(BLACK);
            drawToGridWindow(maskList);
        } else {
            // HEADLESS MODE: Draw to buffer
            clearBuffer();
            drawToBuffer(maskList);
        }
    }
    
    /**
     * Draw visualization to GridWindow (live mode)
     */
    private void drawToGridWindow(boolean[] maskList) {
        for (int i = 0; i < gridWin.length; i++) {
            Cell c = grid.GetAgent(i);

            if (c == null) {
				if (maskList[6]) {
					double oxygenConc = grid.oxygenGrid.Get(i);
					double minScale = 0.0;
					double maxScale = 15000.0;
					double normalized = (oxygenConc - minScale) / (maxScale - minScale);
					normalized = Math.max(0.0, Math.min(1.0, normalized));
//					int warmalized = WarmOxygenGradient(normalized);
//					gridWin.SetPix(i, warmalized);
					gridWin.SetPix(i, HeatMapBGR(normalized, 0, 1));

				}
                continue;
            }

            // Vessel logic
            if (c.type == SimParams.VESSEL) {
                if (c.blockedVessel) {
                    // BLOCKED VESSEL -> show as dark red
                    if (maskList[6]) {
//                        double oxygenConc = grid.oxygenGrid.Get(i);
                        gridWin.SetPix(i, SimParams.COLORLIST[6]);
                    }
                } else { // UNBLOCKED VESSEL -> show as brighter red
                    if (maskList[SimParams.VESSEL]) {
                        gridWin.SetPix(i, SimParams.COLORLIST[5]);
                    }
                }
                continue;
            }

            // Tumor cells
            if (maskList[c.type]) {
                gridWin.SetPix(i, c.color);
            }
        }
    }
    
    /**
     * Draw visualization to offscreen buffer (headless mode)
     */
    private void drawToBuffer(boolean[] maskList) {
        for (int i = 0; i < pixelBuffer.length; i++) {
            Cell c = grid.GetAgent(i);

            if (c == null) {
				if (maskList[6]) {
					double oxygenConc = grid.oxygenGrid.Get(i);
					double minScale = 0.0;
					double maxScale = 15000.0;
					double normalized = (oxygenConc - minScale) / (maxScale - minScale);
					normalized = Math.max(0.0, Math.min(1.0, normalized));
//					pixelBuffer[i] = WarmOxygenGradient(normalized);
					pixelBuffer[i] = HeatMapBGR(normalized, 0, 1);
                }
                continue;
            }

            // Vessel logic
            if (c.type == SimParams.VESSEL) {
                if (c.blockedVessel) {
                    if (maskList[6]) {
                        double oxygenConc = grid.oxygenGrid.Get(i);
                        pixelBuffer[i] = SimParams.COLORLIST[6];
                    }
                } else {
                    if (maskList[SimParams.VESSEL]) {
                        pixelBuffer[i] = SimParams.COLORLIST[5];
                    }
                }
                continue;
            }

            // Tumor cells
            if (maskList[c.type]) {
                pixelBuffer[i] = c.color;
            }
        }
    }
    
    /**
     * Clear the offscreen buffer to black
     */
    private void clearBuffer() {
        for (int i = 0; i < pixelBuffer.length; i++) {
            pixelBuffer[i] = BLACK;
        }
    }
    
    /**
     * Get pixel color at position (for DataLogger export)
     */
    public int GetPix(int x, int y) {
        int index = y * xDim + x;
        
        if (SimParams.PLOT_LIVE_IMAGES) {
            // Read from live window
            return gridWin.GetPix(x, y);
        } else {
            // Read from offscreen buffer
            return pixelBuffer[index];
        }
    }

    /* ------------------------------------------------------------
       Print O2 statistics (always allowed)
       ------------------------------------------------------------ */
    public void displayOxygenStats(int currentDay) {
        double avgO2 = grid.oxygenGrid.GetAvg();
        double minO2 = grid.oxygenGrid.GetMin();
        double maxO2 = grid.oxygenGrid.GetMax();

//        System.out.printf("Day %d | O2: avg=%.2e min=%.2e max=%.2e\n",
//                currentDay, avgO2, minO2, maxO2);
    }
}
