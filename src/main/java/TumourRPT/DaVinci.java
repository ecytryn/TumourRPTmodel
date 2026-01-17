package TumorRPT;

import HAL.Gui.GridWindow;
import HAL.Gui.PlotLine;
import HAL.Gui.PlotWindow;

import java.util.ArrayList;

import static HAL.Util.*;

class DaVinci {
    private Grid grid;
    public GridWindow gridWin;
    public GridWindow gridWinAge;

    public PlotWindow plotWin;

    public PlotLine testLine;

    public DaVinci(Grid grid){
        this.grid = grid;

        // ----------------------------
        // Visualization Toggle
        // ----------------------------
        if (SimParams.VISUALIZATION_ON) {
            this.plotWin    = new PlotWindow("Radioligand versus time",250,250,4,0,0,1,0.000001);
            this.gridWin    = new GridWindow("Tumour and oxygen", SimParams.GRID_SIZE, SimParams.GRID_SIZE, 2);
//            this.gridWinAge = new GridWindow("Age", SimParams.GRID_SIZE, SimParams.GRID_SIZE, 2);

            this.testLine   = new PlotLine(this.plotWin, GREEN);
        } 
        else {
            // HEADLESS MODE
            this.plotWin    = null;
            this.gridWin    = null;
//            this.gridWinAge = null;
            this.testLine   = null;
        }
    }

    /* ------------------------------------------------------------
       PLOTTING
       ------------------------------------------------------------ */
    public void plot(double t, double value){
        if (!SimParams.VISUALIZATION_ON) return;
        if (testLine == null) return;
        testLine.AddSegment(t, value);
    }

    /* ------------------------------------------------------------
       GRID DRAW
       ------------------------------------------------------------ */
    public void gridDraw(boolean[] maskList){
        if (!SimParams.VISUALIZATION_ON) return;
        if (gridWin == null) return;

        gridWin.Clear(BLACK);

        for (int i = 0; i < gridWin.length; i++) {
            Cell c = grid.GetAgent(i);

            if (c == null) {
                if (maskList[6]) {
                    double oxygenConc = grid.oxygenGrid.Get(i);
                    gridWin.SetPix(i, HeatMapBGR(oxygenConc, 0, 
                        1.5 * SimParams.P_O2_VESSEL));
                }
                continue;
            }

            // Vessel logic
            if (c.type == SimParams.VESSEL) {
                if (c.blockedVessel) {
                    // BLOCKED VESSEL → show oxygen instead (NOT vessel marker)
                    if (maskList[6]) {
                        double oxygenConc = grid.oxygenGrid.Get(i);
                        gridWin.SetPix(i, RGB(0.5, 0.0, 0.0));
                    }
                } else {
                    if (maskList[SimParams.VESSEL]) {
                        gridWin.SetPix(i, RGB(1.0, 0.3, 0.3)); // unblocked vessel
                    }
                }
                continue;
            }

            // Tumour or apoptotic/hypoxic/necrotic
            if (maskList[c.type]) {
                gridWin.SetPix(i, c.color);
            }
        }
    }

    /* ------------------------------------------------------------
       AGE DRAW
       ------------------------------------------------------------ */
    public void gridDrawAge(boolean[] maskList){
        if (!SimParams.VISUALIZATION_ON) return;
        if (gridWinAge == null) return;

        gridWinAge.Clear(BLACK);

        for (int i = 0; i < gridWinAge.length; i++) {
            Cell c = grid.GetAgent(i);

            if (c == null) {
                if (maskList[6]) {
                    double oxygenConc = grid.oxygenGrid.Get(i);
                    gridWinAge.SetPix(i, HeatMapBGR(oxygenConc, 0,
                        1.5 * SimParams.P_O2_VESSEL));
                }
                continue;
            }

            if (c.type == SimParams.VESSEL) {
                if (c.blockedVessel) {
                    // draw DEEP BLUE for blocked vessel age view
                    if (maskList[6]) {
                        gridWinAge.SetPix(i, RGB(0, 0, 0.2));
                    }
                } else {
                    if (maskList[SimParams.VESSEL]) {
                        gridWinAge.SetPix(i, RGB(0.4, 0.4, 0.4));
                    }
                }
                continue;
            }

            if (maskList[c.type]) {
                double ageDays = (SimParams.globalTime - c.birthTime)/24.0;
                gridWinAge.SetPix(i, HeatMapRGB(ageDays, 0, 20));
            }
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
