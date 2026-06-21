package TumourRPT;

public class CellBiology {

    public Cell cell;


    public CellBiology(Cell cell){
        this.cell = cell;
    }


    public void DivProbCalc(){
        // this function calculates the division probability of the cell
        // based on several factors like the cell type, oxygen level, etc

	    // Default for vessels and necrotic, apoptotic cells
		this.cell.divisionProb = 0.0;

        if (this.cell.type == SimParams.HYPOXIC){
            this.cell.divisionProb = 3*SimParams.DIVISION_PROB_MAX * this.cell.oxygen/SimParams.P_O2_VESSEL;
        }
        if (this.cell.type == SimParams.NORMAL){
            // The 3* ensures division when oxygen is not suppressed.
            this.cell.divisionProb = 3*SimParams.DIVISION_PROB_MAX * this.cell.oxygen/SimParams.P_O2_VESSEL;
        }
    }



}
