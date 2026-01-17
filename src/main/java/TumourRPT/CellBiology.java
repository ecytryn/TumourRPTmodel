package TumorRPT;

public class CellBiology {

    public Cell cell;


    public CellBiology(Cell cell){
        this.cell = cell;
    }


    public void DivProbCalc(){
        // this function calculates the division probability of the cell
        // based on several factors like the cell type, oxygen level, etc
        this.cell.divisionProb = 1;

        if (this.cell.type == SimParams.HYPOXIC){
            this.cell.divisionProb = SimParams.DIVISON_PROB_MAX * this.cell.oxygen/SimParams.P_O2_VESSEL;
        }
        if (this.cell.type == SimParams.NORMAL){
            // TODO: Magic number here for good visualization
            this.cell.divisionProb = 10*SimParams.DIVISON_PROB_MAX * this.cell.oxygen/SimParams.P_O2_VESSEL;
        }
    }



}
