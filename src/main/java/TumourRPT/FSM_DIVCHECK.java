package TumorRPT;

import java.util.ArrayList;

public class FSM_DIVCHECK {

    // FSM_DIVCHECK will get the following parameters throught its methods
    // @params: N_pl, oxygen, type, survivalProb
    // and will return
    // return: divideFlag, IsAlive, disposeFlag, type

    public Cell cell;

    public FSM_DIVCHECK(Cell cell){
        // constructor of the class
        this.cell = cell;
    }


    public int FSM_Run(){
        // trouble reading the name? This is finite state machine run.
        // Logically, we first need to determine the new cell type (NormalTCell, HypoTCELL, NEC) based on the
        // oxygen parameters. Then, based on the updated type, (which is in fact the type of the cell that the cell was
        // in for the entire duration of the last day)
        if (this.cell.isAlive == true && this.cell.type != SimParams.VESSEL){
            double NecroThreshold = SimParams.P_O2_NECROTIC * this.cell.G.rng.Gaussian(1,0.2);
            if (this.cell.oxygen < NecroThreshold ){
                this.cell.ChangeType(SimParams.NECROTIC);
            }else {
//                double hypoThreshold = SimParams.P_O2_HYPOXIC * Math.abs(this.cell.G.rng.Gaussian(1,0.1));
                double hypoThreshold = SimParams.P_O2_HYPOXIC * this.cell.G.rng.Gaussian(1,0.2);
//                double hypoThreshold = SimParams.P_O2_HYPOXIC  + this.cell.G.rng.Gaussian(0,0.1* SimParams.P_O2_HYPOXIC);
                if (this.cell.oxygen < hypoThreshold) {
                    this.cell.ChangeType(SimParams.HYPOXIC);
                } else {
                    this.cell.ChangeType(SimParams.NORMAL);
                }
				// Check if cell wants to divide
				this.DivideCheck();

				// Only check radiation survival if attempting division
				if (this.cell.divisionFlag) {
					this.cell.GetSurvivalProb();
					if (this.cell.G.rng.Double() < (1 - this.cell.survivalProb)) {
						// Failed survival check during mitosis
						this.cell.ChangeType(SimParams.APOPTOTIC);
					} else {
						// Survived mitosis - reset radiation age for both parent and daughter
						this.cell.birthTime = SimParams.globalTime;
						// divisionFlag remains true, Grid will create daughter cell
					}
				}
				return 0;
            }
        }
        this.DisposeCheck();
        return 0;
    }

    public void DisposeCheck(){
        // if some conditions are satisfied, then sets the dispose flag to be true
        if (this.cell.type == SimParams.APOPTOTIC){
            if(this.cell.G.rng.Double()<SimParams.APOP_REMOVAL_PROB_PER_HOUR){
                this.cell.disposeFlag = true;
            }else{
                if (this.cell.howManyDaysDead >= SimParams.strictRemovalCoeff * SimParams.ApopRemovalTime){
                    this.cell.disposeFlag = true;
                }else{
                    this.cell.disposeFlag = false;
                }
            }
        }
//        if (this.cell.type == SimParams.NECROTIC){
//            if(this.cell.G.rng.Double()<SimParams.NECROTIC_REMOVAL_PROB){
//                this.cell.disposeFlag = true;
//            }
//            else{
//                if (this.cell.howManyDaysDead >= SimParams.strictRemovalCoeff * SimParams.NecroRemovalTime){
//                    this.cell.disposeFlag = true;
//                }else{
//                    this.cell.disposeFlag = false;
//                }
//            }
//        }

    }

    public void DivideCheck(){
        // if some conditions are satisfied, then sets the divide flag to true
        // this.cell.divideFlag = true.
        if (this.cell.type != SimParams.VESSEL) {
//  Added to check range of divProb
//			if (this.cell.G.GetTick() % 240 == 0) {  // Every 10 days
//				System.out.printf("Cell divProb sample: %.3f (type=%d, O2=%.1f)\n", 
//					 this.cell.divisionProb, this.cell.type, this.cell.oxygen);
//			}
            if (this.cell.G.rng.Double()<this.cell.divisionProb){
                this.cell.divisionFlag = true;
            }
            else {
                this.cell.divisionFlag = false;
            }
        }


    }

//    public void DivideRun(){
//        // this is a temporary function to test the output of the code
//        if (this.cell.type == SimParams.NORMAL || this.cell.type == SimParams.HYPOXIC){
//            int emptySites = this.cell.MapEmptyHood(SimParams.divHood);
//            if (emptySites > 0){
//                int selectedSiteIndex = this.cell.G.rng.Int(emptySites);
//                int newAgentSiteIndex = SimParams.divHood[selectedSiteIndex];
//                this.cell.G.NewAgentSQ(newAgentSiteIndex).Init(this.cell.type, );
//            }
//            if (this.cell.G.rng.Double()<0.25){
//                this.cell.Init(SimParams.NECROTIC);
//            }
//        }
//    }
}
