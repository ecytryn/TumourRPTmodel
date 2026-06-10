package TumorRPT;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;

import java.util.ArrayList;


class Cell extends AgentSQ2Dunstackable<Grid> {
    public int type = -1;
    public int color;
    public boolean isAlive;
    public boolean divisionFlag ;
    public boolean disposeFlag;
    public boolean blockedVessel;
    public double survivalProb;
    public double divisionProb;
    public double oxygen;
    public double vesselAmp;
    public int age; // in T units.
    public int birthTime;

    public int howManyHoursDead; // changed from howManyDaysDead with the time step fix


    public FSM_DIVCHECK fsmDiv;
    public CellBiology cellBio;

    void Init(int type, int currentTime){
        // current time is in T units
        this.type = -1;
        this.ChangeType(type);
        this.fsmDiv = new FSM_DIVCHECK(this);
        this.cellBio = new CellBiology(this);
        this.divisionFlag = false;
        this.disposeFlag = false;
        this.howManyHoursDead = 0;
        this.birthTime = currentTime;
        if (type== SimParams.VESSEL && currentTime == 0){

            // Make Less Dense
//           if (this.G.rng.Double() <0.3) {  // turn on initial vessel blocking
           if (this.G.rng.Double() <0.0) {  // turn off initial vessel blocking
               this.blockedVessel = true;
           }else{
               this.blockedVessel = false;
           }

           //Variety in Amp
           this.vesselAmp = Math.abs(this.G.rng.Gaussian(1,0.3));
        }
    }

    void ChangeType(int type){
        // Update population tracking
        if (this.type == -1){
            // Brand new cell - just add to population count
            G.CurrentCellsPops[type] += 1;
        }else{
            // Existing cell changing type - update counts
            G.CurrentCellsPops[this.type] -= 1;
            G.CurrentCellsPops[type] += 1;
        }
        this.type = type;
        this.color = SimParams.COLORLIST[type];

        // Set alive status based on type
        if (type == SimParams.NECROTIC || type == SimParams.APOPTOTIC){
            this.isAlive = false;
        } else {
            this.isAlive = true;  // NORMAL, HYPOXIC, VESSEL, HEALTHY are alive
        }
    }

    void Step(int currentHour, int currentDay){
        this.oxygen = this.G.oxygenGrid.Get(this.Isq());

        this.cellBio.DivProbCalc(); // updates divisionProb
        this.fsmDiv.FSM_Run(); // updates type and flags
        this.CellFate(currentHour, currentDay); // takes the actions

        if (this.isAlive == false){
            // For the forced removal of the dead cells that are not removed because of the
            // high standard deviation of the geometric random variable
            this.howManyHoursDead += 1;
        }
    }



    void CellFate(int currentHour, int currentDay){
        // Using the updated parameters and the flags, will take the appropriate
        // action, i.e. dispose cells, do the division.

        if (divisionFlag == true){
            // creates new cell.

            int emptySites = MapEmptyHood(SimParams.divHood);
            if (emptySites > 0){
                int selectedSiteIndex = G.rng.Int(emptySites);
                int newAgentSiteIndex = SimParams.divHood[selectedSiteIndex];
                G.NewAgentSQ(newAgentSiteIndex).Init(type, SimParams.globalTime);
            }
            this.divisionFlag = false;
        }
        if (disposeFlag == true){
            Dispose();
            G.CurrentCellsPops[this.type] -= 1;
        }
    }

    int GetSurvivalProb(){
        // Query RadioBio directly using birthTime and current type
        // RadioBio maintains ODE states for each birth hour
        this.survivalProb = this.G.radioBio.calculateSF(this.birthTime, this.type);
        return 0;
    }
}
