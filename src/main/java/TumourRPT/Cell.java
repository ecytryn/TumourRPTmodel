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

    public int howManyDaysDead;


    public FSM_DIVCHECK fsmDiv;
    public CellBiology cellBio;

    void Init(int type, int currentTime){
        // current time is in T units
        this.type = -1;
        this.ChangeType(type);
        this.fsmDiv = new FSM_DIVCHECK(this);
        this.cellBio = new CellBiology(this);
        this.isAlive = true;
        this.divisionFlag = false;
        this.disposeFlag = false;
        this.howManyDaysDead = 0;
//        this.age = 0;
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
        if (this.type == -1){
            G.CurrentCellsPops[type] += 1;
        }else{
            G.CurrentCellsPops[this.type] -= 1;
            G.CurrentCellsPops[type] += 1;
        }
        this.type = type;
        this.color = SimParams.COLORLIST[type];

        if (type == SimParams.NECROTIC || type==SimParams.APOPTOTIC){
            this.isAlive = false;
        }
    }

    void Step(int currentHour, int currentDay, ArrayList<double[]> lookupTable){

        this.oxygen = this.G.oxygenGrid.Get(this.Isq());

        this.cellBio.DivProbCalc(); // updates divisionProb
        this.fsmDiv.FSM_Run(lookupTable); // updates type and flags
        this.CellFate(currentHour, currentDay); // takes the actions

        if (this.isAlive == false){
            // For the forced removal of the dead cells that are not removed because of the
            // high standard deviation of the geometric random variable
            this.howManyDaysDead += 1;
        }
//        else{
//            this.age += SimParams.MinorLoop_PER_MajorLoop() * SimParams.T_PER_Minor() + currenHour;
//        }
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

    int GetSurvivalProb(ArrayList<double[]> lookupTable){
		int ageHours = SimParams.globalTime - birthTime;
		int lookUpTableAge = (int)(Math.floorDiv(ageHours, 24));  // Convert hours → days

        if (lookUpTableAge == 0){
            this.survivalProb = lookupTable.get(0)[0];
            // or equivalently we can write this.survivalProb = lookupTable.get(0)[0];
            // this is because at age zero that happens only on day 1 it doesn't matter what
            // is the type of cell.
            return 0;
        }

        if (lookUpTableAge > SimParams.maxLookupAge){
            lookUpTableAge = SimParams.maxLookupAge;
        }

        if (this.type == SimParams.NORMAL) {
            this.survivalProb = lookupTable.get(lookUpTableAge-1)[0];
        }
        if (this.type == SimParams.HYPOXIC){
            this.survivalProb = lookupTable.get(lookUpTableAge-1)[1];
        }

        return 0;
    }
}
