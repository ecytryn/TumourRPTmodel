# PK.java Migration Notes - Your QSS Approach

## Overview

Your existing PK model uses a clever **quasi-steady-state (QSS) approximation** that I initially missed. This is much more efficient than integrating all compartments!

## Key Features of Your Approach

### 1. Timescale Separation
```
FAST (QSS):     C_v, C_int, C_b     (equilibrate quickly with C_BP)
SLOW (ODEs):    C_BP, C_intern      (integrated over time)
```

### 2. Hot/Cold Tracking
- **Hot**: Radioactive Lu-177 (decays with half-life ~6.7 days)
- **Cold**: Decay product (stable, cleared biologically)
- Physical decay converts hot → cold

### 3. Your ODE System
```
dC_BP_hot/dt     = -(λ_bio + λ_decay)·C_BP_hot
dC_BP_cold/dt    = -λ_bio·C_BP_cold + λ_decay·C_BP_hot
dC_intern_hot/dt = k_int·C_b_hot - k_rel·C_intern_hot - λ_decay·C_intern_hot
dC_intern_cold/dt= k_int·C_b_cold - k_rel·C_intern_cold + λ_decay·C_intern_hot
```

Where the QSS compartments are computed algebraically:
```
C_v ≈ C_BP                                (fast equilibration)
C_int ≈ C_BP                              (fast equilibration)
C_b = (R_T/V_int)·C_BP/(C_BP_total + β)  (binding equilibrium)
    with β = (k_off + k_int)/k_on
```

## Migration Changes

### What Changed in the Migrated PK.java

#### 1. **Parameter References**
```java
// OLD (from ParamsPBPK.java)
ParamsPBPK.lambda_bio    // in 1/hour
ParamsPBPK.k_on          // in L/(nmol·hour)
ParamsPBPK.V_BP          // in L

// NEW (from SimParams.java)
SimParams.LAMBDA_BIO     // in 1/s
SimParams.K_ON           // in m³/(mol·s)
SimParams.V_CENTRAL      // in m³
```

#### 2. **No Unit Conversions in Math**
```java
// OLD: Had conversion factors scattered
derivatives[0] = -ParamsPBPK.lambda_bio * C_BP_hot  // C_BP in nmol/L, lambda in 1/h

// NEW: Direct SI units
derivatives[0] = -SimParams.LAMBDA_BIO * C_BP_hot   // C_BP in mol/m³, lambda in 1/s
```

#### 3. **Time Arguments**
```java
// OLD: DoseRateCalc_DEBUG(t_0, t_f, day, hour)
//      t_0, t_f were in HOURS

// NEW: DoseRateCalc(t_0, t_f, day, hour)
//      t_0, t_f are in SECONDS
```

#### 4. **QSS Functions**
The QSS calculation functions now pull geometry directly from Grid:
```java
private double calc_C_b_QSS(double C_BP, double C_BP_total) {
    // Get current tumor state
    int tumorCells = grid.countTumorCells();
    int numVessels = grid.countVessels();
    double receptorMoles = SimParams.computeReceptorMoles(tumorCells, numVessels);
    double tumorVolume = SimParams.computeTumorVolume(tumorCells);
    double V_int = SimParams.computeInterstitialVolume(tumorVolume);
    
    double R_T_tilde = receptorMoles / V_int;  // mol/m³
    double beta = (SimParams.K_OFF + SimParams.K_INT) / SimParams.K_ON;
    
    return R_T_tilde * C_BP / (C_BP_total + beta);
}
```

This replaces the old approach where R_T was updated once per hour in ParamsPBPK.perHourUpdateParams().

## Required Changes in Other Files

### 1. Grid.java Needs These Methods

```java
/**
 * Count tumor cells (excluding vessels)
 */
public int countTumorCells() {
    int count = 0;
    for (Cell cell : this) {
        if (cell != null && cell.type != SimParams.VESSEL) {
            count++;
        }
    }
    return count;
}

/**
 * Count vessels in domain
 */
public int countVessels() {
    int count = 0;
    for (Cell cell : this) {
        if (cell != null && cell.type == SimParams.VESSEL) {
            count++;
        }
    }
    return count;
}
```

### 2. Injection Handling

Your existing injection mechanism in ParamsPBPK.Inject() needs updating:

```java
// OLD units (from ParamsPBPK.java)
currentStateVariables[C_BP_hot_index] += 
    hotFraction[i] * injectionAmounts[i] / V_BP;  // nmol / L = nmol/L

// NEW units (need to implement in new system)
currentStateVariables[C_BP_hot_index] += 
    hotFraction[i] * injectionAmounts[i] / SimParams.V_CENTRAL;  // mol / m³ = mol/m³
```

Where `injectionAmounts[i]` is now in MOLES (not nmol).

### 3. Main.java Time Loop

```java
// OLD: Time in hours
for (int hour = 0; hour < 24; hour++) {
    double t_0 = hour * 3600.0;        // Convert hour to seconds
    double t_f = (hour + 1) * 3600.0;  // Convert hour to seconds
    pk.DoseRateCalc(t_0, t_f, day, hour);
}

// NEW: Same, but be aware t_0 and t_f are in seconds!
for (int hour = 0; hour < 24; hour++) {
    double t_0 = hour * SimParams.TIME_STEP;
    double t_f = (hour + 1) * SimParams.TIME_STEP;
    pk.DoseRateCalc(t_0, t_f, day, hour);
}
```

## State Vector Structure

Your state vector remains the same structure but different units:

```
Index  Variable         OLD Units    NEW Units
-----  ---------------  -----------  ----------
  0    C_BP_hot         nmol/L       mol/m³
  1    C_BP_cold        nmol/L       mol/m³
  2    C_v_hot          nmol/L       mol/m³
  3    C_v_cold         nmol/L       mol/m³
  4    C_int_hot        nmol/L       mol/m³
  5    C_int_cold       nmol/L       mol/m³
  6    C_b_hot          nmol/L       mol/m³
  7    C_b_cold         nmol/L       mol/m³
  8    C_intern_hot     nmol/L       mol/m³
  9    C_intern_cold    nmol/L       mol/m³
 10    A_blob           nmol         mol
```

Note: Indices 2-9 are QSS (not integrated), computed from C_BP.

## Dose Rate Calculation

The dose rate calculation changes slightly:

```java
// OLD approach (from your PK.java)
double A_hot = (C_v_hot * V_v + C_int_hot * V_int + 
                C_b_hot * V_int + C_intern_hot * V_int) * lambda_phys;
                // (nmol/L * L) * (1/h) = nmol/h

double doseRate = A_hot / Activity_PER_Bq() * E_beta / (totalPop * M_cell);
                // Complex unit conversions hidden in Activity_PER_Bq()

// NEW approach (cleaner SI)
double A_hot = (C_v_hot * V_v + C_int_hot * V_int + 
                C_b_hot * V_int + C_intern_hot * V_int) * LAMBDA_DECAY;
                // (mol/m³ * m³) * (1/s) = mol/s

double activityBq = A_hot * AVOGADRO;  // mol/s * (1/mol) = Bq
double energyRate = activityBq * E_BETA;  // Bq * J = J/s
double doseRate_Gy_per_s = energyRate / tumorMass;  // J/s / kg = Gy/s
```

## Validation Strategy

1. **Run both models side-by-side**
   - Keep your old code in one branch
   - Run new code with same injection schedule
   - Compare PK curves (after unit conversion)

2. **Check these quantities match:**
   - C_BP_hot decay curve (exponential with correct half-life)
   - Total radioligand mass balance
   - Dose rate profiles
   - Final tumor response

3. **Unit conversion for comparison:**
   ```
   OLD: nmol/L
   NEW: mol/m³
   Conversion: 1 nmol/L = 1e-9 mol / 1e-3 m³ = 1e-6 mol/m³
   
   So: NEW_value = OLD_value * 1e-6
   ```

## Critical Points

✅ **Preserve Your QSS Approach** - Don't integrate fast compartments
✅ **Keep Hot/Cold Tracking** - Essential for accurate dose calculation
✅ **Use Grid Methods** - Don't store R_T as static variable anymore
✅ **Time is Seconds** - ALL time arguments and rate constants in 1/s
✅ **Validate Carefully** - PK curves should match exactly after unit conversion

## Next Steps

1. Add `countTumorCells()` and `countVessels()` methods to Grid.java
2. Update injection mechanism to use SI units
3. Update Main.java time loop (minimal changes)
4. Test with frozen tumor first
5. Compare PK curves to old model
6. Enable tumor growth once PK validates
