# Tumor Radioligand Therapy (RPT) Model

A hybrid CA/ODE/PDE model simulating tumor response to radioligand therapy, incorporating pharmacokinetics (ODE), oxygen dynamics (PDE), tumour population dynamics including radiation-induced cell death (CA handled by the Hybrid Automaton Library (HAL))

## Overview

This model simulates:
- **Pharmacokinetics**: Radioligand distribution, binding, and internalization
- **Oxygen dynamics**: Diffusion-consumption PDE with vessel delivery and vessel occlusion
- **Radiation effects**: LQ model with Lea-Catcheside G-factor for protracted delivery
- **Tumor dynamics**: Cell division, hypoxia-driven transitions, radiation-induced death

## Current Status (January 2026)

**Working:**
- ✓ PK model validated against analytical solutions
- ✓ Oxygen solver with calibrated boundary conditions
- ✓ Radiobiology with ODE-based G-factor calculation
- ✓ Beta particle escape correction to allow simplified geometry
- ✓ Vessel occlusion with crude mechanical-pressure model

**In Progress:**
- Converting PK from concentration-based to amount-based formulation
- Parameter tuning for realistic treatment response

## Quick Start

### Prerequisites
- Java 17+
- Gradle 7+

### Build
```bash
./gradlew build
```

### Run Single Simulation
```bash
./gradlew run
```

### Run Parameter Sweep - not tested with latest version, status unknown
```bash
./gradlew runSweep
```

## Model Structure

### Core Components
- `Main.java` - Simulation entry point
- `Grid.java` - Spatial domain and agent management
- `Cell.java` - Individual cell agent
- `PK.java` - Pharmacokinetic ODEs
- `RadioBio.java` - Survival probability calculations
- `Oxygen.java` / `OxygenDiffusionSolver.java` - Oxygen PDE solver
- `SimParams.java` - All parameters (SI units)

### Key Features
- **Quasi-steady-state approximation** - Fast compartments eliminated
- **ODE-based G-factor** - O(n) computation vs O(n²)
- **Modular design** - Easy to extend

## Parameters

All parameters in `SimParams.java` use SI units:
- Time: seconds (s)
- Length: meters (m)
- Amount: moles (mol)
- Concentration: mol/m³
- Pressure: Pascals (Pa)

See `SimParams.java` for complete parameter list with sources.

## Documentation

- `CURRENT_STATUS_*.md` - Development status summaries
- `ODE_FORMULATION_LEA_CATCHESIDE.md` - Mathematical derivations
- `FOLDER_STRUCTURE.md` - Project organization
- `TumourRPT_SupplementalMaterial.pdf` - Published model description

## Output

Simulations produce:
- `populations.csv` - Cell counts over time
- `dose.csv` - Dose rate history
- `pkStateVariables.csv` - PK compartment concentrations
- `day_XXX.png` - Spatial visualizations

## Testing

### Frozen Tumor Test
Set `FREEZE_TUMOR = true` in `SimParams.java` to test PK without cell dynamics.

### Validation
- PK curves match analytical solutions
- Oxygen distribution physiologically reasonable
- Dose rates in expected range (0.5-5 Gy per injection)

## Citation

If you use this model, please cite:
[Paper citation to be added]

## License

[License to be determined]

## Contact

[Contact information to be added]

## Acknowledgments

Model developed with assistance from Claude (Anthropic) in January 2026.
