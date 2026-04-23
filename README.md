# Spatially-Resolved RPT Model

Hybrid cellular-automaton ODE/PDE model of radiopharmaceutical therapy (RPT) for solid tumours, coupling radioligand pharmacokinetics to spatially-resolved tumour dynamics. Implemented in Java using the [HAL](https://github.com/MathOnco/HAL) cellular-automaton modelling framework.

Associated manuscript:

> Mollaheydar, Saboury, Rahmim, and Cytrynbaum. *Mechanistic insights for
> radiopharmaceutical therapy: a spatially-resolved model coupling radioligand
> pharmacokinetics to tumour dynamics.* bioRxiv (2025).
> https://doi.org/10.1101/2025.05.28.656722

---

## Requirements

- **Java** 17+, with Gradle (wrapper included)
- **Python** 3.8+ with `numpy`, `matplotlib`, `scipy`, `pandas`, `seaborn`

---

## Quick Start

### 1. Run a single experiment

```bash
./gradlew run --args="NormoxicSmall"
```

Available experiments (defined in `Main.java`):

| Name             | Description                                                                   |
| ---------------- | ----------------------------------------------------------------------------- |
| `WatchGrow`      | No-treatment run; tumour growth, vessel occlusion, cell-type transitions      |
| `NormoxicSmall`  | Normoxic tumour, small (below cure threshold), injection day 5                |
| `NormoxicMedium` | Normoxic tumour, medium (above cure threshold), injection day 5               |
| `NormoxicLarge`  | Normoxic tumour, large (above cure threshold), injection day 5                |
| `HypoxicSmall`   | Hypoxic tumour (40-day pre-sim), small                                        |
| `HypoxicMedium`  | Hypoxic tumour (40-day pre-sim), medium                                       |
| `HypoxicLarge`   | Hypoxic tumour (40-day pre-sim), large — illustrates reoxygenation kill cycle |
| `CustomRun`      | One-off run; edit parameters directly in `Main.java`                          |

Output is written to `results/single_runs/<ExperimentName>_<timestamp>/`.

### 2. Run parameter sweeps

**Dose–receptor sweep** (injected amount vs. receptor density heatmap):
```bash
./gradlew runDoseReceptorSweep
```

For incremental overnight batching (appends to existing results):
```bash
./gradlew runDoseReceptorSweep --args="<timestamp>"
```

**Interval–skew sweep** (injection timing and injected amount distribution):
```bash
./gradlew runIntervalSkewSweep
```

**Tumour size threshold sweep** (cure rate vs. tumour size at injection):
```bash
./gradlew runTumourSizeThresholdSweep
```

Sweep output is written to `results/sweeps/<SweepName>_<timestamp>/`.

---

## Reproducing Figures

Figures 1–4 are snapshots exported directly by `Main.java` during a simulation run.

| Figure          | Script                                                                                   | Notes                                                                                              |
| --------------- | ---------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| 5 (final panel) | `scripts/pop_and_dose_full_viz.py`                                                       | Pass path to a single run directory                                                                |
| 6               | `scripts/tumour_cellcount_vs_cure_plot.py`, `scripts/tumour_size_vs_elimination_time.py` | Pass one or two `TumourSizeThresholdSweep` timestamps; use `--labels` for capillary density labels |
| 7               | `scripts/dose_receptor_sweep_fitting_half_cure.py`                                       | Auto-finds latest `DoseReceptorSweep` directory or pass a timestamp                                |
| 8               | `scripts/dose_receptor_dose_saturation_1panel.py`                                        | Auto-finds latest `DoseReceptorSweep` directory or pass a timestamp                                |
| 9               | `scripts/interval_skew_visualize.py`                                                     | Auto-finds latest `IntervalSkewSweep` directory or pass a timestamp                                |
| 10              | `scripts/interval_skew_tumour_size_determines_cure.py`                                   | Combines IS sweep and tumour size sweep data                                                       |
| 11              | `scripts/pk_comparison_plot.py`                                                          | Runs standalone; no simulation data needed                                                         |
| 12              | Manually drawn                                                                           | —                                                                                                  |
| 13              | `scripts/beta_retention_best_fit_ell_TJ.py`                                              | Runs standalone; no simulation data needed                                                         |

---

## Python Script Reference

### Figure scripts
| Script | Purpose |
|---|---|
| `beta_retention_best_fit_ell_TJ.py` | Fits mean beta range ℓ to Monte Carlo lookup table; plots uniform-sphere EDF vs. MC data (Fig 13) |
| `dose_receptor_dose_saturation_1panel.py` | Captive radioligand vs. time for selected receptor densities across injected amounts (Fig 8) |
| `dose_receptor_sweep_fitting_half_cure.py` | Dose–receptor cure-rate heatmap with fitted half-cure level curve (Fig 7) |
| `interval_skew_tumour_size_determines_cure.py` | Cure rate vs. tumour size at last injection, combining IS and tumour-size sweeps (Fig 10) |
| `interval_skew_visualize.py` | Interval–skew cure-rate heatmap (Fig 9) |
| `pk_comparison_plot.py` | PK model validation plot (Fig 11) |
| `pop_and_dose_full_viz.py` | Cell populations and dose rate over time for a single run (Fig 5) - pass path to data|
| `tumour_cellcount_vs_cure_plot.py` | Cure rate vs. tumour cell count at injection, one curve per capillary density (Fig 6) |
| `tumour_size_vs_elimination_time.py` | Mean time to tumour elimination vs. tumour size at injection (Fig 6) |

### Utility scripts
| Script | Purpose |
|---|---|
| `debug_helper.py` | Interactive helper for identifying and re-running interesting points from a sweep heatmap |
| `sf_estimate.py` | Estimates empirical survival fraction from population time series; useful for validating radiobiology |

### Vessel generation (in `scripts/GenerateVessels/`)
| Script | Purpose |
|---|---|
| `GenerateUniformVesselDensity.py` | Generates vessel configurations with a specified target density using repulsion-based placement |

---

## Java Source Overview

| File | Purpose |
|---|---|
| `Main.java` | Entry point; experiment configurations |
| `SimParams.java` | All model parameters (SI units throughout) |
| `Grid.java` | Spatial domain; agent management and stepping |
| `Cell.java` | Individual cell agent |
| `CellBiology.java` | Cell division logic |
| `FSM_DIVCHECK.java` | Cell state machine (normoxic ↔ hypoxic ↔ necrotic ↔ apoptotic) |
| `PK.java` | Radioligand pharmacokinetics (multi-compartment ODE) |
| `RadioBio.java` | Radiobiology: LQ model with Lea-Catcheside G-factor (ODE-based) |
| `Oxygen.java` | Oxygen dynamics wrapper |
| `OxygenDiffusionSolver.java` | Steady-state O₂ PDE solver (red-black SOR) |
| `BetaRetention.java` | Beta particle dose deposition fraction (lookup table) |
| `DaVinci.java` | Visualisation (tumour and oxygen field images) |
| `DataLogger.java` | CSV output and image export |
| `VesselConfigConvertor.java` | Loads vessel configurations from CSV |
| `MyUtils.java` | General utilities |
| `DoseReceptorSweep.java` | Parameter sweep runner: injected amount × receptor density |
| `IntervalSkewSweep.java` | Parameter sweep runner: injection interval × injected amount skew |
| `TumourSizeThresholdSweep.java` | Parameter sweep runner: cure rate vs. initial tumour radius |
| `DebugRunner.java` | Runs a single parameter set with detailed logging |

---

## Output Structure

```
results/
├── single_runs/
│   └── <ExperimentName>_<timestamp>/
│       ├── parameters.md         # human-readable parameter report
│       ├── parameters.csv        # machine-readable parameters
│       ├── populations.csv       # cell-type counts (hourly)
│       ├── dose.csv              # dose rate over time
│       ├── pkStateVariables.csv  # PK compartment states
│       └── tumour_images/        # spatial snapshots (every 10 days)
└── sweeps/
    └── <SweepName>_<timestamp>/
        ├── sweep_summary.csv     # outcome per parameter combination
        └── <param_combination>/  # per-run subdirectories (as above)
```

---

## Notes

- All parameters in `SimParams.java` are in SI units (metres, seconds, moles).
- The grid is 400 × 400 cells at 10 µm spacing (4 mm × 4 mm domain).
- Reliable tumour radius range: up to ~950 µm initial radius.
- The `OldScripts/` directory contains superseded and exploratory scripts retained for reference.