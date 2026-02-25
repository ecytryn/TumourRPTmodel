#!/bin/bash

# Automated Debug Sweep for DoseReceptorSweep
# USAGE: ./dose_receptor_debug_sweep.sh [sweep_directory]
# Example: ./dose_receptor_debug_sweep.sh results/DoseReceptorSweep/DoseReceptorSweep_2026-02-24_13-00-00

# =============================================================================
# AUTO-DETECT OR USE PROVIDED SWEEP DIRECTORY
# =============================================================================

SWEEP_DIR="$1"

if [ -z "$SWEEP_DIR" ]; then
    # Auto-find most recent DoseReceptorSweep
    SWEEP_DIR=$(ls -td results/DoseReceptorSweep/DoseReceptorSweep_* 2>/dev/null | head -1)
fi

if [ ! -d "$SWEEP_DIR" ]; then
    echo "ERROR: Sweep directory not found: $SWEEP_DIR"
    echo ""
    echo "USAGE: $0 [sweep_directory]"
    echo "Example: $0 results/DoseReceptorSweep/DoseReceptorSweep_2026-02-24_13-00-00"
    echo ""
    echo "Or run without arguments to auto-detect most recent sweep."
    exit 1
fi

echo "============================================================"
echo "  DOSE-RECEPTOR DEBUG SWEEP"
echo "============================================================"
echo ""
echo "Using sweep directory:"
echo "  $SWEEP_DIR"
echo ""

# =============================================================================
# EXTRACT CONFIGURATION FROM PARAMETERS.CSV
# =============================================================================

PARAM_FILE=$(find "$SWEEP_DIR" -name "parameters.csv" -type f | head -1)

if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: No parameters.csv found in $SWEEP_DIR"
    exit 1
fi

echo "Reading configuration from:"
echo "  $PARAM_FILE"
echo ""

# Extract fixed parameters
INITIAL_RADIUS=$(grep "^initial_tumor_radius," "$PARAM_FILE" | cut -d',' -f2)
VESSEL_CONFIG=$(grep "^vessel_density_config," "$PARAM_FILE" | cut -d',' -f2)
HYPOXIA_DEV=$(grep "^hypoxia_dev_days," "$PARAM_FILE" | cut -d',' -f2 2>/dev/null || echo "0")

echo "Configuration from sweep:"
echo "  Initial radius: ${INITIAL_RADIUS} um"
echo "  Vessel config: ${VESSEL_CONFIG}"
echo "  Hypoxia development: ${HYPOXIA_DEV} days"
echo ""

# =============================================================================
# EXTRACT SWEEP PARAMETER RANGES
# =============================================================================

# Get all unique doses and receptor densities from directory names
# Directory format: dose_150_recep_5.00e-19_rep_1

DOSES=($(find "$SWEEP_DIR" -type d -name "dose_*_recep_*_rep_*" | \
         sed 's/.*dose_\([0-9]*\)_recep.*/\1/' | sort -nu))

# For receptors, we need to be more careful with scientific notation
RECEPTORS=($(find "$SWEEP_DIR" -type d -name "dose_*_recep_*_rep_*" | \
            sed 's/.*recep_\([0-9.e-]*\)_rep.*/\1/' | sort -u))

echo "Sweep parameter ranges found:"
echo "  Doses: ${DOSES[@]} nmol"
echo "  Receptor densities: ${#RECEPTORS[@]} values"
echo ""

# =============================================================================
# USER SELECTS ROW OR COLUMN TO DEBUG
# =============================================================================

echo "Select debug mode:"
echo "  1) Dose row (fix receptor density, vary dose)"
echo "  2) Receptor column (fix dose, vary receptor density)"
echo ""
read -p "Choice [1-2]: " MODE_CHOICE

if [ "$MODE_CHOICE" = "1" ]; then
    # Dose row - select which receptor density to fix
    echo ""
    echo "Enter receptor multiplier (e.g., 1.0 for baseline, 1.2 for 120%):"
    read -p "Receptor multiplier: " FIXED_RECEPTOR
    
    SWEEP_VALUES=("${DOSES[@]}")
    SWEEP_PARAM="dose"
    FIXED_PARAM="receptor"
    FIXED_VALUE=$FIXED_RECEPTOR
    
    SWEEP_NAME="dose_row_receptor${FIXED_RECEPTOR}"
    OUTPUT_BASE="results/debug_sweeps/${SWEEP_NAME}_$(date +%Y%m%d_%H%M%S)"
    
elif [ "$MODE_CHOICE" = "2" ]; then
    # Receptor column - select which dose to fix
    echo ""
    echo "Available doses: ${DOSES[@]} nmol"
    read -p "Select dose value (nmol): " FIXED_DOSE
    
    # Create receptor multiplier array (assuming baseline from paramfile)
    BASELINE=$(grep "^receptors_baseline_mol," "$PARAM_FILE" | cut -d',' -f2)
    RECEPTOR_MULTS=()
    for R in "${RECEPTORS[@]}"; do
        MULT=$(echo "scale=2; $R / $BASELINE" | bc -l)
        RECEPTOR_MULTS+=($MULT)
    done
    
    SWEEP_VALUES=("${RECEPTOR_MULTS[@]}")
    SWEEP_PARAM="receptor"
    FIXED_PARAM="dose"
    FIXED_VALUE=$FIXED_DOSE
    
    SWEEP_NAME="dose${FIXED_DOSE}_receptor_col"
    OUTPUT_BASE="results/debug_sweeps/${SWEEP_NAME}_$(date +%Y%m%d_%H%M%S)"
else
    echo "Invalid choice. Exiting."
    exit 1
fi

echo ""
read -p "Number of replicates per point [5]: " NUM_REPLICATES
NUM_REPLICATES=${NUM_REPLICATES:-5}

mkdir -p "$OUTPUT_BASE"

echo ""
echo "Debug sweep configuration:"
echo "  Mode: $MODE_CHOICE (${SWEEP_PARAM} sweep, ${FIXED_PARAM}=${FIXED_VALUE})"
echo "  Values to test: ${SWEEP_VALUES[@]}"
echo "  Replicates: ${NUM_REPLICATES}"
echo "  Output: ${OUTPUT_BASE}"
echo ""
echo "This will run $(( ${#SWEEP_VALUES[@]} * ${NUM_REPLICATES} )) simulations"
echo ""
read -p "Press Enter to continue or Ctrl+C to cancel..."

# =============================================================================
# RUN SIMULATIONS
# =============================================================================

LOGFILE="${OUTPUT_BASE}/sweep_log.txt"
echo "Debug sweep started: $(date)" > "${LOGFILE}"
echo "Source sweep: $SWEEP_DIR" >> "${LOGFILE}"
echo "" >> "${LOGFILE}"

echo ""
echo "Starting simulations..."
echo ""

TOTAL_RUNS=$(( ${#SWEEP_VALUES[@]} * ${NUM_REPLICATES} ))
CURRENT_RUN=0

for REP in $(seq 1 $NUM_REPLICATES); do
    echo "=========================================="
    echo "  REPLICATE ${REP}/${NUM_REPLICATES}"
    echo "=========================================="
    echo ""
    
    for VALUE in "${SWEEP_VALUES[@]}"; do
        CURRENT_RUN=$((CURRENT_RUN + 1))
        
        # Set dose and receptor based on mode
        if [ "$MODE_CHOICE" = "1" ]; then
            DOSE=$VALUE
            RECEPTOR=$FIXED_VALUE
        else
            DOSE=$FIXED_VALUE
            RECEPTOR=$VALUE
        fi
        
        echo "[${CURRENT_RUN}/${TOTAL_RUNS}] Running: dose=${DOSE} nmol, receptor=${RECEPTOR}x, replicate=${REP}"
        
        START_TIME=$(date +%s)
        SEED=$((42 + REP))
        
        # Run with parameter file
        ./gradlew runDebug --args="dose=${DOSE} receptors=${RECEPTOR} seed=${SEED} paramfile=${PARAM_FILE}" \
            2>&1 | tee -a "${LOGFILE}"
        
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        
        # Find most recent debug output
        LATEST_DIR=$(ls -td results/debug_runs/Debug_D*_R*_* 2>/dev/null | head -1)
        
        if [ -z "${LATEST_DIR}" ]; then
            echo "ERROR: Could not find output directory"
            continue
        fi
        
        # Move to organized location
        TARGET="${OUTPUT_BASE}/dose_${DOSE}_receptor_${RECEPTOR}_rep_${REP}"
        mv "$LATEST_DIR" "$TARGET"
        
        echo "  -> Saved to: $TARGET"
        echo "  -> Runtime: ${ELAPSED}s"
        echo ""
    done
done

echo ""
echo "============================================================"
echo "  DEBUG SWEEP COMPLETE"
echo "============================================================"
echo ""
echo "Results saved in: ${OUTPUT_BASE}"
echo ""