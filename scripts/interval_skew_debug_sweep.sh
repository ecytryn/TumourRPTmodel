#!/bin/bash

# Automated Debug Sweep for IntervalSkewSweep
# USAGE: ./interval_skew_debug_sweep.sh [sweep_directory]
# Example: ./interval_skew_debug_sweep.sh results/IntervalSkewSweep/IntervalSkewSweep_2026-02-24_13-00-00

# =============================================================================
# AUTO-DETECT OR USE PROVIDED SWEEP DIRECTORY
# =============================================================================

SWEEP_DIR="$1"

if [ -z "$SWEEP_DIR" ]; then
    # Auto-find most recent IntervalSkewSweep
    SWEEP_DIR=$(ls -td results/IntervalSkewSweep/IntervalSkewSweep_* 2>/dev/null | head -1)
fi

if [ ! -d "$SWEEP_DIR" ]; then
    echo "ERROR: Sweep directory not found: $SWEEP_DIR"
    echo ""
    echo "USAGE: $0 [sweep_directory]"
    echo "Example: $0 results/IntervalSkewSweep/IntervalSkewSweep_2026-02-24_13-00-00"
    echo ""
    echo "Or run without arguments to auto-detect most recent sweep."
    exit 1
fi

echo "============================================================"
echo "  INTERVAL-SKEW DEBUG SWEEP"
echo "============================================================"
echo ""
echo "Using sweep directory:"
echo "  $SWEEP_DIR"
echo ""

# =============================================================================
# EXTRACT CONFIGURATION FROM PARAMETERS.CSV
# =============================================================================

# Find a sample parameters.csv (any one will do - all runs have same fixed params)
PARAM_FILE=$(find "$SWEEP_DIR" -name "parameters.csv" -type f | head -1)

if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: No parameters.csv found in $SWEEP_DIR"
    echo "Make sure the sweep completed successfully."
    exit 1
fi

echo "Reading configuration from:"
echo "  $PARAM_FILE"
echo ""

# Extract fixed parameters (these are the same across all runs in the sweep)
INITIAL_RADIUS=$(grep "^initial_tumor_radius," "$PARAM_FILE" | cut -d',' -f2)
VESSEL_CONFIG=$(grep "^vessel_density_config," "$PARAM_FILE" | cut -d',' -f2)
HYPOXIA_DEV=$(grep "^hypoxia_dev_days," "$PARAM_FILE" | cut -d',' -f2 2>/dev/null || echo "0")
ALPHA_HYPOXIC=$(grep "^alpha_hypoxic," "$PARAM_FILE" | cut -d',' -f2)
BETA_HYPOXIC=$(grep "^beta_hypoxic," "$PARAM_FILE" | cut -d',' -f2)

echo "Configuration from sweep:"
echo "  Initial radius: ${INITIAL_RADIUS} um"
echo "  Vessel config: ${VESSEL_CONFIG}"
echo "  Hypoxia development: ${HYPOXIA_DEV} days"
echo "  Hypoxic α/β: ${ALPHA_HYPOXIC} / ${BETA_HYPOXIC}"
echo ""

# =============================================================================
# EXTRACT SWEEP PARAMETER RANGES
# =============================================================================

# Get all unique intervals and skews from directory names
INTERVALS=($(find "$SWEEP_DIR" -type d -name "interval_*_skew_*_rep_*" | \
             sed 's/.*interval_\([0-9]*\)_skew.*/\1/' | sort -nu))

SKEWS=($(find "$SWEEP_DIR" -type d -name "interval_*_skew_*_rep_*" | \
         sed 's/.*skew_\([-0-9]*\)_rep.*/\1/' | sort -nu))

echo "Sweep parameter ranges found:"
echo "  Intervals: ${INTERVALS[@]}"
echo "  Skews: ${SKEWS[@]}"
echo ""

# =============================================================================
# USER SELECTS ROW OR COLUMN TO DEBUG
# =============================================================================

echo "Select debug mode:"
echo "  1) Row (fix skew, vary interval)"
echo "  2) Column (fix interval, vary skew)"
echo ""
read -p "Choice [1-2]: " MODE_CHOICE

if [ "$MODE_CHOICE" = "1" ]; then
    # Row mode - select which skew to fix
    echo ""
    echo "Available skews: ${SKEWS[@]}"
    read -p "Select skew value: " FIXED_SKEW
    
    SWEEP_VALUES=("${INTERVALS[@]}")
    SWEEP_PARAM="interval"
    FIXED_PARAM="skew"
    FIXED_VALUE=$FIXED_SKEW
    
    SWEEP_NAME="interval_row_skew${FIXED_SKEW}"
    OUTPUT_BASE="results/debug_sweeps/${SWEEP_NAME}_$(date +%Y%m%d_%H%M%S)"
    
elif [ "$MODE_CHOICE" = "2" ]; then
    # Column mode - select which interval to fix
    echo ""
    echo "Available intervals: ${INTERVALS[@]}"
    read -p "Select interval value: " FIXED_INTERVAL
    
    SWEEP_VALUES=("${SKEWS[@]}")
    SWEEP_PARAM="skew"
    FIXED_PARAM="interval"
    FIXED_VALUE=$FIXED_INTERVAL
    
    SWEEP_NAME="interval${FIXED_INTERVAL}_skew_col"
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
echo "Configuration file: $PARAM_FILE" >> "${LOGFILE}"
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
        
        # Set interval and skew based on mode
        if [ "$MODE_CHOICE" = "1" ]; then
            INTERVAL=$VALUE
            SKEW=$FIXED_VALUE
        else
            INTERVAL=$FIXED_VALUE
            SKEW=$VALUE
        fi
        
        echo "[${CURRENT_RUN}/${TOTAL_RUNS}] Running: interval=${INTERVAL}, skew=${SKEW}, replicate=${REP}"
        
        START_TIME=$(date +%s)
        SEED=$((42 + REP))  # Seeds: 43, 44, 45, 46, 47 for reps 1-5
        
        # Run with parameter file so DebugRunner loads fixed params
        ./gradlew runDebug --args="interval=${INTERVAL} skew=${SKEW} seed=${SEED} paramfile=${PARAM_FILE}" \
            2>&1 | tee -a "${LOGFILE}"
        
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        
        # Find most recent debug output directory
        LATEST_DIR=$(ls -td results/debug_runs/Debug_I${INTERVAL}_S${SKEW}_* 2>/dev/null | head -1)
        
        if [ -z "${LATEST_DIR}" ]; then
            echo "ERROR: Could not find output directory"
            echo "ERROR: interval=${INTERVAL}, skew=${SKEW}, rep=${REP} - no output" >> "${LOGFILE}"
            continue
        fi
        
        # Move to organized location
        TARGET="${OUTPUT_BASE}/interval_${INTERVAL}_skew_${SKEW}_rep_${REP}"
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
echo "To visualize, run:"
if [ "$MODE_CHOICE" = "1" ]; then
    echo "  python compare_sweep_sections.py --mode row --value ${FIXED_VALUE}"
else
    echo "  python compare_sweep_sections.py --mode col --value ${FIXED_VALUE}"
fi
echo ""