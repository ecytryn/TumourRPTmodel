#!/bin/bash

# Automated Debug Sweep Runner - DOSE-RECEPTOR VERSION
# Runs simulations across a dose or receptor parameter row with multiple replicates
# Intelligently interleaves replicates so you can see trends early

# =============================================================================
# CONFIGURATION - Edit these to match your needs
# =============================================================================

# Parameter space to explore - CHOOSE ONE:

# OPTION A: Dose row (fixed receptor density)
DOSES=(50 75 100 125 150 175 200 225 250)  # nmol
RECEPTOR=1.0                                 # Fixed at baseline (100%)

# OPTION B: Receptor column (fixed dose) - UNCOMMENT TO USE
# DOSES=(150)                                # Fixed dose in nmol
# RECEPTORS=(0.76 0.84 0.92 1.0 1.08 1.16 1.24 1.32 1.4 1.48)  # Multipliers of baseline

# Number of replicates per point
NUM_REPLICATES=5

# Output organization
SWEEP_NAME="dose_row_receptor1.0"
OUTPUT_BASE="results/debug_sweeps/${SWEEP_NAME}_$(date +%Y%m%d_%H%M%S)"

# Specific days to extract images
IMAGE_DAYS=(0 5 10 20 30 40 50 60 80 100)

# =============================================================================
# SETUP
# =============================================================================

echo "============================================================"
echo "  AUTOMATED DEBUG SWEEP: ${SWEEP_NAME}"
echo "============================================================"
echo ""
echo "Configuration:"

# Determine sweep type
if [ ${#DOSES[@]} -gt 1 ]; then
    echo "  Sweep type: DOSE ROW"
    echo "  Doses: ${DOSES[@]} nmol"
    echo "  Receptor: ${RECEPTOR} (fixed)"
    PARAM_ARRAY=("${DOSES[@]}")
    SWEEP_TYPE="dose"
elif [ ${#RECEPTORS[@]} -gt 1 ]; then
    echo "  Sweep type: RECEPTOR COLUMN"
    echo "  Dose: ${DOSES[0]} nmol (fixed)"
    echo "  Receptors: ${RECEPTORS[@]}"
    PARAM_ARRAY=("${RECEPTORS[@]}")
    SWEEP_TYPE="receptor"
else
    echo "ERROR: Must have either multiple doses or multiple receptors"
    exit 1
fi

echo "  Replicates: ${NUM_REPLICATES}"
echo "  Output: ${OUTPUT_BASE}"
echo ""
echo "This will run $(( ${#PARAM_ARRAY[@]} * ${NUM_REPLICATES} )) simulations"
echo "Estimated time: ~$(( ${#PARAM_ARRAY[@]} * ${NUM_REPLICATES} * 10 )) minutes (10 min/sim)"
echo ""
read -p "Press Enter to continue or Ctrl+C to cancel..."
echo ""

# Create output directory
mkdir -p "${OUTPUT_BASE}"

# Log file
LOGFILE="${OUTPUT_BASE}/sweep_log.txt"
echo "Sweep started: $(date)" > "${LOGFILE}"
if [ "$SWEEP_TYPE" = "dose" ]; then
    echo "Configuration: doses=${DOSES[@]}, receptor=${RECEPTOR}, replicates=${NUM_REPLICATES}" >> "${LOGFILE}"
else
    echo "Configuration: dose=${DOSES[0]}, receptors=${RECEPTORS[@]}, replicates=${NUM_REPLICATES}" >> "${LOGFILE}"
fi
echo "" >> "${LOGFILE}"

# =============================================================================
# RUN SIMULATIONS
# =============================================================================

echo "Starting simulations..."
echo "Strategy: Run all parameter values for replicate 1, then rep 2, etc."
echo "This gives you full parameter coverage quickly!"
echo ""

TOTAL_RUNS=$(( ${#PARAM_ARRAY[@]} * ${NUM_REPLICATES} ))
CURRENT_RUN=0

for REP in $(seq 1 ${NUM_REPLICATES}); do
    echo "=========================================="
    echo "  REPLICATE ${REP}/${NUM_REPLICATES}"
    echo "=========================================="
    echo ""
    
    for PARAM in "${PARAM_ARRAY[@]}"; do
        CURRENT_RUN=$((CURRENT_RUN + 1))
        
        # Set parameters based on sweep type
        if [ "$SWEEP_TYPE" = "dose" ]; then
            DOSE=$PARAM
            REC=$RECEPTOR
            echo "[${CURRENT_RUN}/${TOTAL_RUNS}] Running: dose=${DOSE} nmol, receptor=${REC}, replicate=${REP}"
        else
            DOSE=${DOSES[0]}
            REC=$PARAM
            echo "[${CURRENT_RUN}/${TOTAL_RUNS}] Running: dose=${DOSE} nmol, receptor=${REC}, replicate=${REP}"
        fi
        
        # Run the simulation
        START_TIME=$(date +%s)
        
        SEED=$((42 + REP))  # Seeds: 43, 44, 45, 46, 47 for reps 1-5
        ./gradlew runDebug --args="dose=${DOSE} receptors=${REC} seed=${SEED}" 2>&1 | tee -a "${LOGFILE}"
        
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        
        # Find the most recent debug output directory
        # Format: Debug_D150_R1.00_TIMESTAMP
        LATEST_DIR=$(ls -td results/debug_runs/Debug_D*_R*_* 2>/dev/null | head -1)
        
        if [ -z "${LATEST_DIR}" ]; then
            echo "ERROR: Could not find output directory for dose=${DOSE}, receptor=${REC}"
            echo "ERROR: dose=${DOSE}, receptor=${REC}, rep=${REP} - no output found" >> "${LOGFILE}"
            continue
        fi
        
        # Organize output: move to sweep directory with clear naming
        if [ "$SWEEP_TYPE" = "dose" ]; then
            TARGET_DIR="${OUTPUT_BASE}/dose_${DOSE}_rep_${REP}"
        else
            TARGET_DIR="${OUTPUT_BASE}/receptor_${REC}_rep_${REP}"
        fi
        
        mv "${LATEST_DIR}" "${TARGET_DIR}"
        
        echo "  -> Saved to: ${TARGET_DIR}"
        echo "  -> Runtime: ${ELAPSED}s"
        echo ""
        echo "Run ${CURRENT_RUN}/${TOTAL_RUNS}: dose=${DOSE}, receptor=${REC}, rep=${REP}, time=${ELAPSED}s, dir=${TARGET_DIR}" >> "${LOGFILE}"
    done
    
    echo ""
    echo "Replicate ${REP} complete! You can now review results while next replicate runs."
    echo ""
    
    # Automatically generate comparison grid for completed replicates
#    echo "Generating comparison grid for replicates 1-${REP}..."
#    python3 create_image_grid.py "${OUTPUT_BASE}" --output "${OUTPUT_BASE}/comparison_grid_rep1-${REP}.pdf" 2>&1 | grep -E "(ERROR|saved|Page added)"    
#    if [ -f "${OUTPUT_BASE}/comparison_grid_rep1-${REP}.pdf" ]; then
#        echo "  -> Grid saved: ${OUTPUT_BASE}/comparison_grid_rep1-${REP}.pdf"
#        echo "  -> Open this file to see results so far!"
#    fi
#    echo ""
done

echo ""
echo "============================================================"
echo "  ALL SIMULATIONS COMPLETE"
echo "============================================================"
echo "Sweep finished: $(date)" >> "${LOGFILE}"
echo ""
echo "Results saved in: ${OUTPUT_BASE}"
echo ""

# =============================================================================
# GENERATE SUMMARY
# =============================================================================

echo "Generating summary..."

SUMMARY_FILE="${OUTPUT_BASE}/sweep_summary.txt"

echo "Debug Sweep Summary" > "${SUMMARY_FILE}"
echo "===================" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Sweep: ${SWEEP_NAME}" >> "${SUMMARY_FILE}"
echo "Date: $(date)" >> "${SUMMARY_FILE}"

if [ "$SWEEP_TYPE" = "dose" ]; then
    echo "Doses: ${DOSES[@]} nmol" >> "${SUMMARY_FILE}"
    echo "Receptor: ${RECEPTOR} (fixed)" >> "${SUMMARY_FILE}"
else
    echo "Dose: ${DOSES[0]} nmol (fixed)" >> "${SUMMARY_FILE}"
    echo "Receptors: ${RECEPTORS[@]}" >> "${SUMMARY_FILE}"
fi

echo "Replicates: ${NUM_REPLICATES}" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Results:" >> "${SUMMARY_FILE}"
echo "--------" >> "${SUMMARY_FILE}"

for PARAM in "${PARAM_ARRAY[@]}"; do
    echo "" >> "${SUMMARY_FILE}"
    
    if [ "$SWEEP_TYPE" = "dose" ]; then
        echo "Dose = ${PARAM} nmol:" >> "${SUMMARY_FILE}"
        PREFIX="dose_${PARAM}"
    else
        echo "Receptor = ${PARAM}:" >> "${SUMMARY_FILE}"
        PREFIX="receptor_${PARAM}"
    fi
    
    for REP in $(seq 1 ${NUM_REPLICATES}); do
        DIR="${OUTPUT_BASE}/${PREFIX}_rep_${REP}"
        
        if [ -f "${DIR}/populations.csv" ]; then
            # Get final tumor count (sum of normoxic + hypoxic)
            FINAL_LINE=$(tail -1 "${DIR}/populations.csv")
            NORMOXIC=$(echo "$FINAL_LINE" | cut -d',' -f2)
            HYPOXIC=$(echo "$FINAL_LINE" | cut -d',' -f3)
            FINAL_POP=$(echo "$NORMOXIC + $HYPOXIC" | bc)
            
            # Determine outcome
            if (( $(echo "$FINAL_POP < 10" | bc -l) )); then
                OUTCOME="CURE"
            else
                OUTCOME="FAILURE"
            fi
            
            echo "  Rep ${REP}: ${OUTCOME} (final cells: ${FINAL_POP})" >> "${SUMMARY_FILE}"
        else
            echo "  Rep ${REP}: ERROR - no data" >> "${SUMMARY_FILE}"
        fi
    done
done

cat "${SUMMARY_FILE}"

echo ""
echo "Summary saved to: ${SUMMARY_FILE}"
echo ""

# =============================================================================
# NEXT STEPS
# =============================================================================

echo "============================================================"
echo "  NEXT: Create Image Comparison Grid"
echo "============================================================"
echo ""
echo "To create a visual comparison of all simulations, run:"
echo ""
echo "  python create_image_grid.py ${OUTPUT_BASE}"
echo ""
echo "This will extract key timepoints and create a comparison PDF."
echo ""
