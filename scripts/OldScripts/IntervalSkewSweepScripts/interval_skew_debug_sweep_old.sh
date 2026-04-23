#!/bin/bash

# Automated Debug Sweep Runner
# Runs simulations across a parameter row with multiple replicates
# Intelligently interleaves replicates so you can see trends early

# =============================================================================
# CONFIGURATION - Edit these to match your needs
# =============================================================================

# Make sure the lines with ********** right above them are set for your sweep.

# Parameter space to explore
# **********
INTERVALS=(20 22 24 26 28 30)  # Your interval values
#INTERVALS=20               # Fixed interval value (column in heatmap)
# **********
SKEW=5                      # Fixed skew value (row in heatmap)
#SKEW=(-20 -10 0 10 10)      # Your skew values

# Number of replicates per point
# **********
NUM_REPLICATES=5

# Output organization
# **********
SWEEP_NAME="interval_row_skew${SKEW}"
# SWEEP_NAME="interval_${INTERVAL}_skew_column"

OUTPUT_BASE="results/debug_sweeps/${SWEEP_NAME}_$(date +%Y%m%d_%H%M%S)"

# Specific days to extract images (adjust based on your injection schedule)
# For interval row, you probably want: injection day, peak response, endpoint
# **********
IMAGE_DAYS=(0 5 10 15 20 25 30 35 40 45 50 55)

# =============================================================================
# SETUP
# =============================================================================

echo "============================================================"
echo "  AUTOMATED DEBUG SWEEP: ${SWEEP_NAME}"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  Intervals: ${INTERVALS[@]}"
echo "  Skew: ${SKEW} nmol"
echo "  Replicates: ${NUM_REPLICATES}"
echo "  Output: ${OUTPUT_BASE}"
echo ""
echo "This will run $(( ${#INTERVALS[@]} * ${NUM_REPLICATES} )) simulations"
echo "Estimated time: ~$(( ${#INTERVALS[@]} * ${NUM_REPLICATES} * 10 )) minutes (10 min/sim)"
echo ""
read -p "Press Enter to continue or Ctrl+C to cancel..."
echo ""

# Create output directory
mkdir -p "${OUTPUT_BASE}"

# Log file
LOGFILE="${OUTPUT_BASE}/sweep_log.txt"
echo "Sweep started: $(date)" > "${LOGFILE}"
echo "Configuration: intervals=${INTERVALS[@]}, skew=${SKEW}, replicates=${NUM_REPLICATES}" >> "${LOGFILE}"
echo "" >> "${LOGFILE}"

# =============================================================================
# RUN SIMULATIONS
# Strategy: Interleave replicates so you get full row coverage quickly
# =============================================================================

echo "Starting simulations..."
echo "Strategy: Run all intervals for replicate 1, then rep 2, etc."
echo "This gives you full parameter coverage quickly!"
echo ""

TOTAL_RUNS=$(( ${#INTERVALS[@]} * ${NUM_REPLICATES} ))
CURRENT_RUN=0

for REP in $(seq 1 ${NUM_REPLICATES}); do
    echo "=========================================="
    echo "  REPLICATE ${REP}/${NUM_REPLICATES}"
    echo "=========================================="
    echo ""
    
    for INTERVAL in "${INTERVALS[@]}"; do
        CURRENT_RUN=$((CURRENT_RUN + 1))
        
        echo "[${CURRENT_RUN}/${TOTAL_RUNS}] Running: interval=${INTERVAL}, skew=${SKEW}, replicate=${REP}"
        
        # Run the simulation
        START_TIME=$(date +%s)
        
		SEED=$((42 + REP))  # Seed 43, 44, 45, 46, 47 for reps 1-5
		./gradlew runDebug --args="interval=${INTERVAL} skew=${SKEW} seed=${SEED}" 2>&1 | tee -a "${LOGFILE}"
        
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        
        # Find the most recent debug output directory
        LATEST_DIR=$(ls -td results/debug_runs/Debug_I${INTERVAL}_S${SKEW}_* 2>/dev/null | head -1)
        
        if [ -z "${LATEST_DIR}" ]; then
            echo "ERROR: Could not find output directory for interval=${INTERVAL}, skew=${SKEW}"
            echo "ERROR: interval=${INTERVAL}, skew=${SKEW}, rep=${REP} - no output found" >> "${LOGFILE}"
            continue
        fi
        
        # Organize output: move to sweep directory with clear naming
        TARGET_DIR="${OUTPUT_BASE}/interval_${INTERVAL}_rep_${REP}"
        mv "${LATEST_DIR}" "${TARGET_DIR}"
        
        echo "  -> Saved to: ${TARGET_DIR}"
        echo "  -> Runtime: ${ELAPSED}s"
        echo ""
        echo "Run ${CURRENT_RUN}/${TOTAL_RUNS}: interval=${INTERVAL}, skew=${SKEW}, rep=${REP}, time=${ELAPSED}s, dir=${TARGET_DIR}" >> "${LOGFILE}"
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
echo "Intervals: ${INTERVALS[@]}" >> "${SUMMARY_FILE}"
echo "Skew: ${SKEW} nmol" >> "${SUMMARY_FILE}"
echo "Replicates: ${NUM_REPLICATES}" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "Results:" >> "${SUMMARY_FILE}"
echo "--------" >> "${SUMMARY_FILE}"

for INTERVAL in "${INTERVALS[@]}"; do
    echo "" >> "${SUMMARY_FILE}"
    echo "Interval = ${INTERVAL} days:" >> "${SUMMARY_FILE}"
    
    for REP in $(seq 1 ${NUM_REPLICATES}); do
        DIR="${OUTPUT_BASE}/interval_${INTERVAL}_rep_${REP}"
        
        if [ -f "${DIR}/populations.csv" ]; then
            # Get final tumor count
            FINAL_POP=$(tail -1 "${DIR}/populations.csv" | cut -d',' -f2)
            
            # Determine outcome
            if [ "${FINAL_POP}" -lt 10 ]; then
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
