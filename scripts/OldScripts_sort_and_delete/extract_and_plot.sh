#!/bin/bash
# Extract CSV data from simulation output and generate plots

if [ $# -lt 1 ]; then
    echo "Usage: $0 <simulation_output.txt> [output_prefix]"
    echo ""
    echo "Extracts CSV_DATA lines from simulation output and creates plots"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_PREFIX="${2:-lookup_table}"

# Extract CSV data
echo "Extracting CSV data from $INPUT_FILE..."
grep "^CSV_DATA:" "$INPUT_FILE" | sed 's/^CSV_DATA://' > "${OUTPUT_PREFIX}.csv"

# Count how many data points
NUM_LINES=$(wc -l < "${OUTPUT_PREFIX}.csv")
echo "Extracted $NUM_LINES data lines"

if [ "$NUM_LINES" -eq 0 ]; then
    echo "No CSV data found in output!"
    echo "Make sure the simulation output contains CSV_DATA: lines"
    exit 1
fi

# Generate plots using Python
echo "Generating plots..."
python3 scripts/plot_lookup_table.py "${OUTPUT_PREFIX}.csv" "$OUTPUT_PREFIX"

echo ""
echo "Done! Check these files:"
echo "  - ${OUTPUT_PREFIX}.csv (extracted data)"
echo "  - ${OUTPUT_PREFIX}_lookup.png (4-panel plot)"
echo "  - ${OUTPUT_PREFIX}_comparison.png (comparison plot)"
