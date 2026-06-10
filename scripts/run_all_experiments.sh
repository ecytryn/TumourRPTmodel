#!/bin/bash

# Array of experiment names
experiments=(
    "WatchGrow"
    "NormoxicSmall"
    "NormoxicMedium"
    "NormoxicLarge"
#    "HypoxicSmall"
#    "HypoxicMedium"
    "HypoxicLarge"
)

# Run each experiment
for exp in "${experiments[@]}"; do
    echo "========================================="
    echo "Running experiment: $exp"
    echo "========================================="
    ./gradlew run --args="$exp"
    
    # Check if successful
    if [ $? -eq 0 ]; then
        echo "✓ $exp completed successfully"
    else
        echo "✗ $exp failed!"
        exit 1
    fi
    echo ""
done

echo "All experiments complete!"