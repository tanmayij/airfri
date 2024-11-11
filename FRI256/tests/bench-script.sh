#!/bin/bash

total_time=0
num_runs=100
successful_runs=0

echo "Running the program $num_runs times..."

for ((i=1; i<=num_runs; i++))
do
    echo "Running iteration $i/$num_runs..."

    # Run the program and capture its output
    output=$(./fri_with_merkle_gpu 2>&1)

    # Debug: Print the full output to check what is being captured
    echo "Output from iteration $i:"
    echo "$output"

    # Extract time using grep and awk. Adjust the regex if needed.
    time_taken=$(echo "$output" | grep -oP 'fri took \K[0-9.]+(?= seconds)')

    # Check if the time was found
    if [[ ! -z "$time_taken" ]]; then
        echo "Time for iteration $i: $time_taken seconds"
        total_time=$(echo "$total_time + $time_taken" | bc)
        ((successful_runs++))
    else
        echo "Error in iteration $i: Time not found or program failed."
    fi
done

# Calculate average time if there are successful runs
if [[ $successful_runs -gt 0 ]]; then
    average_time=$(echo "scale=6; $total_time / $successful_runs" | bc)
    echo "Average time for $successful_runs successful runs: $average_time seconds"
else
    echo "No successful runs were recorded."
fi