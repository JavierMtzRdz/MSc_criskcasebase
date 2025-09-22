#!/bin/bash

# A shell script to run an R script multiple times with a different index each time.
# This version is enhanced to pass multiple arguments to the R script.

echo "Starting analysis loop..."

# --- Configuration ---
# Define the start and end of your index loop.
START_INDEX=1
END_INDEX=15

# Define the path to your R script.
R_SCRIPT_NAME="paper/code/fit-curve-lambda.R"

# --- Additional Arguments ---
# Define static P_ARG.
P_ARG=1000

# Calculate the dynamic NUM_TRUE value based on a formula derived from P_ARG.
# The formula is: round(P_ARG * 0.04166667) * 4
TEMP_VAL=$(echo "$P_ARG * 0.04166667" | bc)
ROUNDED_VAL=$(printf "%.0f" "$TEMP_VAL")
CALCULATED_NUM_TRUE=$((ROUNDED_VAL * 4))

# Define the array of NUM_TRUE values to be evaluated in the loop.
# It will test the static value 20 and the dynamically calculated value.
NUM_TRUE_VALUES=(20 $CALCULATED_NUM_TRUE)

# --- Loop Execution ---
# Outermost loop to iterate through the specified NUM_TRUE values.
# Note: If P_ARG=120, the calculated value is also 20, so the loop may run twice for the same value.
# This is expected behavior based on the request.
for nt in "${NUM_TRUE_VALUES[@]}"
do
  echo "================================="
  echo "=== Processing for NUM_TRUE: $nt ==="
  echo "================================="

  # Middle loop to iterate through settings from 1 to 5.
  for s in $(seq 1 5)
  do
    echo "#################################"
    echo "### Processing for SETTING: $s ###"
    echo "#################################"

    # Innermost 'for' loop to iterate from START_INDEX to END_INDEX for each setting.
    for i in $(seq $START_INDEX $END_INDEX)
    do
      # Print which iteration is currently running.
      echo "---------------------------------"
      echo "Running iteration with index: $i (NUM_TRUE: $nt, Setting: $s)"
      echo "---------------------------------"

      # Execute the R script, passing the current NUM_TRUE ($nt) from the outer loop.
      # The order of arguments matters and must match what the R script expects.
      Rscript "$R_SCRIPT_NAME" "$P_ARG" "$nt" "$s" "$i"

      # Optional: Check if the R script ran successfully.
      if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully."
    #  else
    #    echo "Error in iteration $i. Aborting loop."
    #    exit 1 # Exit the shell script with an error code.
      fi
    done
  done
done

echo "---------------------------------"
echo "All iterations completed."
echo "---------------------------------"

