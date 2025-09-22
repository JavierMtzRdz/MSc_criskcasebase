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
# Define static arguments that will be passed to the R script on every iteration.
# These could be file paths, configuration settings, etc.
P_ARG=120
NUM_TRUE=20

# Create the output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# --- Loop Execution ---
# Outer loop to iterate through settings from 1 to 5.
for s in $(seq 1 5)
do
  echo "#################################"
  echo "### Processing for SETTING: $s ###"
  echo "#################################"

  # Inner 'for' loop to iterate from START_INDEX to END_INDEX for each setting.
  for i in $(seq $START_INDEX $END_INDEX)
  do
    # Print which iteration is currently running.
    echo "---------------------------------"
    echo "Running iteration with index: $i (Setting: $s)"
    echo "---------------------------------"

    # Execute the R script, passing the index and the other configured arguments.
    # The order of arguments matters and must match what the R script expects.
    # The setting variable '$s' is now passed instead of a static value.
    Rscript "$R_SCRIPT_NAME" "$P_ARG" "$NUM_TRUE" "$s" "$i"

    # Optional: Check if the R script ran successfully.
    if [ $? -eq 0 ]; then
      echo "Iteration $i completed successfully."
  #  else
  #    echo "Error in iteration $i. Aborting loop."
  #    exit 1 # Exit the shell script with an error code.
    fi
  done
done

echo "---------------------------------"
echo "All iterations completed."
echo "---------------------------------"
