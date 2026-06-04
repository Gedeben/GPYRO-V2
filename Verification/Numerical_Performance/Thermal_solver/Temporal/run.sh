#!/bin/bash

# Automatically locate and source the shared utility script
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
UTILS_FILE="run_utils.sh"

while [ "$CURRENT_DIR" != "/Verification_Gpyro_V08" ]; do
  if [ -f "$CURRENT_DIR/$UTILS_FILE" ]; then
    source "$CURRENT_DIR/$UTILS_FILE"
    break
  fi
  CURRENT_DIR="$(dirname "$CURRENT_DIR")"
done

if ! declare -f parse_common_args > /dev/null; then
  echo "Error: Unable to locate $UTILS_FILE in Verification."
  exit 1
fi

# Parse arguments and check executable
parse_common_args "$@"
set -- "${PARSED_ARGS[@]}"
check_gpyro_executable "$1"

# Define the expensive cases
SKIP_DIRS=("dt_0001s")


# Script configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

GPYRO_CASE="./radiation_1d.data"
ANALYTIQUE_DATA="$SCRIPT_DIR/analytical_solution.py"
POST_TRAITEMENT="$SCRIPT_DIR/plot_results.py"

export OMP_NUM_THREADS=1



# Run main simulations
# Loop through all subdirectories for mesh refinement studies
echo "---- RUNNING TIME STEP CONVERGENCE CASES ----"
for dir in */ ; do
    # Skip if not a directory with input file
    if [ -d "$dir" ] && [ -f "$dir/$GPYRO_CASE" ]; then

        # Check if this is a costly case and FULL_RUN is not set
        if [[ " ${SKIP_DIRS[@]} " =~ " ${dir%/} " ]] && [ "$FULL_RUN" = false ]; then
            echo ">> Skipping costly case: $dir (use --full to include)"
            continue
        fi

        echo ">> Entering $dir"
        (cd "$dir" && "$GPYRO" "$GPYRO_CASE")
        echo ""
    fi
done

# Optional: analytical solution
if [ "$RUN_ANALYTICAL" = true ]; then
  echo "----CALCULATING ANALYTICAL DATA----"
  "$PYTHON_INTERPRETER" "$ANALYTIQUE_DATA"
  echo ""
fi


# Post-processing
echo "----POST-TREATING RESULTS----"
"$PYTHON_INTERPRETER" "$POST_TRAITEMENT"
POST_EXIT_CODE=$?
echo ""


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
echo ""
echo ""
exit $POST_EXIT_CODE











