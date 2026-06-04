#!/bin/bash

# Automatically locate and source the shared utility script
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
UTILS_FILE="run_utils.sh"

while [ "$CURRENT_DIR" != "/Verification" ]; do
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

# Script configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


GPYRO_CASE="./sensible_enthalpy.data"

POST_TRAITEMENT="$SCRIPT_DIR/plot_results.py"

export OMP_NUM_THREADS=1


# Run main simulation
cd "$SCRIPT_DIR"
echo "----LAUNCHING GPYRO----"
"$GPYRO" "$GPYRO_CASE"
echo ""

# Post-processing
echo "----POST-TREATING RESULTS----"
"$PYTHON_INTERPRETER" "$POST_TRAITEMENT"
POST_EXIT_CODE=$?
echo ""

# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
exit $POST_EXIT_CODE




