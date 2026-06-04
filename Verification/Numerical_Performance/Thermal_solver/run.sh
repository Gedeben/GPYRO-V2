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
RUN_PDF=true

# Script configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


DIR_A="$SCRIPT_DIR/Spatial"
DIR_B="$SCRIPT_DIR/Temporal"

LATEX_FILE="Thermal_solver_convergence"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"
export OMP_NUM_THREADS=1

echo "----LAUNCHING GPYRO----"

echo ">> Running Spatial convergence"
run_case "$DIR_A" "run.sh"
RET_A=$?

echo ">> Running Temporal convergence"
run_case "$DIR_B" "run.sh"
RET_B=$?

# Check results
if [ $RET_A -eq 0 ] && [ $RET_B -eq 0 ]; then
    POST_EXIT_CODE=0
else
    POST_EXIT_CODE=1
fi

# Optional: PDF report generation
compile_pdf_report "simple"

exit $POST_EXIT_CODE
