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
POST_SCRIPT="plot_results.py"

CAS="rad.data"

TEMPERATURE="analytical_temperature.py"

LATEX_FILE="rad"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"
export OMP_NUM_THREADS=1

echo "----LAUNCHING GPYRO----"
$GPYRO "$CAS"
echo ""

# Run analytical solution if not skipped
if [ "$RUN_ANALYTICAL" = true ]; then
	echo "----CALCULATING ANALYTICAL SOLUTION----"
	"$PYTHON_INTERPRETER" "$TEMPERATURE"
	echo ""
fi

echo "----POST-TREATING RESULTS----"

POST_EXIT_CODE=0
if [ -f "$POST_SCRIPT" ]; then
  "$PYTHON_INTERPRETER" "$POST_SCRIPT"
  EXIT_CODE=$?
  if [ $EXIT_CODE -ne 0 ]; then
    POST_EXIT_CODE=1
  fi
else
  echo "Warning: No 'plot_results.py' found in $SCRIPT_DIR"
  POST_EXIT_CODE=1
fi
echo ""


# Optional: PDF report generation
compile_pdf_report "bibtex"


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
exit $POST_EXIT_CODE



