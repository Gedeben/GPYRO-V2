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

# Check Smokeview availability
check_smokeview_executable

# Script configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

GPYRO_CASE="./Reference_CC_3D_Conv.data"
POST_TRAITEMENT="$SCRIPT_DIR/plot_results.py"
LATEX_FILE="ref_cc_3D_conv"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"

export OMP_NUM_THREADS=4

echo "----LAUNCHING GPYRO----"
(
"$GPYRO" "$GPYRO_CASE"
)

echo "----POST-TREATING RESULTS----"
if [ "$SKIP_SMOKEVIEW" = false ]; then
  "$smokeview" -runscript reference_cc_3D_conv_01 > /dev/null
  "$PYTHON_INTERPRETER" "$POST_TRAITEMENT"
   POST_EXIT_CODE=$?
echo ""

else
  echo "Skipping visual post-processing due to missing or invalid Smokeview."
fi


# Optional: PDF report generation

compile_pdf_report "simple"


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
echo ""
echo ""
exit $POST_EXIT_CODE



