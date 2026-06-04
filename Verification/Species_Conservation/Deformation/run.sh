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

CAS_A="$SCRIPT_DIR/deformation_pure_condesed"
CAS_B="$SCRIPT_DIR/no_deformation"
CAS_C="$SCRIPT_DIR/shrinking"
CAS_D="$SCRIPT_DIR/swelling"
CAS_E="$SCRIPT_DIR/non_charring"
LATEX_FILE="deformation"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"
POST_TRAITEMENT="$SCRIPT_DIR/plot_results.py"
export OMP_NUM_THREADS=1


echo "----LAUNCHING GPYRO----"
echo ">> Running sub-case 1/5"
run_case "$CAS_A" "run.sh"
echo ""

echo ">> Running sub-case 2/5"
run_case "$CAS_B" "run.sh"
echo ""

echo ">> Running sub-case 3/5"
run_case "$CAS_C" "run.sh"
echo ""

echo ">> Running sub-case 4/5"
run_case "$CAS_D" "run.sh"
echo ""

echo ">> Running sub-case 5/5"
run_case "$CAS_E" "run.sh"
echo ""

# Post-processing
echo "----POST-TREATING RESULTS----"
"$PYTHON_INTERPRETER" "$POST_TRAITEMENT"
POST_EXIT_CODE=$?
echo ""


# Optional: PDF report generation
compile_pdf_report "simple"


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
exit $POST_EXIT_CODE



