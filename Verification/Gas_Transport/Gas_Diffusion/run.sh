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


LATEX_FILE="gas_species_diffusion"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"


CAS_A="$SCRIPT_DIR/hm_01"
CAS_B="$SCRIPT_DIR/hm_0001"



echo ""
echo ""
echo ""
echo "----LAUNCHING GPYRO----"
echo ">> Runnig sub-case 1/2"
run_case "$CAS_A" "run.sh"

echo ""
echo ""
echo ""
echo ">> Runnig sub-case 2/2"

run_case "$CAS_B" "run.sh"




# PDF report generation
compile_pdf_report "simple"

exit 0
