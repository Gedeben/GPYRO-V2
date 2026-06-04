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

ARGS=()

if [ "$RUN_ANALYTICAL" = true ]; then
  ARGS+=(--compute-analytical)
fi

if [ "$RUN_PDF" = false ]; then
  ARGS+=( --no-pdf-report)
fi

LATEX_FILE="verification_report_species_conservation"
LATEX_CMD="\\def\\SpecsCons{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"

CAS_B="$SCRIPT_DIR/Simple_TGA"
CAS_C="$SCRIPT_DIR/Deformation"
CAS_D="$SCRIPT_DIR/Drying"
CAS_E="$SCRIPT_DIR/Mixture_Law"
CAS_F="$SCRIPT_DIR/T_Dependant_Variables"
CAS_G="$SCRIPT_DIR/Oxydation"
CAS_H="$SCRIPT_DIR/Reaction_Enthalpy"

echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#          SPECIES CONSERVATION CASES          #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "-------------Simple TGA--------------"
echo "====================================="
run_case "$CAS_B" "run.sh"
echo "END OF: ATG PMMA ---------------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-------------Deformation-------------"
echo "====================================="
run_case "$CAS_C" "run.sh"
echo "END OF: Deformation -----------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "--------------Drying 0D--------------"
echo "====================================="
run_case "$CAS_D" "run.sh"
echo "END OF: Drying 0D -------------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-------------Mixture Law-------------"
echo "====================================="
run_case "$CAS_E" "run.sh"
echo "END OF: Mixture Law -----------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "------T Dependant Variables----------"
echo "====================================="
run_case "$CAS_F" "run.sh"
echo "END OF: T Dependant Variables -------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-------------Oxydation---------------"
echo "====================================="
run_case "$CAS_G" "run.sh"
echo "END OF: Oxydation -------------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "---------Reaction Enthalpy-----------"
echo "====================================="
run_case "$CAS_H" "run.sh"
echo "END OF: Reaction Enthalpy -----------"
echo "====================================="

# PDF report generation
compile_pdf_report "simple"

echo ""
echo ""
echo ""
echo "################################################"
echo "#     DONE ALL SPECIES CONSERVATION CASES      #"
echo "################################################"
echo ""
exit 0
