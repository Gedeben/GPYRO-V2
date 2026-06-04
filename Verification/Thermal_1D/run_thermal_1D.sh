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

LATEX_FILE="verification_report_thermal_1D"
LATEX_CMD="\\def\\ThermalOneD{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"



CAS_A="$SCRIPT_DIR/convective_cooling"
CAS_B="$SCRIPT_DIR/heat_conduction"
CAS_C="$SCRIPT_DIR/insulated_steel_plate"
CAS_D="$SCRIPT_DIR/heat_conduction_kc"
CAS_E="$SCRIPT_DIR/radiation_convection"
CAS_F="$SCRIPT_DIR/radiation_loss"
CAS_G="$SCRIPT_DIR/fixed_temp_sides"
CAS_H="$SCRIPT_DIR/radiation_1d"
CAS_I="$SCRIPT_DIR/radiation_depth"
CAS_J="$SCRIPT_DIR/rad_conv_depth"


echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#               THERMAL 1D CASES               #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "---------CONVECTIVE COOLING----------"
echo "====================================="

run_case "$CAS_A" "run.sh"
echo "END OF: CONVECTIVE COOLING-----------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "-----------HEAT CONDUCTION-----------"
echo "====================================="

run_case "$CAS_B" "run.sh"
echo "END OF: HEAT CONDUCTION--------------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "--------INSULATED STEEL PLATE--------"
echo "====================================="

run_case "$CAS_C" "run.sh"
echo "END OF: INSULATED STEEL PLATE -------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "------HEAT CONDUCTION K,C (T)--------"
echo "====================================="

run_case "$CAS_D" "run.sh"
echo "END OF: HEAT CONDUCTION K,C (T)------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "--------RADIATION CONVECTION---------"
echo "====================================="

run_case "$CAS_E" "run.sh"
echo "END OF: RADIATION CONVECTION---------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "--------RADIATION LOSS---------------"
echo "====================================="
run_case "$CAS_F" "run.sh"
echo "END OF: RADIATION LOSS---------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-----FIXED BOUNDARY TEMPERATURES-----"
echo "====================================="
run_case "$CAS_G" "run.sh"
echo "END OF: FIXED BOUNDARY TEMPERATURES--"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "----------SURFACE RADIATION----------"
echo "====================================="
run_case "$CAS_H" "run.sh"
echo "END OF: SURFACE RADIATION------------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "--------RADIATION ABSORPTION---------"
echo "====================================="
run_case "$CAS_I" "run.sh"
echo "END OF: RADIATION ABSORPTION---------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-----RADIATION ABSORPTION + CONV-----"
echo "====================================="
run_case "$CAS_J" "run.sh"
echo "END OF: RADIATION ABSORPTION + CONV--"
echo "====================================="

# PDF report generation
compile_pdf_report "bibtex"


echo "################################################"
echo "#        DONE ALL THERMAL 1D CASES             #"
echo "################################################"
echo ""
exit 0
