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

ARGS=()

if [ "$RUN_ANALYTICAL" = true ]; then
  ARGS+=(--compute-analytical)
fi

if [ "$RUN_PDF" = false ]; then
  ARGS+=( --no-pdf-report)
fi

LATEX_FILE="verification_report_thermal_2D_3D"
LATEX_CMD="\\def\\ThermalTreeD{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"



CAS_A="$SCRIPT_DIR/Thermal_equilibrium_3D"
CAS_B="$SCRIPT_DIR/I_beam"



echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#             THERMAL 2D/3D CASES              #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "-------THERMAL EQUILBRIUM 3D---------"
echo "====================================="
run_case "$CAS_A" "run.sh"
echo "END OF: THERMAL EQUILBRIUM 3D--------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-------------- I BEAM ---------------"
echo "====================================="
run_case "$CAS_B" "run.sh"
echo "END OF:I BEAM -----------------------"
echo "====================================="


# PDF report generation
compile_pdf_report "simple"


echo "################################################"
echo "#      DONE ALL THERMAL 2D/3D CASES            #"
echo "################################################"
echo ""
exit 0
