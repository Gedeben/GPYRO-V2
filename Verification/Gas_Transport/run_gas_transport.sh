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

LATEX_FILE="verification_gas_transport"
LATEX_CMD="\\def\\GasTransport{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"

CAS_A="$SCRIPT_DIR/Instantaneous_gas_release"
CAS_B="$SCRIPT_DIR/darcy_Pfixed"
CAS_C="$SCRIPT_DIR/darcy_mfixed"
CAS_D="$SCRIPT_DIR/Gas_Diffusion"
CAS_E="$SCRIPT_DIR/YJ_Advection"
CAS_F="$SCRIPT_DIR/Thermal_equilibrium"

echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#              GAS TRANSPORT CASES             #"
echo "#                                              #"
echo "################################################"

echo ""
echo ""
echo ""
echo "====================================="
echo "-----Instantaneous gas release-------"
echo "====================================="
run_case "$CAS_A" "run.sh"
echo "END OF: Instantaneous gas release ---"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "-----Darcy with fixed Pressure-------"
echo "====================================="
run_case "$CAS_B" "run.sh"
echo "END OF: Darcy with fixed Pressure ---"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "-----Darcy with Imposed Flux --------"
echo "====================================="
run_case "$CAS_C" "run.sh"
echo "END OF: Darcy with Imposed Flux -----"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "------ Gas Species Diffusion --------"
echo "====================================="
run_case "$CAS_D" "run.sh"
echo "END OF: Gas Species Diffusion -------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "------ Gas Species Advection --------"
echo "====================================="
run_case "$CAS_E" "run.sh"
echo "END OF: Gas Species Advetion -------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "---- Gas-Solid Thermal exchange -----"
echo "====================================="
run_case "$CAS_F" "run.sh"
echo "END OF: Gas-Solid Thermal exchange --"
echo "====================================="


# PDF report generation
compile_pdf_report "simple"

echo ""
echo ""
echo ""
echo "################################################"
echo "#         DONE ALL GAS TRANSPORT CASES         #"
echo "################################################"
echo ""
exit 0
