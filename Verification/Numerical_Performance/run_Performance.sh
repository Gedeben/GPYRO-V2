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

LATEX_FILE="verification_numerical_performance"
LATEX_CMD="\\def\\NumPerf{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"



CAS_A="$SCRIPT_DIR/Species_solver"
CAS_B="$SCRIPT_DIR/Thermal_solver"
CAS_C="$SCRIPT_DIR/OMP_parallelization"
CAS_D="$SCRIPT_DIR/Cone_Calorimeter"



echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#             NUMERICAL PERFORMACE             #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "-----------SPECIES SOLVER------------"
echo "====================================="

run_case "$CAS_A" "run.sh"
echo "END OF: SPECIES SOLVER --------------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "----------THERMAL SOLVER 1D----------"
echo "====================================="

run_case "$CAS_B" "run.sh"
echo "END OF: THERMAL SOLVER 1D------------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "--------OMP PARALELLIZATION----------"
echo "====================================="
run_case "$CAS_C" "run.sh"
echo "END OF: OMP PARALELLIZATION ---------"
echo "====================================="
echo ""

echo ""
echo ""
echo ""
echo "====================================="
echo "----CONE CALORIMETER CONVERGENCE-----"
echo "====================================="
run_case "$CAS_D" "run.sh"
echo "END OF:CONE CALORIMETER CONVERGENCE---"
echo "====================================="


# PDF report generation
compile_pdf_report "simple"


echo "################################################"
echo "#    DONE ALL NUMERICAL PERFORMACE CASES       #"
echo "################################################"
echo ""
exit 0
