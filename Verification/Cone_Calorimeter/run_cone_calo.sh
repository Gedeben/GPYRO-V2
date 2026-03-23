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

if [ "$FULL_RUN" = true ]; then
  ARGS+=( --full)
fi

# Check Smokeview availability
check_smokeview_executable

LATEX_FILE="verification_report_cone_calo"
LATEX_CMD="\\def\\ConeCalo{}  \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"


CAS_A="$SCRIPT_DIR/Reference_Cone_Calorimmeter"
CAS_B="$SCRIPT_DIR/Time_Step_Convergence"
CAS_C="$SCRIPT_DIR/Spatial_Convergence"
CAS_D="$SCRIPT_DIR/Extention_Cone_Calo_2D_3D"
CAS_E="$SCRIPT_DIR/3D_Lateral_convection"
CAS_F="$SCRIPT_DIR/OMP_parallelization"

echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#            CONE CALORIMETER CASES            #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "----REFERENCE CONE CALORIMETER-------"
echo "====================================="
run_case "$CAS_A" "run.sh"
echo "END OF: REFERENCE CONE CALORIMETER---"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "--------TIME STEP CONVERGENCE--------"
echo "====================================="
run_case "$CAS_B" "run.sh"
echo "END OF: TIME STEP CONVERGENCE--------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "---------SPATIAL CONVERGENCE---------"
echo "====================================="
run_case "$CAS_C" "run.sh"
echo "END OF: SPATIAL CONVERGENCE----------"
echo "====================================="

echo ""
echo ""
echo ""
echo "====================================="
echo "-----------EXTENDED 2D/3D------------"
echo "====================================="
run_case "$CAS_D" "run.sh"
echo "END OF: EXTENDED 2D/3D --------------"
echo "====================================="


echo ""
echo ""
echo ""
echo "====================================="
echo "-------3D LATERAL CONVECTION---------"
echo "====================================="
run_case "$CAS_E" "run.sh"
echo "END OF: 3D LATERAL CONVECTION -------"
echo "====================================="
echo ""


echo ""
echo ""
echo ""
echo "====================================="
echo "--------OMP PARALELLIZATION----------"
echo "====================================="
run_case "$CAS_F" "run.sh"
echo "END OF: OMP PARALELLIZATION ---------"
echo "====================================="
echo ""






# PDF report generation
compile_pdf_report "bibtex"



echo "################################################"
echo "#        DONE ALL CONE CALORIMETER CASES       #"
echo "################################################"
echo ""
exit 0















