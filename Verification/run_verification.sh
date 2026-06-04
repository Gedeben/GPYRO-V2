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

cd $SCRIPT_DIR
rm -f status_report.txt



echo ""
echo ""
echo "==============================================================================="
echo "*************************  RUNNING GPYRO VERIFICATION  ************************"
echo "==============================================================================="

printf "%-25s : %s\n" "GPYRO executable"        "$GPYRO"
printf "%-25s : %s\n" "Smokeview executable"          "$smokeview"
printf "%-25s : %s\n" "Python Interpreter"       "$PYTHON_INTERPRETER"
printf "%-25s : %s\n" "LaTeX Interpreter"        "$LATEX_INTERPRETER"

echo "-------------------------------------------------------------------------------"
printf "%-40s : %s\n" "Re-compute analytical solution"             "$RUN_ANALYTICAL"
printf "%-40s : %s\n" "Generate PDF report"                        "$RUN_PDF"
printf "%-40s : %s\n" "Run all time-expensive cases"               "$FULL_RUN"
echo "==============================================================================="
echo ""


#"#=================================================#"
# "------------------- Validation -------------------"
# "#================================================#"
run_case "$SCRIPT_DIR/Validation" "run_validation.sh"



#"#=================================================#"
# "----------- Verification Thermal 1D --------------"
# "#================================================#"
run_case "$SCRIPT_DIR/Thermal_1D" "run_thermal_1D.sh"


#"#=================================================#"
# "-------- Verification Thermal 2D/3D --------------"
# "#================================================#"
run_case "$SCRIPT_DIR/Thermal_2D_3D" "run_thermal_2D_3D.sh"

#"#================================================#"
#"------- Verification Species Conservation --------"
#"#================================================#"
run_case "$SCRIPT_DIR/Species_Conservation" "run_Species_Conservation.sh"


#"#================================================#"
# "----------- Verification Multiphyiscs -----------"
# "#================================================#"
run_case "$SCRIPT_DIR/Multiphyiscs" "run_multiphyics.sh"



#"#================================================#"
# "--------- Verification Gas Transport ----------"
# "#================================================#"
run_case "$SCRIPT_DIR/Gas_Transport" "run_gas_transport.sh"


#"#================================================#"
# "------ Verification Numerical Performance --------"
# "#================================================#"
run_case "$SCRIPT_DIR/Numerical_Performance" "run_Performance.sh"

check_status_report "$SCRIPT_DIR/status_report.txt"

# PDF report generation
LATEX_FILE="verification_report"
LATEX_CMD="\\def\\Verification{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"

compile_pdf_report "book"




















