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


DIR_A="$SCRIPT_DIR/heat_conduction_a"
DIR_B="$SCRIPT_DIR/heat_conduction_b"
DIR_C="$SCRIPT_DIR/heat_conduction_c"
DIR_D="$SCRIPT_DIR/heat_conduction_d"
CAS_A="heat_conduction_a.data"
CAS_B="heat_conduction_b.data"
CAS_C="heat_conduction_c.data"
CAS_D="heat_conduction_d.data"
LATEX_FILE="heat_conduction"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"
export OMP_NUM_THREADS=1

echo "----LAUNCHING GPYRO----"
echo ">> Runnig sub-case 1/4"
cd $DIR_A
$GPYRO $CAS_A
cd $SCRIPT_DIR
echo ""

echo ">> Runnig sub-case 2/4"
cd $DIR_B
"$GPYRO" "$CAS_B"
cd $SCRIPT_DIR
echo ""

echo ">> Runnig sub-case 3/4"
cd $DIR_C
"$GPYRO" "$CAS_C"
cd $SCRIPT_DIR/
echo ""

echo ">> Runnig sub-case 4/4"
cd $DIR_D
"$GPYRO" "$CAS_D"
cd $SCRIPT_DIR/
echo ""

# Run analytical solution if not skipped
if [ "$RUN_ANALYTICAL" = true ]; then
  echo "----CALCULATING ANALYTICAL DATA----"
	for DIR in "$DIR_A" "$DIR_B" "$DIR_C" "$DIR_D"; do
	    if [ -f "$DIR/analytical_solution.py" ]; then
		"$PYTHON_INTERPRETER" "$DIR/analytical_solution.py"
	    else
		echo "Warning: No 'analytical_solution.py' found in $DIR"
	    fi
	    echo ""
	done
  "$PYTHON_INTERPRETER" "$ANALYTIQUE_DATA"
  echo ""
fi

echo "----POST-TREATING RESULTS----"
POST_EXIT_CODE=0

for DIR in "$DIR_A" "$DIR_B" "$DIR_C" "$DIR_D"; do
  POST_SCRIPT="$DIR/plot_results.py"
  
  if [ -f "$POST_SCRIPT" ]; then
    "$PYTHON_INTERPRETER" "$POST_SCRIPT"
    EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
      POST_EXIT_CODE=1
    fi
  else
    echo "Warning: No 'plot_results.py' found in $DIR"
    POST_EXIT_CODE=1
  fi
  echo ""
done

# Optional: PDF report generation
compile_pdf_report "simple"


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
exit $POST_EXIT_CODE



