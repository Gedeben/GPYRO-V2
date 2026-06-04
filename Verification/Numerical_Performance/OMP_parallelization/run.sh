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



# List of thread counts to test
THREAD_COUNTS=(1 2 3 4 6 8)

# Detect available number of processors
MAX_PROC=$(grep -c ^processor /proc/cpuinfo)

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
GPYRO_CASE="ref_CC_NZ5000.data"
POST_TRAITEMENT="$SCRIPT_DIR/plot_results.py"
LATEX_FILE="OMP_parallelization"
LATEX_CMD="\\def\\standalone{} \\input{${SCRIPT_DIR}/${LATEX_FILE}.tex}"


echo "----LAUNCHING GPYRO----"

# Loop through thread counts and launch simulations
for NTHREADS in "${THREAD_COUNTS[@]}"; do
  if [ "$NTHREADS" -ge "$MAX_PROC" ]; then
    #echo "Reached system limit: $MAX_PROC cores available. Stopping tests."
    break
  fi

  (DIR="OMP_${NTHREADS}"
  echo ">> Runnig sub-case with $NTHREADS thread(s)"

  mkdir -p "$DIR"
  cp "$GPYRO_CASE" "$DIR/"
  cd $DIR 

  export OMP_NUM_THREADS="$NTHREADS"
  "$GPYRO" "$GPYRO_CASE"
)
done

echo "----POST-TREATING RESULTS----"
"$PYTHON_INTERPRETER" "$POST_TRAITEMENT"
POST_EXIT_CODE=$?


# Optional: PDF report generation
compile_pdf_report "simple"


# Write status report
write_status_report "$SCRIPT_DIR" "$POST_EXIT_CODE"
echo ""
echo ""
exit $POST_EXIT_CODE











