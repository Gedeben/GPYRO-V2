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




CAS_A="$SCRIPT_DIR/MaCFP_PMMA"

echo ""
echo ""
echo ""
echo "################################################"
echo "#          LAUNCHING VERIFICATION OF:          #"
echo "#                  VALIDATION                  #"
echo "#                                              #"
echo "################################################"


echo ""
echo ""
echo ""
echo "====================================="
echo "--------------MaCFP_PMMA-------------"
echo "====================================="
run_case "$CAS_A" "run.sh"
echo "END OF: MaCFP_PMMA-------------------"
echo "====================================="















