#!/bin/bash

# Default values (as per --help)
FILTER_ADC="0"
CHANNEL=(1)
MIN_EVENTS="10000"
LOG_FILE=""
SAMPLING="2500"
FOLDER=""
TRG="NIM"
SHIFT_WAVEFORMS=false
VBIAS=""
FILENAME=""

# Parse input arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --filter_ADC)
      FILTER_ADC=$2
      shift 2
      ;;
    --channel)
      shift
      CHANNEL=()  # Reset default
      while [[ $# -gt 0 && $1 != --* ]]; do
        CHANNEL+=("$1")
        shift
      done
      ;;
    --min_events)
      MIN_EVENTS=$2
      shift 2
      ;;
    --log_file)
      LOG_FILE=$2
      shift 2
      ;;
    --sampling)
      SAMPLING=$2
      shift 2
      ;;
    --vbias)
      VBIAS=$2
      shift 2
      ;;
    --folder)
      FOLDER=$2
      shift 2
      ;;
    --trg)
      TRG=$2
      shift 2
      ;;
    --shift_waveforms)
      SHIFT_WAVEFORMS=true
      shift
      ;;
    --filename)
      FILENAME=$2
      shift 2
      ;;
    -h|--help)
      python run_dgz.py --help
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Set default for filename if not provided
if [[ -z "$FILENAME" ]]; then
  FILENAME="TEMP"
  echo "[INFO] No --filename provided. Using default: \"$FILENAME\""
fi

# Assemble command
CMD="python run_dgz.py"
[[ -n "$FILTER_ADC" ]] && CMD+=" --filter_ADC $FILTER_ADC"
[[ ${#CHANNEL[@]} -gt 0 ]] && CMD+=" --channel ${CHANNEL[*]}"
[[ -n "$MIN_EVENTS" ]] && CMD+=" --min_events $MIN_EVENTS"
[[ -n "$LOG_FILE" ]] && CMD+=" --log_file $LOG_FILE"
[[ -n "$SAMPLING" ]] && CMD+=" --sampling $SAMPLING"
[[ -n "$VBIAS" ]] && CMD+=" --vbias $VBIAS" || { echo "Error: --vbias is required"; exit 1; }
[[ -n "$FOLDER" ]] && CMD+=" --folder $FOLDER"
[[ -n "$TRG" ]] && CMD+=" --trg $TRG"
$SHIFT_WAVEFORMS && CMD+=" --shift_waveforms"
[[ -n "$FILENAME" ]] && CMD+=" --filename $FILENAME"

# Print and execute
echo "Executing: $CMD"
eval "$CMD"





