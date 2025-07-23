#!/bin/bash
set -e

# === ARGUMENT PARSING ===
while [[ $# -gt 0 ]]; do
  case $1 in
    --filename) FILENAME="$2"; shift 2 ;;
    --vbias) VBIAS="$2"; shift 2 ;;
    --sampling) SAMPLING="$2"; shift 2 ;;
    --folder) FOLDER="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# === CHECK ARGS ===
if [[ -z "$FILENAME" || -z "$VBIAS" || -z "$SAMPLING" || -z "$FOLDER" || -z "$CHANNEL" ]]; then
  echo "Missing required arguments."; exit 1
fi

INPUT_NPZ="../data/${FOLDER}/${FILENAME}_ch${CHANNEL}.npz"
OUTBASE="../data/${FOLDER}/vbias_${VBIAS}_${FILENAME}_ch${CHANNEL}"

echo "🌑 DARK analysis on: $INPUT_NPZ (ch $CHANNEL)"

python ../analysis/script/get_transitions.py \
  --npz "$INPUT_NPZ" \
  --output "${OUTBASE}.get_transitions.txt" \
  --scanthr \
  --range 0-20 \
  --sign 1 \
  --voltage_step_mV 0.25 \
  --x_start 100 \
  --x_end 987

echo -e "\e[33m>>> DARK analysis complete.\e[0m"