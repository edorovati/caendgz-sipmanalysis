#!/bin/bash
set -e

show_help() {
  echo ""
  echo "Usage: ./analyze_dark.sh --filename PATH --vbias VALUE --sampling RATE --folder DIR --channel N"
  echo ""
  echo "Required arguments:"
  echo "  --filename   Full path to the .npz file to analyze"
  echo "  --vbias      Bias voltage in volts (e.g. 53.0)"
  echo "  --sampling   Sampling frequency (e.g. 5000)"
  echo "  --folder     Folder name for output organization (same as original script)"
  echo "  --channel    Channel number to analyze (e.g. 0,1,...)"
  echo ""
  echo "Example:"
  echo "  ./analyze_dark.sh --filename ../data/test/run-1_ch1.npz --vbias 53.0 --sampling 5000 --folder test --channel 1"
  echo ""
  echo "This script runs the dark-count analysis via dark_transition.py."
  echo "Output files are saved as:"
  echo "  <folder>/vbias_<vbias>_<filename>_ch<channel>.dark_transition.txt"
  echo ""
}

# === HELP OPTION ===
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  show_help
  exit 0
fi

# === ARGUMENT PARSING ===
while [[ $# -gt 0 ]]; do
  case $1 in
    --filename) FILENAME="$2"; shift 2 ;;
    --vbias) VBIAS="$2"; shift 2 ;;
    --sampling) SAMPLING="$2"; shift 2 ;;
    --folder) FOLDER="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; show_help; exit 1 ;;
  esac
done

# === CHECK ARGS ===
if [[ -z "$FILENAME" || -z "$VBIAS" || -z "$SAMPLING" || -z "$FOLDER" || -z "$CHANNEL" ]]; then
  echo "Missing required arguments."
  show_help
  exit 1
fi

# === MAIN LOGIC ===
INPUT_NPZ="$FILENAME"  # usiamo il percorso completo senza modifiche

# Rimuoviamo .npz dal basename
BASENAME="$(basename "${INPUT_NPZ%.npz}")"

# Evitiamo duplicazioni di _ch se già presente
if [[ "$BASENAME" == *_ch* ]]; then
  OUTBASE="../data/${FOLDER}/vbias_${VBIAS}_${BASENAME}"
else
  OUTBASE="../data/${FOLDER}/vbias_${VBIAS}_${BASENAME}_ch${CHANNEL}"
fi

# Normalizziamo il percorso per evitare doppi slash e ../
OUTBASE="$(realpath -m "$OUTBASE")"

echo "🌑 DARK analysis on: $INPUT_NPZ (ch $CHANNEL)"

python ../analysis/script/dark_transition.py \
  --npz "$INPUT_NPZ" \
  --output "${OUTBASE}.dark_transition.txt" \
  --scanthr \
  --range 0-200 \
  --sign 1 \
  --lowpass 200e6 \
  --step_voltage_mV 0.5 \
  --x_start 100 \
  --x_end 987

echo -e "\e[33m>>> DARK analysis complete.\e[0m"
