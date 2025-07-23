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

echo "📈 TAU analysis on: $INPUT_NPZ (ch $CHANNEL)"

# === Median waveform ===
echo -e "\e[32m>>> median_wf.py\e[0m"
python ../analysis/script/median_wf.py "$INPUT_NPZ" \
  --output "${OUTBASE}.npz" \
  --amplitude_threshold 5.0 \
  --unit mV \
  --sampling "$SAMPLING" \
  --save_npz &
pid_median=$!

wait $pid_median

# === Tau analysis ===
python ../analysis/script/tau_direct.py "${OUTBASE}_filtered.npz" \
  --output "${OUTBASE}_tau_direct_all.txt" \
  --sampling "$SAMPLING"

python ../analysis/script/tau_direct.py "${OUTBASE}_median_wf.npz" \
  --output "${OUTBASE}_tau_direct_median.txt" \
  --sampling "$SAMPLING"

python ../analysis/script/tau_fit.py "${OUTBASE}_filtered.npz" \
  --output "${OUTBASE}_tau_fit_all.txt" \
  --sampling "$SAMPLING"

python ../analysis/script/tau_fit.py "${OUTBASE}_median_wf.npz" \
  --output "${OUTBASE}_tau_fit_median.txt" \
  --sampling "$SAMPLING"

echo -e "\e[33m>>> TAU analysis complete.\e[0m"