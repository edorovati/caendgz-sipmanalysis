#!/bin/bash

set -e

# Parse input arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --filename)
      FILENAME="$2"
      shift 2
      ;;
    --vbias)
      VBIAS="$2"
      shift 2
      ;;
    --sampling)
      SAMPLING="$2"
      shift 2
      ;;
    --folder)
      FOLDER="$2"
      shift 2
      ;;
    --channel)
      CHANNEL="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --filename <run-X> --vbias <voltage> --sampling <MHz> --folder <folder> --channel <channel>"
      exit 1
      ;;
  esac
done

# Check mandatory arguments
if [[ -z "$FILENAME" || -z "$VBIAS" || -z "$SAMPLING" || -z "$FOLDER" || -z "$CHANNEL" ]]; then
  echo "Missing required arguments."
  exit 1
fi

# Construct file path
INPUT_NPZ="../data/${FOLDER}/${FILENAME}_ch${CHANNEL}.npz"

# Output name base
OUTBASE="../data/${FOLDER}/vbias_${VBIAS}_${FILENAME}_ch${CHANNEL}"

echo "📈 Analyzing file: $INPUT_NPZ (ch $CHANNEL)"

# === Median waveform ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] median_wf.py\e[0m"
python ../analysis/script/median_wf.py "$INPUT_NPZ" \
  --output "${OUTBASE}.npz" \
  --amplitude_threshold 5.0 \
  --unit mV \
  --sampling "$SAMPLING" \
  --save_npz &
pid_median=$!

wait $pid_median

# === Tau from filtered waveforms (direct) ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] tau_direct from filtered\e[0m"
python ../analysis/script/tau_direct.py "${OUTBASE}_filtered.npz" \
  --output "${OUTBASE}_tau_direct_all.txt" \
  --sampling "$SAMPLING"

# === Tau from median waveform (direct) ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] tau_direct from median\e[0m"
python ../analysis/script/tau_direct.py "${OUTBASE}_median_wf.npz" \
  --output "${OUTBASE}_tau_direct_median.txt" \
  --sampling "$SAMPLING"

# === Tau from filtered (fit) ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] tau_fit from filtered\e[0m"
python ../analysis/script/tau_fit.py "${OUTBASE}_filtered.npz" \
  --output "${OUTBASE}_tau_fit_all.txt" \
  --sampling "$SAMPLING"

# === Tau from median (fit) ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] tau_fit from median\e[0m"
python ../analysis/script/tau_fit.py "${OUTBASE}_median_wf.npz" \
  --output "${OUTBASE}_tau_fit_median.txt" \
  --sampling "$SAMPLING"

# === Get transitions ===
echo -e "\e[32m>>> [${FILENAME} | ch${CHANNEL}] get_transitions.py\e[0m"
python ../analysis/script/get_transitions.py \
  --npz "$INPUT_NPZ" \
  --output "${OUTBASE}.get_transitions.txt" \
  --scanthr \
  --range 0-50 \
  --sign 1 \
  --x_start 100 \
  --x_end 987 &
pid_get=$!

wait $pid_get

echo -e "\e[33m>>> [${FILENAME} | ch${CHANNEL}] analysis complete.\e[0m"