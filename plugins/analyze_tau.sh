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


# === CHECK REQUIRED ARGS IN ONE LOOP ===
for ARG_NAME in FILENAME VBIAS SAMPLING CHANNEL FOLDER; do
    if [[ -z "${!ARG_NAME}" ]]; then
        echo "❌ Missing required argument: --${ARG_NAME,,}"  # ,, trasforma in minuscolo
        exit 1
    fi
done


# === NORMALIZE FILENAME ===
BASENAME=$(basename "$FILENAME" .npz)

# --- Input path handling ---
if [[ "$BASENAME" == *_ch* ]]; then
    INPUT_NPZ="$FILENAME"
    [[ "$INPUT_NPZ" != *.npz ]] && INPUT_NPZ="${INPUT_NPZ}.npz"
else
    if [[ "$FILENAME" == *"/"* ]]; then
        INPUT_NPZ="${FILENAME}_ch${CHANNEL}.npz"
    else
        if [[ -n "$FOLDER" ]]; then
            INPUT_NPZ="${FOLDER%/}/${BASENAME}_ch${CHANNEL}.npz"
        else
            INPUT_NPZ="./${BASENAME}_ch${CHANNEL}.npz"
        fi
    fi
fi

# --- Output directory ---
if [[ -n "$FOLDER" ]]; then
    OUTDIR="${FOLDER%/}"  # remove trailing slash
else
    OUTDIR=$(dirname "$INPUT_NPZ")
fi
mkdir -p "$OUTDIR"

# --- OUTBASE handling ---
if [[ "$BASENAME" == *_ch* ]]; then
    OUTBASE="${OUTDIR}/vbias_${VBIAS}_${BASENAME}"
else
    OUTBASE="${OUTDIR}/vbias_${VBIAS}_${BASENAME}_ch${CHANNEL}"
fi

echo "🌟 Starting TAU analysis"
echo "📈 TAU analysis on: $INPUT_NPZ (ch $CHANNEL)"
echo "Output base: $OUTBASE"

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
