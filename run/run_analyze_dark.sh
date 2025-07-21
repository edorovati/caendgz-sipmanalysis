#!/bin/bash

# Default values
N_RUNS=""
VBIAS=""
SAMPLING=""
FOLDER=""
CHANNELS=()  # canali da analizzare
host="aimtti-plh120p-00"
port=9221

N_WAVES=10000
DELETE_NPZ=true

# Parse input arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --n_runs)
      N_RUNS="$2"
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
    --channels)
      shift
      echo "📥 Parsing channels..."
      while [[ $# -gt 0 && $1 != --* ]]; do
        echo "  ➕ Aggiungo canale: $1"
        CHANNELS+=("$1")
        shift
      done
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --n_runs <number> --vbias <voltage> --sampling <rate> --folder <foldername> [--channels 1 2 3 ...]"
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$N_RUNS" || -z "$VBIAS" || -z "$SAMPLING" || -z "$FOLDER" ]]; then
  echo "Error: --n_runs, --vbias, --sampling and --folder are required."
  exit 1
fi

# Powering on the bias
echo "Setting V1 to $VBIAS V..."
echo "V1 $VBIAS" | nc -w1 -W1 "$host" "$port"
echo "Switching output ON..."
echo "OP1 1" | nc -w1 -W1 "$host" "$port"
sleep 10
echo "Output is ON and V1 is set to $VBIAS V."

# === RUN LOOP ===
for ((i=0; i<N_RUNS; i++)); do
  FILENAME="run-$i"
  echo ">>> Starting acquisition $i with filename: $FILENAME"

  SECONDS=0

  echo -e "\n🐞 DEBUG: Chiamo run_dgz.py con:"
  echo -n "python run_dgz.py "
  echo -n "--vbias $VBIAS "
  echo -n "--filter_ADC 0 "
  echo -n "--filename $FILENAME "
  echo -n "--sampling $SAMPLING "
  echo -n "--trg laser "
  echo -n "--channel 0 ${CHANNELS[*]} "
  echo -n "--folder $FOLDER "
  echo -n "--min_events $N_WAVES"
  echo -e "\n"

  python run_dgz.py \
    --vbias "$VBIAS" \
    --filter_ADC 0 \
    --filename "$FILENAME" \
    --sampling "$SAMPLING" \
    --trg laser \
    --channel 0 "${CHANNELS[@]}" \
    --folder "$FOLDER" \
    --min_events "$N_WAVES"

  echo ">>> Acquisition $i completed."

  sleep 5

  # # === LOOP SU CANALI (tranne 0) ===
  # for CH in "${CHANNELS[@]}"; do
  #   echo -e "\e[32m>>> [RUN $i | ch$CH] Analisi in corso...\e[0m"
  #   ./analyze_channel.sh \
  #     --filename "$FILENAME" \
  #     --vbias "$VBIAS" \
  #     --sampling "$SAMPLING" \
  #     --folder "$FOLDER" \
  #     --channel "$CH" &
  #   pids+=($!)
  # done

  # # Aspetta tutte le analisi del run
  # for pid in "${pids[@]}"; do
  #   wait $pid
  # done
  # unset pids

  duration=$SECONDS
  echo ">>> Time taken for run $i: $((duration / 60)) min $((duration % 60)) sec"
  echo ">>> Run $i completed."

done

# === Power OFF ===
echo "Switching output OFF..."
echo "OP1 0" | nc -w1 -W1 "$host" "$port"
echo "✅ Output is OFF. All runs completed."