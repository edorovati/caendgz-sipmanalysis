#!/bin/bash

# Default values
N_RUNS=""
VBIAS=""
SAMPLING=""
FOLDER=""
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
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --n_runs <number> --vbias <voltage> [--sampling <rate>] [--folder <foldername>]"
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$N_RUNS" || -z "$VBIAS" ]]; then
  echo "Error: --n_runs and --vbias are required. --sampling and --folder are optional."
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
  FILENAME_CH="run-${i}_ch1.npz"
  echo ">>> Starting acquisition $i with filename: $FILENAME"

  SECONDS=0  # Timer start

  python run_dgz.py \
    --vbias "$VBIAS" \
    --filter_ADC 0 \
    --filename "$FILENAME" \
    --sampling "$SAMPLING" \
    --trg laser \
    --channel 0 1 \
    --folder "$FOLDER" \
    --min_events "$N_WAVES"

  echo ">>> Acquisition $i completed. Output saved to $FILENAME_CH"

  sleep 5

  echo -e "\e[32m>>> [RUN $i] Launching median_wf.py\e[0m"
  python ../analysis/script/median_wf.py ../data/${FOLDER}/${FILENAME_CH} \
    --output ../data/${FOLDER}/vbias_${VBIAS}_run-${i}.npz \
    --amplitude_threshold 5.0 \
    --unit mV \
    --num_waveforms ${N_WAVES} \
    --sampling $SAMPLING \
    --save_npz &
  pid_median=$!

  wait $pid_median

  echo -e "\e[32m>>> [RUN $i] Analyzing tau from filtered waveforms (tau_direct)\e[0m"
  python ../analysis/script/tau_direct.py ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_filtered.npz \
    --output ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_tau_direct_all.txt \
    --sampling $SAMPLING

  echo -e "\e[32m>>> [RUN $i] Analyzing tau from median waveform (tau_direct)\e[0m"
  python ../analysis/script/tau_direct.py ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_median_wf.npz \
    --output ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_tau_direct_median.txt \
    --sampling $SAMPLING

  echo -e "\e[32m>>> [RUN $i] Extracting tau as fit constant from all filtered waveforms (tau_fit)\e[0m"
  python ../analysis/script/tau_fit.py ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_filtered.npz \
    --output ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_tau_fit_all.txt \
    --sampling $SAMPLING

  echo -e "\e[32m>>> [RUN $i] Extracting tau as fit constant from median waveform (tau_fit)\e[0m"
  python ../analysis/script/tau_fit.py ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_median_wf.npz \
    --output ../data/${FOLDER}/vbias_${VBIAS}_run-${i}_tau_fit_median.txt \
    --sampling $SAMPLING

  echo -e "\e[32m>>> [RUN $i] Launching get_transitions.py\e[0m"
  python ../analysis/script/get_transitions.py \
    --npz ../data/${FOLDER}/${FILENAME_CH} \
    --output ../data/${FOLDER}/vbias_${VBIAS}.run_${i}.get_transitions.txt \
    --scanthr \
    --range 0-50 \
    --sign 1 \
    --x_start 100 \
    --x_end 987 \
    --num_waveforms ${N_WAVES} &
  pid_get=$!

  wait $pid_get
  echo -e "\e[33m>>> [RUN $i] get_transitions.py finished\e[0m"

  duration=$SECONDS
  echo ">>> Time taken for run $i: $((duration / 60)) min $((duration % 60)) sec"

  # Delete waveform if needed
  if [[ "$DELETE_NPZ" == true && -f ../data/${FOLDER}/${FILENAME_CH} ]]; then
    rm ../data/${FOLDER}/${FILENAME_CH}
  fi

  # Rename files (optional)
  if [[ -f "../data/${FOLDER}/${FILENAME_CH}" ]]; then
    mv "../data/${FOLDER}/${FILENAME_CH}" "../data/${FOLDER}/vbias_${VBIAS}.${FILENAME_CH}"
  fi

  echo ">>> Analysis of run $i completed."

  # ─── Pretty progress bar ─────────────────────────
  BAR_WIDTH=30
  PROGRESS=$(( (i + 1) * BAR_WIDTH / N_RUNS ))
  BAR=$(printf "%-${BAR_WIDTH}s" "#" | cut -c1-$PROGRESS)
  SPACE=$(printf "%-${BAR_WIDTH}s" " " | cut -c1-$((BAR_WIDTH - PROGRESS)))
  printf "\r\e[36mProgress: [%s%s] %d/%d\e[0m\n" "$BAR" "$SPACE" $((i + 1)) $N_RUNS
done

# === MERGE FINAL RESULTS ===
echo ">>> Merging DCR and Tau data"

# Merge DCR
python ../analysis/script/code/merge.py \
  --merge ../data/${FOLDER}/vbias_${VBIAS}.*get_transitions.txt \
  --mode sum --columns 2 --sum 2 \
  --output ../data/${FOLDER}/merged.vbias_${VBIAS}.get_transitions.txt

# Merge tau: direct on median
python ../analysis/script/code/merge.py \
  --merge ../data/${FOLDER}/vbias_${VBIAS}_run-*_tau_direct_median.txt \
  --mode append --columns 3 \
  --output ../data/${FOLDER}/merged.vbias_${VBIAS}.tau_direct_median.txt

# Merge tau: direct on all
python ../analysis/script/code/merge.py \
  --merge ../data/${FOLDER}/vbias_${VBIAS}_run-*_tau_direct_all.txt \
  --mode append --columns 3 \
  --output ../data/${FOLDER}/merged.vbias_${VBIAS}.tau_direct_all.txt

# Merge tau: fit on median
python ../analysis/script/code/merge.py \
  --merge ../data/${FOLDER}/vbias_${VBIAS}_run-*_tau_fit_median.txt \
  --mode append --columns 6 \
  --output ../data/${FOLDER}/merged.vbias_${VBIAS}.tau_fit_median.txt

# Merge tau: fit on all
python ../analysis/script/code/merge.py \
  --merge ../data/${FOLDER}/vbias_${VBIAS}_run-*_tau_fit_all.txt \
  --mode append --columns 6 \
  --output ../data/${FOLDER}/merged.vbias_${VBIAS}.tau_fit_all.txt

# Turn output OFF
echo "Switching output OFF..."
echo "OP1 0" | nc -w1 -W1 "$host" "$port"
echo "Output is OFF. All runs completed."