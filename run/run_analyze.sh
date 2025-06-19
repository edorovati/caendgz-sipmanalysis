#!/bin/bash

# Default values
N_RUNS=""
VBIAS=""
SAMPLING=""
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
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --n_runs <number> --vbias <voltage>"
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$N_RUNS" || -z "$VBIAS" ]]; then
  echo "Error: --n_runs and --vbias are required. --sampling is optional."
  echo "Usage: $0 --n_runs <number> --vbias <voltage>"
  exit 1
fi

# powering on the bias
echo "Setting V1 to $VBIAS V..."
echo "V1 $VBIAS" | nc -w1 -W1 "$host" "$port"
# Turn output ON
echo "Switching output ON..."
echo "OP1 1" | nc -w1 -W1 "$host" "$port"
sleep 10  # let output stabilize and give message
echo "Output is ON and V1 is set to $VBIAS V."

# Run loop
for ((i=0; i<N_RUNS; i++)); do
  FILENAME="run-$i.npz"
  echo ">>> Starting acquisition $i with filename: $FILENAME"

  SECONDS=0  # reset timer

  ./run_dgz.sh \
    --vbias "$VBIAS" \
    --filename "$FILENAME" \
    --sampling "$SAMPLING"
  
  #analysis
  echo ">>> Acquisition $i completed. Output saved to $FILENAME. Analyzing data..."
  echo "Measuring temperature"
  /home/eic/TIFPA/Keithley/taking_temperature.sh --namefile /home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/data/$VBIAS.temperature

  sleep 5
  echo -e "\e[32m>>> [RUN $i] Launching median_wf.py\e[0m"
  python ../analysis/script/median_wf.py ../data/$FILENAME --output ../data/vbias_${VBIAS}_run-${i} --amplitude_threshold 4.0 --unit mV --num_waveforms ${N_WAVES} --sampling $SAMPLING --save_npz 2>&1 &
  pid_median=$!
  
  # # #tau_calculations!
  
  wait $pid_median
  echo -e "\e[32m>>> [RUN $i] Analyzing tau as pure constant from all the filtered waveforms\e[0m"x
  python ../analysis/script/tau_direct.py ../data/vbias_${VBIAS}_run-${i}_filtered.npz --output ../data/vbias_${VBIAS}_run-${i}_tau_direct_all.txt --sampling $SAMPLING 2>&1 &
  wait $pid_median
  echo -e "\e[32m>>> [RUN $i] Analyzing tau as pure constant from the median waveform\e[0m"x
  python ../analysis/script/tau_direct.py ../data/vbias_${VBIAS}_run-${i}_median_wf.npz --output ../data/vbias_${VBIAS}_run-${i}_tau_direct_median.txt --sampling $SAMPLING 2>&1 &

  echo -e "\e[32m>>> [RUN $i] Analyzing tau with fit\e[0m"
  # fit all filtered wf
  wait $pid_median
  echo -e "\e[32m>>> [RUN $i] Extracting tau as fit constant from all the filtered waveforms\e[0m"
  python ../analysis/script/tau_fit.py ../data/vbias_${VBIAS}_run-${i}_filtered.npz --output ../data/vbias_${VBIAS}_run-${i}_tau_fit_all.txt --sampling $SAMPLING 2>&1 &
  # fit the median wf
  wait $pid_median
  echo -e "\e[32m>>> [RUN $i] Extracting tau as fit constant from the median waveform\e[0m"x
  python ../analysis/script/tau_fit.py ../data/vbias_${VBIAS}_run-${i}_median_wf.npz --output ../data/vbias_${VBIAS}_run-${i}_tau_fit_median.txt --sampling $SAMPLING 2>&1 &
  # pid_tau_fit=$!
  

  # GET TRANSITIONS
  echo -e "\e[32m>>> [RUN $i] Launching get_transitions.py\e[0m"
  python ../analysis/script/get_transitions.py --npz ../data/$FILENAME --output ../data/vbias_${VBIAS}.run_${i}.get_transitions.txt --scanthr --range 0-50 --sign 1 --x_start 100 --x_end 987 --num_waveforms ${N_WAVES} 2>&1 &
  pid_get=$!
  






  # wait $pid_get
  # echo -e "\e[33m>>> [RUN $i] get_transitions.py finished\e[0m"
  
  wait $pid_get
  echo -e "\e[33m>>> [RUN $i] tau.py finished\e[0m"
  
  echo ">>> Analysis of $i-run completed."

  duration=$SECONDS
  echo ">>> Time taken for run $i: $((duration / 60)) min $((duration % 60)) sec"

  # delete waveforms if requested
  [[ "$DELETE_NPZ" == true ]] && rm ../data/$FILENAME

  ### this is a bit messy, I am renaming the NPZ to avoid overwriting for different bias voltages
  mv ../data/${FILENAME} ../data/vbias_${VBIAS}.${FILENAME}
  
done

# merge dcr data and tau data
# echo ">>> Merging DCR and Tau data"
echo ">>> Merging DCR and Tau data"

# Merge DCR
python ../analysis/script/code/merge.py \
  --merge ../data/vbias_${VBIAS}.*get_transitions.txt \
  --mode sum --columns 2 --sum 2 \
  --output ../data/merged.vbias_${VBIAS}.get_transitions.txt

# Merge tau: direct on median
python ../analysis/script/code/merge.py \
  --merge ../data/vbias_${VBIAS}_run-*_tau_direct_median.txt \
  --mode append --columns 3 \
  --output ../data/merged.vbias_${VBIAS}.tau_direct_median.txt

# Merge tau: direct on all
python ../analysis/script/code/merge.py \
  --merge ../data/vbias_${VBIAS}_run-*_tau_direct_all.txt \
  --mode append --columns 3 \
  --output ../data/merged.vbias_${VBIAS}.tau_direct_all.txt

# Merge tau: fit on median
python ../analysis/script/code/merge.py \
  --merge ../data/vbias_${VBIAS}_run-*_tau_fit_median.txt \
  --mode append --columns 6 \
  --output ../data/merged.vbias_${VBIAS}.tau_fit_median.txt

# Merge tau: fit on all
python ../analysis/script/code/merge.py \
  --merge ../data/vbias_${VBIAS}_run-*_tau_fit_all.txt \
  --mode append --columns 6 \
  --output ../data/merged.vbias_${VBIAS}.tau_fit_all.txt


# Turn output OFF after all runs
echo "Switching output OFF..."
echo "OP1 0" | nc -w1 -W1 "$host" "$port"
echo "Output is OFF. All runs completed."
