#!/bin/bash
set -e

# === ARGUMENT PARSING ===
while [[ $# -gt 0 ]]; do
  case $1 in
    --filename) FILENAME="$2"; shift 2 ;;
    --vbias) VBIAS="$2"; shift 2 ;;  # non usato ma tenuto per consistenza
    --sampling) SAMPLING="$2"; shift 2 ;;
    --folder) FOLDER="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;  # non usato ma mantenuto
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# === CHECK ARGS ===
if [[ -z "$FILENAME" || -z "$SAMPLING" || -z "$FOLDER" ]]; then
  echo "Missing required arguments."; exit 1
fi

echo "🔦 LASER analysis on: $FILENAME (folder: $FOLDER)"

python3 /home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/script/rooter-laser.py \
  --laser "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/data/${FOLDER}/${FILENAME}_ch0.npz" \
  --soglia_laser 0 \
  --ch1   "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/data/${FOLDER}/${FILENAME}_ch1.npz" \
  --sampling "$SAMPLING" \
  --output "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/data/${FOLDER}/rooted_${FILENAME}.root"

echo -e "\e[33m>>> LASER analysis complete.\e[0m"