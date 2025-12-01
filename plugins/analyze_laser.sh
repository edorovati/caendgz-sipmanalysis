#!/bin/bash
set -e

show_help() {
  echo "Usage: $0 --filename PATH --sampling RATE --folder PATH [--vbias VALUE] [--channel N]"
  echo
  echo "Description:"
  echo "  Robust launcher for rooter-laser.py (handles absolute/relative paths and filename cleanup)."
  echo
  echo "Required arguments:"
  echo "  --filename   Path or base name of the acquisition (can include _chX.npz)."
  echo "  --sampling   Sampling rate in MS/s."
  echo "  --folder     Output folder path."
  echo
  exit 0
}

[[ $# -eq 0 ]] && show_help

# --- Parse args ---
while [[ $# -gt 0 ]]; do
  case $1 in
    --filename) FILENAME="$2"; shift 2 ;;
    --vbias) VBIAS="$2"; shift 2 ;;
    --sampling) SAMPLING="$2"; shift 2 ;;
    --folder) FOLDER="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;
    -h|--help) show_help ;;
    *) echo "Unknown option: $1"; show_help ;;
  esac
done

# --- Check required ---
if [[ -z "$FILENAME" || -z "$SAMPLING" || -z "$FOLDER" ]]; then
  echo "❌ Missing required arguments."
  echo
  show_help
fi

echo "🔦 LASER analysis on: $FILENAME (folder: $FOLDER) at $SAMPLING MS/s"
echo

# --- Normalize paths safely ---
FILENAME=$(realpath -m "$FILENAME")   # canonical path (even if partially invalid)
FOLDER=$(realpath -m "$FOLDER")

# --- Clean up possible doubled suffixes like _ch0.npz_ch1.npz ---
CLEAN_BASE=$(basename "$FILENAME")

# remove ALL trailing channel markers + .npz, even doubled
CLEAN_BASE=$(echo "$CLEAN_BASE" | sed -E 's@(_ch[0-9]+\.npz)+$@@')

# remove accidental duplicate extensions (e.g. .npz.npz)
CLEAN_BASE=${CLEAN_BASE%.npz}

# Rebuild absolute base path
BASE_DIR=$(dirname "$FILENAME")
BASE_PATH="${BASE_DIR}/${CLEAN_BASE}"

# --- Debug info ---
echo "🧩 Using base: $BASE_PATH"
echo "    → Laser: ${BASE_PATH}_ch0.npz"
echo "    → Ch1:   ${BASE_PATH}_ch1.npz"
echo "    → Output: ${FOLDER}/rooted_${CLEAN_BASE}.root"
echo

mkdir -p "$FOLDER"

python3 /home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/script/rooter-laser.py \
  --laser   "${BASE_PATH}_ch0.npz" \
  --ch1     "${BASE_PATH}_ch1.npz" \
  --soglia_laser 40 \
  --sampling "$SAMPLING" \
  --output  "${FOLDER}/data.vbias_{$VBIAS}.root"

echo -e "\e[33m>>> LASER analysis complete.\e[0m"
