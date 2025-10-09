#!/bin/bash
set -e

# === ARGUMENT PARSING ===
DO_TAU=false
DO_DARK=false
DO_LASER=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --filename) FILENAME="$2"; shift 2 ;;
    --vbias) VBIAS="$2"; shift 2 ;;
    --sampling) SAMPLING="$2"; shift 2 ;;
    --folder) FOLDER="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;
    --tau) DO_TAU=true; shift ;;
    --dark) DO_DARK=true; shift ;;
    --laser) DO_LASER=true; shift ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --filename <run-X> --vbias <voltage> --sampling <MHz> --folder <folder> --channel <channel> [--tau] [--dark] [--laser]"
      exit 1
      ;;
  esac
done

# === CHECK REQUIRED ARGS ===
for ARG_NAME in FILENAME VBIAS SAMPLING CHANNEL FOLDER; do
    if [[ -z "${!ARG_NAME}" ]]; then
        echo "❌ Missing required argument: --${ARG_NAME,,}"  # ,, trasforma in minuscolo
        exit 1
    fi
done

# === CHECK IF ANY ANALYSIS WAS SELECTED ===
if ! $DO_TAU && ! $DO_DARK && ! $DO_LASER; then
  echo "❌ No analysis selected. Please specify at least one of --tau, --dark, or --laser."
  exit 1
fi

echo -e "\e[34m🌟 Starting analysis for $FILENAME | ch$CHANNEL\e[0m"

# === RUN SELECTED ANALYSES IN PARALLELO ===
if $DO_TAU; then
  echo -e "\e[36m>>> TAU analysis\e[0m"
  ./../plugins/analyze_tau.sh \
    --filename "$FILENAME" \
    --vbias "$VBIAS" \
    --sampling "$SAMPLING" \
    --folder "$FOLDER" \
    --channel "$CHANNEL" &
  pid_tau=$!
fi

if $DO_DARK; then
  echo -e "\e[36m>>> DARK analysis\e[0m"
  ./../plugins/analyze_dark.sh \
    --filename "$FILENAME" \
    --vbias "$VBIAS" \
    --sampling "$SAMPLING" \
    --folder "$FOLDER" \
    --channel "$CHANNEL" &
  pid_dark=$!
fi

if $DO_LASER; then
  echo -e "\e[36m>>> LASER analysis\e[0m"
  ./../plugins/analyze_laser.sh \
    --filename "$FILENAME" \
    --vbias "$VBIAS" \
    --sampling "$SAMPLING" \
    --folder "$FOLDER" \
    --channel "$CHANNEL" &
  pid_laser=$!
fi

# === WAIT FOR ALL SELECTED ANALYSES ===
[[ $DO_TAU == true ]] && wait $pid_tau
[[ $DO_DARK == true ]] && wait $pid_dark
[[ $DO_LASER == true ]] && wait $pid_laser

echo -e "\e[32m✅ [${FILENAME} | ch${CHANNEL}] SELECTED ANALYSIS COMPLETE\e[0m"
