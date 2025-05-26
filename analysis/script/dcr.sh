#!/usr/bin/env bash
set -euo pipefail    # Strict mode to catch errors early

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION ✱ Set paths and directories
#   • DATA_DIR      ⇒ Directory where data is stored
#   • SCRIPT_PATH   ⇒ Path to the script for getting transitions
#   • SCAN_DIR      ⇒ Directory for scanning files
# ──────────────────────────────────────────────────────────────────────────────
DATA_DIR="Data"
SCRIPT_PATH="code/get_transitions.py"
SCAN_DIR="Scan"

# ──────────────────────────────────────────────────────────────────────────────
# READ COMMAND-LINE ARGUMENTS ✱ Accept input from user
#   • ans1 ⇒ First answer (yes/no) activate threshold scan
#   • ans2 ⇒ Second answer (yes/no) activate produce root with data
#   • ans3 ⇒ Third answer (yes/no) activate plotting
# ──────────────────────────────────────────────────────────────────────────────
ans1=${1:-}
ans2=${2:-}
ans3=${3:-}

# ──────────────────────────────────────────────────────────────────────────────
# VALIDATE INPUT ✱ Ensure required arguments are passed and valid
#   • Shows usage if any argument is missing
#   • Disallows the forbidden combination: yes no yes
# ──────────────────────────────────────────────────────────────────────────────
if [[ -z "$ans1" || -z "$ans2" || -z "$ans3" ]]; then
  echo -e "\033[1;31m✘ Error: Missing arguments.\033[0m"
  echo -e "Usage: $0 yes|no yes|no yes|no"
  echo -e "Example: $0 yes yes no"
  exit 1
fi

if [[ "$ans1" == "yes" && "$ans2" == "no" && "$ans3" == "yes" ]]; then
  echo -e "\033[1;31m✘ Error: Combination 'yes no yes' is not allowed.\033[0m"
  exit 1
fi

# ──────────────────────────────────────────────────────────────────────────────
# CONDA ENVIRONMENT SETUP ✱ Check for Conda and activate environment
#   • Initializes conda shell if available
#   • Activates 'root-env' where needed
# ──────────────────────────────────────────────────────────────────────────────
if command -v conda &> /dev/null; then
  echo -e "\033[1;34mℹ Conda detected. Initializing...\033[0m"
  set +u
  eval "$(conda shell.bash hook)"
  set -u
  CONDA_AVAILABLE=true
else
  echo -e "\033[1;33m⚠ Conda not found: skipping environment activation.\033[0m"
  CONDA_AVAILABLE=false
fi

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 1 ✱ Initial scan (placeholder)
#   • Scanning file: python3 "$SCRIPT_PATH" --npz "$file" --scanthr --range 0-1 --sign 1 --x_start 100 --x_end 900
#   • You can edit --range 0-1 which is the range in mV of scan in thr
#   • You can edit --x_start 100 --x_end 900 and use less or more in terms of range
#   • You can add --lowpass 200e6 to filter signal and help in analysis
#   • Please refer to read me and comments on code/get_transitions.py to found other functionalities
#   • Moving txt produced in folder
# ──────────────────────────────────────────────────────────────────────────────
echo -e "\033[1;36m▶ PHASE 1: Running initial scan... (placeholder)\033[0m"
if [[ "$ans1" == "yes" ]]; then
    mkdir -p "$SCAN_DIR"
    for file in "$DATA_DIR"/*.npz; do
        echo -e "\033[1;36mAnalyze: $file\033[0m"
        python3 "$SCRIPT_PATH" --npz "$file" --scanthr --range 0-20 --sign 1 --x_start 100 --x_end 987 --num_waveforms 10000
    done
    mv -- *.txt "$SCAN_DIR"/
fi



# ──────────────────────────────────────────────────────────────────────────────
# PHASE 2 ✱ Execute dcr.cpp if ans2 is 'yes'
#   • Conditionally runs ROOT macro with optional conda env
#   • Produces a root file with TDirectory which are the overvoltage avaialable
#   • Produces a TTree which contains main information from analysis of peaks
# ──────────────────────────────────────────────────────────────────────────────
if [[ "${ans2}" == "yes" ]]; then
  echo -e "\033[1;32m▶ PHASE 2: Executing dcr.cpp...\033[0m"
  if [[ "$CONDA_AVAILABLE" == true ]]; then
    set +u
    conda activate root-env
    set -u
  fi
  root -l -q code/dcr.cpp
fi

# ──────────────────────────────────────────────────────────────────────────────
# PHASE 3 ✱ Execute dcr_plot.cpp if ans3 is 'yes'
#   • Conditionally runs ROOT plotting macro
#   • Activates Conda env if it wasn’t activated in Phase 2
#   • Produces main plots from the TTree info, obtained by previous phase
#   • Not start it unless you have a output.root file please
# ──────────────────────────────────────────────────────────────────────────────
if [[ "${ans3}" == "yes" ]]; then
  echo -e "\033[1;32m▶ PHASE 3: Executing dcr_plot.cpp...\033[0m"
  if [[ "$CONDA_AVAILABLE" == true && "${ans2}" != "yes" ]]; then
    set +u
    conda activate root-env
    set -u
  fi
  root -l code/dcr_plot.cpp
fi

# ──────────────────────────────────────────────────────────────────────────────
# FINAL MESSAGE ✱ Wrap up and exit
# ──────────────────────────────────────────────────────────────────────────────
echo -e "\033[1;34m✓ All tasks completed successfully.\033[0m"
