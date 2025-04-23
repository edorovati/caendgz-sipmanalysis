#!/bin/bash

# === CONFIG ===
DATA_DIR="Data"
SCRIPT_PATH="code/get_transitions.py"
SCAN_DIR="Scan"

# === LEGGI ARGOMENTI DA LINEA DI COMANDO ===
ans1=$1
ans2=$2
ans3=$3

# === VERIFICA INPUT ===
if [[ -z "$ans1" || -z "$ans2" || -z "$ans3" ]]; then
    echo "Devi passare esattamente 3 argomenti (yes/no)."
    echo "Esempio: ./dcr.sh yes yes yes"
    exit 1
fi

# === CHECK COMBINAZIONE VIETATA ===
if [[ "$ans1" == "yes" && "$ans2" == "no" && "$ans3" == "yes" ]]; then
    echo "Combinazione sbagliata: yes no yes"
    exit 1
fi

# === ESECUZIONE PYTHON ===
if [[ "$ans1" == "yes" ]]; then
    echo "Inizio scan sugli .npz..."
    mkdir -p "$SCAN_DIR"
    for file in "$DATA_DIR"/*.npz; do
        echo "Analizzo: $file"
        python3 "$SCRIPT_PATH" --npz "$file" --scanthr --range 0-1 --sign 1 --x_start 100 --x_end 900
    done
    echo "Sposto i file .txt in '$SCAN_DIR'..."
    mv -- *.txt "$SCAN_DIR"/
fi

# === ROOT: dcr.cpp ===
if [[ "$ans2" == "yes" ]]; then
    echo "Attivo ambiente conda e lancio dcr.cpp..."
    conda activate root-env
    root -l -q code/dcr.cpp
fi

# === ROOT: dcr_plot.cpp ===
if [[ "$ans3" == "yes" ]]; then
    echo "Lancio dcr_plot.cpp..."
    # Attiva l'ambiente solo se non già fatto sopra
    if [[ "$ans2" != "yes" ]]; then
        conda activate root-env
    fi
    root -l code/dcr_plot.cpp
fi

echo "Script completato."
