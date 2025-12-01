#!/usr/bin/env bash

# Lista dei valori di vbias
vbias_list=(39)
#vbias_list=(35.5) 
# Percorso dello script vbias
VBIAS_SCRIPT="../plugins/vbias.sh"
LASER_SCRIPT="../plugins/laser.sh"

echo "-----------------------------------------"
echo "Accendo laser"
"$LASER_SCRIPT" --on
for bias in "${vbias_list[@]}"; do
    echo "-----------------------------------------"
    echo "Imposto vbias a $bias"
    "$VBIAS_SCRIPT" --vbias "$bias" --on

    echo "Attendo 10 secondi per stabilizzazione..."
    sleep 10

    echo "Eseguo run_analyze con vbias $bias"
    ./run_analyze_copy_copy_copy.sh \
        --n_runs 1 \
        --vbias "$bias" \
        --sampling 5000 \
        --folder "../data/test/laser_BSI/2025_12_01/BSI_sn2_A1/GB" \
        --laser \
        --channels 1
done

echo "-----------------------------------------"
echo "Spengo vbias"
"$VBIAS_SCRIPT" --vbias 0 --off
echo "-----------------------------------------"
echo "Spegno laser"
"$LASER_SCRIPT" --off
