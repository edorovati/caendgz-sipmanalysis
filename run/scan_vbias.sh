#!/usr/bin/env bash

# Lista dei valori di vbias
vbias_list=(55)
#vbias_list=(33.5 34 34.5 35 35.5 36 36.5 37 37.5 38 39 40 41)
#vbias_list=(32.5 33 33.5 34 34.5 35 35.5) 
#vbias_list=(55)
# Percorso dello script vbias
VBIAS_SCRIPT="../plugins/vbias.sh"
LASER_SCRIPT="../plugins/laser.sh"

echo "-----------------------------------------"
echo "Accendo laser"
"$LASER_SCRIPT" --off
for bias in "${vbias_list[@]}"; do
    echo "-----------------------------------------"
    echo "Imposto vbias a $bias"
    "$VBIAS_SCRIPT" --vbias "$bias" --on

    echo "Attendo 10 secondi per stabilizzazione..."
    sleep 10

    echo "Eseguo run_analyze con vbias $bias"
    ./run_analyze.sh \
        --n_runs 5 \
        --vbias "$bias" \
        --sampling 750 \
        --folder "../data/test/laser_test/alcor_testing_dark" \
        --dark \
        --channels 1
done

echo "-----------------------------------------"
echo "Spengo vbias"
"$VBIAS_SCRIPT" --vbias 0 --off
echo "-----------------------------------------"
echo "Spegno laser"
"$LASER_SCRIPT" --off
