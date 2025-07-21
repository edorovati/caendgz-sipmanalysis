#!/bin/bash
# rsleep 22:45
source ~/RICCARDO/myenv/bin/activate

# Lista dei valori di bias da usare
vbias_list="52 52.5 53 53.5 54 54.5 55 56 57 58"

# Parametri fissi
n_runs=2
sampling=2500

# Temporary changes
# n_runs=1

# Loop sui voltaggi
for vbias in ${vbias_list}; do
    echo ">>> Inizio run per VBIAS = $vbias V"
    ./run_analyze.sh --n_runs "$n_runs" --vbias "$vbias" --sampling "$sampling"
    echo ">>> Completato VBIAS = $vbias V"
    echo
done
