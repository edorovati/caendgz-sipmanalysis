#!/bin/bash
source ~/RICCARDO/myenv/bin/activate

if [ $# -ne 1 ]; then
    echo "Usage: $0 <cartella_base>"
    exit 1
fi

base_path="$1"

vbias_list="56"
# vbias_list="52.5 53 53.5 54 54.5 55 56 57 58"

for vbias in $vbias_list; do
    vbias_label=${vbias/./_}
    input_file="${base_path}/vbias_${vbias}.run-0.npz"
    output_file="${base_path}/tau_${vbias_label}-run0.txt"

    echo "Processing Vbias=${vbias}..."
    python tau_median.py "$input_file" \
        --num_waveforms 10000 \
        --amplitude_threshold 4.0 \
        --unit mV \
        --sampling 2500 \
        --output "$output_file"\
        --save_npz \

done
