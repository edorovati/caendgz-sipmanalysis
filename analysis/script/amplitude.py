import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

# === Local imports ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from utils import Utils
from waveform_analysis import waveform_analysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Persistence plot from filtered waveform set.")
    parser.add_argument("path", help="Path to .npz file")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--num_waveforms", type=int, default=None)
    parser.add_argument("--amplitude_threshold", type=float, default=4.0)
    parser.add_argument("--sampling", type=float, required=True, help="Sampling rate in MHz")
    parser.add_argument("--align", default="False", help="Alignment method: False, 'max', or 'auto'")
    parser.add_argument("--output", type=str, required=True, help="Output filename prefix (no extension)")
    parser.add_argument("--save_npz", action="store_true", help="Save filtered waveforms to output_name.npz")
    parser.add_argument("--save_txt", action="store_true", help="Save peak values to output_name_peaks.txt")

    args = parser.parse_args()

    npz_path = args.path
    sampling = args.sampling
    unit = args.unit
    align = args.align
    num_waveforms = args.num_waveforms
    amplitude_threshold = args.amplitude_threshold

    wf_data = Utils.load_waveforms(npz_path)
    info = Utils.get_info(npz_path)
    waveforms = wf_data.get("waveforms", None)
    
    # # === Analyze waveforms ===
    baseline_range = (49, 973)
    waveforms, rising_edges = waveform_analysis.analyze_waveforms(
        wf_data=wf_data,
        fs=sampling,
        baseline_range=baseline_range,
        amplitude_threshold=3.8,
        max_waveforms=num_waveforms,
        peak_index_range=(300, 600),
        unit=unit,
        get_edges=True
    )

    if not waveforms:
        print("⚠️ No waveforms passed the filter.")
        exit()

    # # === Alignment ===
    # if align and align != "False":
    #     waveforms = waveform_analysis.align_waveforms(waveforms, rising_edges, mode=align)

    # if not waveforms:
    #     print("⚠️ No valid waveforms to plot.")
    #     exit()

    
    if waveforms is None or len(waveforms) == 0:
        print("❌ No waveforms found in the file.")
        exit(1)

    # === Compute amplitudes ===
    output_txt = args.output + "_peaks.txt" if args.output else None
    peak_vals = waveform_analysis.amplitudes(
        waveforms=waveforms,
        sampling=args.sampling,
        output_filename=output_txt if args.save_txt else None
    )

    # === Plot histogram ===
    plt.figure(figsize=(8, 5))
    plt.hist(peak_vals, bins=1000, color="steelblue", edgecolor="black", alpha=0.8)
    plt.title("Distribuzione delle ampiezze di picco")
    plt.xlabel("Ampiezza [mV]")
    #set range
    #set y log scale...
    # plt.yscale("log")
    plt.xlim(left = 0, right = 20)
    plt.ylabel("Conteggio")
    plt.grid(True)

    # if args.save_txt and args.output:
    #     out_plot = args.output + "_hist.png"
    #     plt.savefig(out_plot)
    #     print(f"📊 Istogramma salvato in: {out_plot}")

    plt.show()