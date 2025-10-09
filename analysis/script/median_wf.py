import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

# === Local imports ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from utils import Utils
from waveform_analysis import waveform_analysis
from plotting import Plotting

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute tau and overlay waveforms from an .npz file.")
    parser.add_argument("npz_path", help="Path to the .npz file.")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--amplitude_threshold", type=float, default=5.0)
    parser.add_argument("--num_waveforms", type=int, default=10000)
    parser.add_argument("--color", default="green")
    parser.add_argument("--align", default="False")
    parser.add_argument("--sampling", type=float, required=True)
    parser.add_argument("--output", type=str, required=True, help="Base name for output files (no extension)")
    parser.add_argument("--save_npz", action="store_true", help="Save filtered and median waveforms to .npz files")
    parser.add_argument("--save_plot", action="store_true", help="Save waveform overlay to output.png")

    args = parser.parse_args()
    print(f"DEBUG: args.output = '{args.output}'")

    # Remove any extension from output base name
    output_base = os.path.splitext(args.output)[0]

    wf_data = Utils.load_waveforms(args.npz_path)
    info = Utils.get_info(args.npz_path)

    single_peak_waveforms, avg_wf = waveform_analysis.store_avg(
        npz_path=args.npz_path,
        unit=args.unit,
        amplitude_threshold=args.amplitude_threshold,
        num_waveforms=args.num_waveforms,
        align=args.align,
        sampling=args.sampling
    )

    if not single_peak_waveforms:
        print("⚠️ No valid waveforms to analyze.")
        exit(1)

    # === Save NPZ (optional) ===
    if args.save_npz:
        filtered_npz = output_base + "_filtered.npz"
        median_npz = output_base + "_median_wf.npz"

        np.savez(filtered_npz, waveforms=np.array(single_peak_waveforms))
        np.savez(median_npz, median_waveform=np.array([avg_wf]))


        print(f"💾 Filtered waveforms saved to: {filtered_npz}")
        print(f"💾 Median waveform saved to: {median_npz}")
