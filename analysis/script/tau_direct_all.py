import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

# === Local imports ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis
from plotting import Plotting

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute tau from a .npz file of waveforms (no filtering).")
    parser.add_argument("npz_path", help="Path to the .npz file with waveforms")
    parser.add_argument("--sampling", type=float, required=True, help="Sampling rate in GS/s")
    parser.add_argument("--output", type=str, required=True, help="Base name for output files (no extension)")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--save_plot", action="store_true", help="Save overlay plot of waveforms")

    args = parser.parse_args()

    output_base = os.path.splitext(args.output)[0]
    data = np.load(args.npz_path)

    # Load waveforms: either "waveforms" or "median_waveform"
    if "waveforms" in data:
        waveforms = data["waveforms"]
    else:
        print("❌ .npz file does not contain 'waveforms' or 'median_waveform'")
        exit(1)

    if len(waveforms) == 0:
        print("⚠️ No waveform found in input npz.")
        exit(1)

    median_wf = np.median(waveforms, axis=0)

    # === Compute tau ===
    waveform_analysis.compute_tau(waveforms, median_wf, args.sampling, output_base + ".txt")
    print(f"📄 Tau saved to: {output_base}.txt")
