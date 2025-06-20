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
    

    # === Load waveform(s) ===
    if "waveforms" in data:
        waveforms = data["waveforms"]
    elif "median_waveform" in data:
        waveforms = data["median_waveform"]
    else:
        print("❌ .npz file must contain either 'waveforms' or 'median_waveform'")
        exit(1)

    if waveforms is None or waveforms.size == 0:
        print("⚠️ No waveform found in input npz.")
        exit(1)

    
    print(f"📂 Loaded {len(waveforms)} waveforms from {args.npz_path}")
    # === Compute tau using FIT! ===> then fit distribution etc etc
    taus, tau_errors, amps, amp_errors, offsets, offset_errors, fitted_waveforms = waveform_analysis.fit_waveform_tau(
        waveforms,
        sampling=args.sampling,
        filename=args.output
    )

    
    
    print(f"📄 Tau-fit saved to: {output_base}.txt")
