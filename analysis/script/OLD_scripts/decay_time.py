# import os
# import argparse
# #let's import the path of the waveform_analysis class in ../script/
import numpy as np
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))

import argparse
from waveform_analysis import waveform_analysis

def main():
    parser = argparse.ArgumentParser(description="Decay time analysis from waveform data (.npz).")

    parser.add_argument("npz_path", help="Path to the input .npz file")
    parser.add_argument("--unit", default="mV", choices=["mV", "ADC"],
                        help="Unit to use for waveform amplitude (default: mV)")
    parser.add_argument("--amplitude_threshold", type=float, default=4.0,
                        help="Minimum peak amplitude to accept a waveform [default: 4.0]")
    parser.add_argument("--num_waveforms", type=int, default=100,
                        help="Maximum number of waveforms to analyze [default: 100]")
    parser.add_argument("--color", default="green",
                        help="Color for plotting (e.g., green, red, yellow) [default: green]")
    parser.add_argument("--single-pe-peak", action="store_true",
                        help="Filter to select only waveforms with a single peak ≥ threshold")
    parser.add_argument("--align", default="False",
                        help="Alignment mode: 'False', 'auto', or index value [default: False]")
    parser.add_argument("--sampling", type=float, required=True,
                        help="Sampling rate in MHz (required)")

    # Nuovi parametri
    parser.add_argument("--fit-window", type=int, default=50,
                        help="Max number of points to fit after peak [default: 50]")
    parser.add_argument("--tail-sigma-cut", type=float, default=2.0,
                        help="Stop fit when signal drops below baseline + Nσ [default: 2.0]")
    parser.add_argument("--baseline-window", type=int, default=50,
                        help="Number of points to compute baseline before peak [default: 50]")

    args = parser.parse_args()

    result = waveform_analysis.decay_time(
        npz_path=args.npz_path,
        unit=args.unit,
        amplitude_threshold=args.amplitude_threshold,
        num_waveforms=args.num_waveforms,
        color=args.color,
        align=args.align if args.align != "False" else False,
        single_pe_filter=args.single_pe_peak,
        sampling=args.sampling,
        fit_window=args.fit_window,
        tail_sigma_cut=args.tail_sigma_cut,
        baseline_window=args.baseline_window
    )

    if result:
        taus, tau_errors, amps, amp_errors, offsets, offset_errors = result
        print(f"✅ Decay time fit successful: {len(taus)} waveforms fitted")
        print(f"Average tau: {np.mean(taus):.2f} ns ± {np.std(taus):.2f} ns")
        print(f"Average A:   {np.mean(amps):.2f} ± {np.std(amps):.2f}")
        print(f"Average C:   {np.mean(offsets):.2f} ± {np.std(offsets):.2f}")
    else:
        print("⚠️  No successful decay time fits.")

if __name__ == "__main__":
    main()
