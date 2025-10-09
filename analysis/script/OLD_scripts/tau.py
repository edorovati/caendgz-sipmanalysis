import numpy as np
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))

import argparse
from waveform_analysis import waveform_analysis
from utils import Utils

def main():
    parser = argparse.ArgumentParser(description="Estimate tau using 1/e decay method from waveform data (.npz).")

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
    parser.add_argument("--output", type=str, required=True,
                        help="Output file name")

    args = parser.parse_args()

    wf_data = Utils.load_waveforms(args.npz_path)
    info = Utils.get_info(args.npz_path)
    baseline_range = (49, 973)

    selected_waveforms, _ = waveform_analysis.analyze_waveforms(
        wf_data=wf_data,
        fs=args.sampling,
        baseline_range=baseline_range,
        amplitude_threshold=args.amplitude_threshold,
        max_waveforms=args.num_waveforms,
        unit=args.unit,
        get_edges=True
    )

    if not selected_waveforms:
        print("⚠️  No waveforms passed the P2P threshold!")
        return

#    txt_output = os.path.join(
#        os.path.dirname(args.npz_path),
#        os.path.basename(args.npz_path).replace(".npz", "_direct_tau.txt")
#    )

    txt_output = args.output
    
    taus = waveform_analysis.estimate_tau_direct(
        waveforms=selected_waveforms,
        sampling=args.sampling,
        single_pe_filter=args.single_pe_peak,
        threshold=args.amplitude_threshold,
        filename=args.npz_path,
        txt_output=txt_output
    )

    if taus:
        print(f"✅ Direct tau estimation successful: {len(taus)} waveforms processed")
        print(f"Average tau: {np.mean(taus):.2f} ns ± {np.std(taus):.2f} ns")
    else:
        print("⚠️  No valid tau estimations.")

if __name__ == "__main__":
    main()
