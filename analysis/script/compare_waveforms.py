import os
import argparse
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis
from utils import Utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare aligned and normalized waveforms from multiple .npz files.")
    parser.add_argument("files", nargs='+', help="List of .npz files with full or relative paths.")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--num_waveforms", type=int, default=100)
    parser.add_argument("--amplitude_threshold", type=float, default=4.5)
    parser.add_argument("--sampling", type=float, required=True, help="Sampling rate in GHz")

    args = parser.parse_args()

    full_paths = args.files  # usa i percorsi così come li passi in input

    analyzer = waveform_analysis()
    analyzer.compare_waveforms(
        npz_paths=full_paths,
        unit=args.unit,
        amplitude_threshold=args.amplitude_threshold,
        num_waveforms=args.num_waveforms,
        sampling=args.sampling
    )
