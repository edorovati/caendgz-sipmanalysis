import os
import argparse
import sys

# aggiungo il path corretto per waveform_analysis
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis
from utils import Utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estimate tau from .npz file(s).")
    parser.add_argument("path", help="Single .npz file or folder inside ./data/")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--num_waveforms", type=int, default=100)
    parser.add_argument("--amplitude_threshold", type=float, default=4.5)
    parser.add_argument("--color", default="green")
    parser.add_argument("--sampling", type=float, required=True, help="Sampling rate in MHz")
    parser.add_argument("--output_filename", type=str, default=None, help="Nome file di output TXT per i risultati tau")

    args = parser.parse_args()
    full_path = args.path

    analyzer = waveform_analysis()

    if os.path.isfile(full_path) and full_path.endswith(".npz"):
        print(f"Processing single file: {args.path}")
        analyzer.tau_direct_mean(
            full_path,
            unit=args.unit,
            amplitude_threshold=args.amplitude_threshold,
            num_waveforms=args.num_waveforms,
            color=args.color,
            sampling=args.sampling,
            output_filename=args.output_filename
        )
    elif os.path.isdir(full_path):
        for filename in os.listdir(full_path):
            if filename.endswith(".npz"):
                file_path = os.path.join(full_path, filename)
                print(f"Processing file: {filename}")
                analyzer.tau_direct_mean(
                    file_path,
                    unit=args.unit,
                    amplitude_threshold=args.amplitude_threshold,
                    num_waveforms=args.num_waveforms,
                    color=args.color,
                    sampling=args.sampling,
                    output_filename=args.output_filename
                )
    else:
        print(f"'{full_path}' is not a valid .npz file or folder.")