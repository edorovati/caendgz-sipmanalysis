import os
import argparse
#let's import the path of the waveform_analysis class in ../script/
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis
from utils import Utils  # Importa la classe Utils per trovare la cartella 'data'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Overlay filtered waveforms above the P2P threshold for .npz file(s).")
    parser.add_argument("path", help="Single .npz file or folder inside ./data/")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--num_waveforms", type=int, default=100)
    parser.add_argument("--p2p_threshold", type=float, default=3.5)
    parser.add_argument("--color", default="green")
    args = parser.parse_args()

    data_dir = Utils.find_data_directory()
    full_path = os.path.join(data_dir, args.path)

    analyzer = waveform_analysis()

    if os.path.isfile(full_path) and full_path.endswith(".npz"):
        print(f"Processing single file: {args.path}")
        analyzer.persistence(
            full_path,
            unit=args.unit,
            p2p_threshold=args.p2p_threshold,
            num_waveforms=args.num_waveforms,
            color=args.color
        )
    elif os.path.isdir(full_path):
        for filename in os.listdir(full_path):
            if filename.endswith(".npz"):
                file_path = os.path.join(full_path, filename)
                print(f"Processing file: {filename}")
                analyzer.persistence(
                    file_path,
                    unit=args.unit,
                    p2p_threshold=args.p2p_threshold,
                    num_waveforms=args.num_waveforms,
                    color=args.color
                )
    else:
        print(f"'{full_path}' is not a valid .npz file or folder.")
