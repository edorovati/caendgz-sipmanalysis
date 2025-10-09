import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from histograms import Histograms

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Plot histograms of timing distributions from .npz file(s).")
#     parser.add_argument("path", help="File .npz singolo o cartella contenente .npz")
#     parser.add_argument("--threshold", type=float, default=5.0, help="Soglia in mV o ADC")
#     parser.add_argument("--unit", choices=["mV", "ADC"], default="mV", help="Unità del segnale")
#     parser.add_argument("--mode", choices=["times", "differences", "peaks"], default="times", help="Tipo di distribuzione da plottare")
#     parser.add_argument("--subplot", action="store_true", help="Plottare in subplot invece che in un singolo plot")
#     parser.add_argument("--save", action="store_true", help="Salvare il grafico in PDF")
#     args = parser.parse_args()

#     files = []
#     if os.path.isfile(args.path) and args.path.endswith(".npz"):
#         files = [args.path]
#     elif os.path.isdir(args.path):
#         files = [os.path.join(args.path, f) for f in os.listdir(args.path) if f.endswith(".npz")]
#     else:
#         print(f"Errore: '{args.path}' non è un file .npz o una cartella valida.")
#         exit(1)

#     Histograms.plot_multiple(files, threshold=args.threshold, unit=args.unit, mode=args.mode, subplot=args.subplot, save=args.save)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms of tau distributions from .txt file(s).")
    parser.add_argument("path", help="File .txt singolo o cartella contenente .txt")
    parser.add_argument("--output_dir", default=".", help="Cartella di output per i grafici")
    parser.add_argument("--bins", type=int, default=100, help="Numero di bin per l'istogramma")
    args = parser.parse_args()
    
    files = []
    if os.path.isfile(args.path) and args.path.endswith(".txt"):
        files = [args.path]
    elif os.path.isdir(args.path):
        files = [os.path.join(args.path, f) for f in os.listdir(args.path) if f.endswith(".txt")]
    else:
        print(f"Errore: '{args.path}' non è un file .txt o una cartella valida.")
        exit(1)
        
    for file in files:
        print(f"Processing file: {file}")
        Histograms.plot_histograms_from_txt(file, output_dir=args.output_dir, bins=args.bins)
        print(f"Saved histograms to {args.output_dir}")
        
        