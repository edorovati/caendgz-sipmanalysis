import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))

import argparse
from waveform_analysis import waveform_analysis  # dove sta la funzione estimate_tau_direct_debug
from utils import Utils  # per caricare waveform e info

def main():
    parser = argparse.ArgumentParser(description="Stima tau con debug e plot numerico.")
    parser.add_argument("npz_path", help="Percorso del file .npz contenente le waveform")
    parser.add_argument("--unit", default="mV", choices=["mV", "ADC"],
                        help="Unità di misura dell'ampiezza (default: mV)")
    parser.add_argument("--amplitude_threshold", type=float, default=4.0,
                        help="Soglia minima di picco per accettare una waveform [default: 4.0]")
    parser.add_argument("--num_waveforms", type=int, default=100,
                        help="Numero massimo di waveform da analizzare [default: 100]")
    parser.add_argument("--sampling", type=float, required=True,
                        help="Frequenza di campionamento in MHz")
    parser.add_argument("--output", type=str, required=True,
                        help="File di output per salvare i tau")
    args = parser.parse_args()

    # Carica le waveform
    wf_data = Utils.load_waveforms(args.npz_path)
    info = Utils.get_info(args.npz_path)
    baseline_range = (49, 973)  # default range per baseline

    # Seleziona waveform filtrate
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
        print("⚠️  Nessuna waveform supera la soglia P2P!")
        return

    # Stima dei tau + plot con debug numerico
    taus, amplitudes = waveform_analysis.estimate_tau_direct_debug(
        waveforms=selected_waveforms,
        sampling=args.sampling,
        single_pe_filter=True,
        threshold=args.amplitude_threshold,
        filename=args.npz_path,
        txt_output=args.output
    )

    if taus:
        print(f"✅ Stima completata: {len(taus)} waveform processate.")
        print(f"Media tau: {np.mean(taus):.3f} ns ± {np.std(taus):.3f} ns")
    else:
        print("⚠️  Nessuna stima di tau valida.")

if __name__ == "__main__":
    main()