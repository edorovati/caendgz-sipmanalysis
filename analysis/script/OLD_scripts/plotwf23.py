import os
import argparse
import sys

# Aggiungi il path corretto alla cartella 'python'
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis

def main():
    parser = argparse.ArgumentParser(description="Estrai waveform media da un singolo file .npz")
    parser.add_argument("npz_file", type=str, help="Path al file .npz")
    parser.add_argument("--amplitude_threshold", type=float, required=True, help="Soglia minima per il picco")
    parser.add_argument("--sampling", type=float, required=True, help="Frequenza di campionamento (MHz)")
    parser.add_argument("--num_waveforms", type=int, required=True, help="Numero massimo di waveform da processare")
    parser.add_argument("--unit", type=str, choices=["mV", "ADC"], default="mV", help="Unità di misura")
    parser.add_argument("--batch", action="store_true", help="Esegui in modalità batch (opzionale)")
    parser.add_argument("--output", type=str, help="Path per il file di output (opzionale)")
    args = parser.parse_args()

    analyzer = waveform_analysis()
    wf = analyzer.store_avg(
        npz_path=args.npz_file,
        amplitude_threshold=args.amplitude_threshold,
        sampling=args.sampling,
        num_waveforms=args.num_waveforms,
        unit=args.unit,
        batch=True
    )

    if wf is not None and len(wf) > 0:
        print("✅ Waveform media estratta con successo.")
        # Volendo, puoi salvarla o passarla ad altro codice qui
    else:
        print("⚠️ Nessuna waveform valida trovata.")

if __name__ == "__main__":
    main()
