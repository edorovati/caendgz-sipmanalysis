import os
import argparse
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import re

# Aggiungi il path corretto alla cartella 'python'
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis

def extract_bias_float(filename):
    match = re.search(r'vbias[_-]?(\d+\.?\d*)', filename)
    if match:
        return float(match.group(1))
    return None

def main():
    parser = argparse.ArgumentParser(description="Plot average-like waveforms da più file .npz (ordinati o etichettati)")
    parser.add_argument("input_path", type=str, help="Cartella contenente file .npz")
    parser.add_argument("--amplitude_threshold", type=float, required=True, help="Soglia minima per il picco")
    parser.add_argument("--sampling", type=float, required=True, help="Frequenza di campionamento (MHz)")
    parser.add_argument("--num_waveforms", type=int, required=True, help="Numero massimo di waveforms da processare per file")
    parser.add_argument("--unit", type=str, choices=["mV", "ADC"], default="mV", help="Unità di misura")
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    if not os.path.isdir(input_path):
        print(f"❌ Percorso non valido: {input_path}")
        return

    analyzer = waveform_analysis()
    entries_with_bias = []
    entries_no_bias = []

    # Loop ordinato sui file .npz
    npz_files = sorted(glob.glob(os.path.join(input_path, "*.npz")))
    for npz_file in npz_files:
        filename = os.path.basename(npz_file)
        bias = extract_bias_float(filename)
        label = f"{bias:.2f} V" if bias is not None else filename

        if bias is None:
            print(f"⚠️  Bias non trovato nel filename: {filename}")
        else:
            print(f"📂 Processing file: {filename} (bias = {label})")

        wf = analyzer.store_avg(
            npz_file,
            amplitude_threshold=args.amplitude_threshold,
            sampling=args.sampling,
            num_waveforms=args.num_waveforms,
            unit=args.unit,
            batch=False
        )
        if wf is not None and len(wf) > 0:
            if bias is not None:
                entries_with_bias.append((bias, wf, label))
            else:
                entries_no_bias.append((filename, wf, label))

    # Se ci sono bias, ordina per valore numerico
    if entries_with_bias:
        entries_with_bias.sort(key=lambda x: x[0])

    final_entries = entries_with_bias + entries_no_bias
    if not final_entries:
        print("❌ Nessuna waveform valida trovata.")
        return

    print(f"📊 Plotting {len(final_entries)} waveform(s)")

    # --- Asse temporale fisso a 1024 campioni, centrato automaticamente sul segnale ---
    nsamples = 1024
    margin = 10  # numero di sample da lasciare prima del picco

    # Usa la prima waveform valida per stimare la posizione del picco
    first_wf = final_entries[0][1][:nsamples]
    peak_index = np.argmax(first_wf)
    time_shift = peak_index - margin
    x_time = (np.arange(nsamples) - time_shift) / args.sampling * 1e3  # ns

    # --- Plot originale ---
    plt.figure(figsize=(10, 6))
    for _, wf, label in final_entries:
        plt.plot(x_time, wf[:nsamples], alpha=0.7, label=label)

    plt.xlabel("Time (ns)")
    plt.ylabel(f"Amplitude ({args.unit})")
    plt.title("Overlay waveforms")
    plt.grid(True)
    plt.legend(title="Label", loc="best")
    plt.tight_layout()
    plt.savefig("overlay_waveforms_ordered.png")
    print("✅ Plot salvato in 'overlay_waveforms_ordered.png'")

    # --- Plot normalizzato ---
    plt.figure(figsize=(10, 6))
    for _, wf, label in final_entries:
        segment = wf[:nsamples]
        norm = np.max(np.abs(segment))
        norm_wf = segment / norm if norm > 0 else segment
        plt.plot(x_time, norm_wf, alpha=0.7, label=label)

    plt.xlabel("Time (ns)")
    plt.ylabel("Normalized amplitude")
    plt.title("Overlay waveforms (normalized)")
    plt.grid(True)
    plt.legend(title="Label", loc="best")
    plt.tight_layout()
    plt.savefig("overlay_waveforms_ordered_normalized.png")
    print("✅ Plot salvato in 'overlay_waveforms_ordered_normalized.png'")

if __name__ == "__main__":
    main()
