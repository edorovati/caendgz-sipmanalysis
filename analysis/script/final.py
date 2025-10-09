import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import re

def extract_vbias(filename):
    match = re.search(r'tau_(\d+)(?:_(\d+))?-run\d+\.txt', filename)
    if match:
        integer_part = int(match.group(1))
        decimal_part = match.group(2)
        if decimal_part:
            vbias = float(f"{integer_part}.{decimal_part}")
        else:
            vbias = float(integer_part)
        return vbias
    else:
        raise ValueError(f"Nome file non riconosciuto: {filename}")

def process_folder(folder):
    search_pattern = os.path.join(folder, 'tau_*-run0.txt')
    file_list = sorted(glob.glob(search_pattern))

    if not file_list:
        print(f"Nessun file trovato nella cartella {folder}")
        return None, None

    vbias_values = []
    tau_values = []

    for file in file_list:
        try:
            vbias = extract_vbias(os.path.basename(file))
            # Usa genfromtxt per ignorare header, prendi solo la prima riga dati
            data = np.genfromtxt(file, skip_header=1, max_rows=1)

            print(f"DEBUG: file={file}, data={data}")

            if data.size == 0:
                print(f"Attenzione: file {file} vuoto dopo header, salto")
                continue

            # Se data è 1D array o singolo valore, prendi la prima colonna
            if data.ndim == 0:
                tau = data
            else:
                tau = data[0]

            vbias_values.append(vbias)
            tau_values.append(tau)

        except Exception as e:
            print(f"Errore con il file {file}: {e}")

    if not vbias_values:
        return None, None

    vbias_values = np.array(vbias_values)
    tau_values = np.array(tau_values)
    sorted_indices = np.argsort(vbias_values)

    return vbias_values[sorted_indices], tau_values[sorted_indices]


def main(folders):
    plt.figure(figsize=(10, 6))
    
    for folder in folders:
        vbias, tau = process_folder(folder)
        if vbias is None or tau is None:
            continue
        
        label = os.path.basename(os.path.normpath(folder))
        plt.plot(vbias, tau, marker='o', linestyle='-', label=label)
    
    plt.ylim(0, 70)
    plt.xlabel('Vbias (V)')
    plt.ylabel('Tau')
    plt.title('Tau vs Vbias')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./tau_vs_vbias.png')
    print("Plot salvato in 'tau_vs_vbias.png'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot tau vs vbias from text files in multiple folders')
    parser.add_argument('folders', nargs='+', type=str, help='Directories containing the txt files')
    args = parser.parse_args()
    main(args.folders)
