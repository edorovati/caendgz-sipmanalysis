import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Compare histograms of selected columns from multiple txt files.")
parser.add_argument("input", type=str, help="Input txt file or glob pattern (e.g. data/*.txt)")
parser.add_argument("--labels", nargs="+", type=str, required=True, help="Labels for each selected column")
parser.add_argument("--cols", nargs="+", type=int, required=True, help="Columns to plot (zero-based)")
parser.add_argument("--bins", type=int, default=30, help="Number of bins for the histograms")
parser.add_argument("--save", type=str, help="Output image filename (optional)")
args = parser.parse_args()

files = sorted(glob.glob(args.input))
if not files:
    raise FileNotFoundError(f"No files match pattern: {args.input}")

if len(args.labels) != len(args.cols):
    raise ValueError("Number of labels must match number of selected columns.")

# Read and organize data
data_per_file = []
for file in files:
    try:
        data = np.loadtxt(file)
        data_per_file.append(data)
    except Exception as e:
        print(f"Skipping {file}: {e}")

# Plot
fig, axs = plt.subplots(len(args.cols), 1, figsize=(7, 4 * len(args.cols)), sharex=False)
if len(args.cols) == 1:
    axs = [axs]

colors = plt.cm.get_cmap('Set2', len(files))  # colori distinguibili

for idx_col, (col_idx, label) in enumerate(zip(args.cols, args.labels)):
    ax = axs[idx_col]

    # Calcolo il range globale per tutti i file su questa colonna
    all_values = np.concatenate([data[:, col_idx] for data in data_per_file])
    vmin = np.min(all_values)
    vmax = np.max(all_values)
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5
    bin_edges = np.linspace(vmin, vmax, args.bins + 1)


    for i, file_data in enumerate(data_per_file):
        if col_idx >= file_data.shape[1]:
            raise ValueError(f"Column {col_idx} out of bounds in file {files[i]}")
        values = file_data[:, col_idx]
        ax.hist(values, bins=bin_edges, density=True, alpha=0.6, label=f"{files[i]}", color=colors(i))

    ax.set_title(f"Normalized Histogram: {label}")
    ax.set_xlabel(label)
    ax.set_ylabel("Normalized Counts")
    ax.legend()
    ax.grid(True)

plt.tight_layout()

if args.save:
    plt.savefig(args.save)
    print(f"[INFO] Plot saved as: {args.save}")
else:
    plt.show()
