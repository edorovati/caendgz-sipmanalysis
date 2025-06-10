import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

# === Local imports ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from utils import Utils
from waveform_analysis import waveform_analysis
from plotting import Plotting

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute tau and overlay waveforms from an .npz file.")
    parser.add_argument("npz_path", help="Path to the .npz file.")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--amplitude_threshold", type=float, default=4.0)
    parser.add_argument("--num_waveforms", type=int, default=100)
    parser.add_argument("--color", default="green")
    parser.add_argument("--align", default="False")
    parser.add_argument("--sampling", type=float, required=True)
    parser.add_argument("--output", type=str, required=True, help="Output filename prefix (no extension)")
    parser.add_argument("--save_npz", action="store_true", help="Save filtered waveforms to output_name.npz")

    args = parser.parse_args()

    wf_data = Utils.load_waveforms(args.npz_path)
    info = Utils.get_info(args.npz_path)

    single_peak_waveforms, avg_wf = waveform_analysis.store_avg(
        npz_path=args.npz_path,
        unit=args.unit,
        amplitude_threshold=args.amplitude_threshold,
        num_waveforms=args.num_waveforms,
        align=args.align,
        sampling=args.sampling
    )

    if not single_peak_waveforms:
        print("⚠️ No valid waveforms to plot.")
        exit(1)

    waveform_analysis.compute_tau(single_peak_waveforms, avg_wf, args.sampling, args.output)

    # === Plot ===
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('#000c1f')

    Plotting.scope_plot(
        ax,
        single_peak_waveforms,
        unit=args.unit,
        npz_path=args.npz_path,
        color=(1.0, 1.0, 0.0),
        align=args.align,
        sampling=args.sampling,
        plot_avg=True
    )

    plt.tight_layout()
    fig.subplots_adjust(right=0.81)

    bias_value = info["bias"]
    filename_info = f'Bias voltage: {bias_value} V\n'
    filename_info += f"CAEN DGZ DT5742B\n"
    filename_info += f"Sampling rate: {args.sampling} GS\n"
    filename_info += f"Persistence: {"Auto" if args.align == "auto" else args.align}\n"
    fig.text(0.82, 0.5, filename_info, ha='left', va='center', fontsize=8, color='white', fontfamily='monospace')

    # === SAVE PNG ===
    output_png = args.output + ".png"
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=1200)
    print(f"✅ Plot saved to: {output_png}")

    # === SAVE NPZ ===
    if args.save_npz:
        output_npz = args.output + ".npz"
        np.savez(output_npz, waveforms=np.array(single_peak_waveforms), average=avg_wf)
        print(f"💾 Saved filtered waveforms to: {output_npz}")
