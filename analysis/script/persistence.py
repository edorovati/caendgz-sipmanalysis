import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

# === Local imports ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from waveform_analysis import waveform_analysis
from utils import Utils
from plotting import Plotting

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Persistence plot from filtered waveform set.")
    parser.add_argument("path", help="Path to .npz file")
    parser.add_argument("--unit", choices=["mV", "ADC"], default="mV")
    parser.add_argument("--num_waveforms", type=int, default=100)
    parser.add_argument("--amplitude_threshold", type=float, default=4.0)
    parser.add_argument("--sampling", type=float, required=True, help="Sampling rate in MHz")
    parser.add_argument("--align", default="False", help="Alignment method: False, 'max', or 'auto'")
    parser.add_argument("--output", type=str, required=True, help="Output filename prefix (no extension)")
    parser.add_argument("--save_npz", action="store_true", help="Save filtered waveforms to output_name.npz")

    args = parser.parse_args()

    npz_path = args.path
    sampling = args.sampling
    unit = args.unit
    align = args.align
    num_waveforms = args.num_waveforms
    amplitude_threshold = args.amplitude_threshold

    wf_data = Utils.load_waveforms(npz_path)
    info = Utils.get_info(npz_path)

    # === Analyze waveforms ===
    baseline_range = (49, 973)
    waveforms, rising_edges = waveform_analysis.analyze_waveforms(
        wf_data=wf_data,
        fs=sampling,
        baseline_range=baseline_range,
        amplitude_threshold=amplitude_threshold,
        max_waveforms=num_waveforms,
        unit=unit,
        get_edges=True
    )

    if not waveforms:
        print("⚠️ No waveforms passed the filter.")
        exit()

    # === Alignment ===
    if align and align != "False":
        waveforms = waveform_analysis.align_waveforms(waveforms, rising_edges, mode=align)

    if not waveforms:
        print("⚠️ No valid waveforms to plot.")
        exit()

    # === Plot ===
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('#000c1f')

    alpha_decay = np.linspace(0.05, 0.8, len(waveforms))
    Plotting.scope_plot(ax, waveforms, unit=unit, npz_path=npz_path, color=(1.0, 1.0, 0.0),
                        align=align, sampling=sampling, plot_avg=True)

    plt.tight_layout()
    Plotting.add_inset_zoom(ax, waveforms, alpha_decay, default_target_idx=150 if align else None)
    fig.subplots_adjust(right=0.81)

    bias_value = info["bias"]
    filename_info = f"Bias voltage: {bias_value} V\n"
    filename_info += f"CAEN DGZ DT5742B\n"
    filename_info += f"Sampling rate: {sampling} GS\n"
    filename_info += f"Persistence: {'Auto' if align == 'auto' else align}\n"
    fig.text(0.82, 0.5, filename_info, ha='left', va='center', fontsize=8, color='white', fontfamily='monospace')

    # === SAVE PNG ===
    output_png = args.output + ".png"
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=1200)
    print(f"✅ Plot saved to: {output_png}")

    # === SAVE NPZ ===
    if args.save_npz:
        output_npz = args.output + ".npz"
        np.savez(output_npz, waveforms=np.array(waveforms))
        print(f"💾 Saved filtered waveforms to: {output_npz}")
