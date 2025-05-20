import numpy as np
import matplotlib.pyplot as plt
from utils import Utils
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class Plotting:

    @staticmethod
    def scope_plot(ax, waveforms, unit, npz_path, color=(1.0, 1.0, 0.0), align='auto'):
        """
        Plot waveforms in oscilloscope-style overlay with time axis.

        Parameters:
        - ax: matplotlib axis object
        - waveforms: list of waveforms to plot
        - unit: string, e.g. 'mV'
        - npz_path: path to .npz file (used to extract sampling rate)
        - color: base color for waveform lines
        - align: alignment mode for title
        """
        if not waveforms:
            print("⚠️ No waveforms to plot in scope_plot.")
            return

        alpha_decay = np.linspace(0.05, 0.8, len(waveforms))
        x_time = Utils.cell_to_seconds(np.arange(len(waveforms[0])), npz_path) * 1e9

        for wf, alpha in zip(waveforms, alpha_decay):
            ax.plot(x_time, wf, color=color, alpha=alpha, linewidth=1)

        ax.set_xlabel("Time [ns]")
        ax.set_ylabel(f"Amplitude [{unit}]")
        ax.set_title(f'Oscilloscope-style Overlay of {len(waveforms)} {"Aligned" if align else "Raw"} Waveforms ({unit})')
        ax.grid(color='white', alpha=0.2)
        
    
    @staticmethod
    def add_inset_zoom(ax, waveforms, alpha_decay, base_color=(1.0, 1.0, 0.0), default_target_idx=150):
        # Use "central" frames for signal research (100–850)
        peak_idxs = [np.argmax(wf) for wf in waveforms]
        valid_peaks = [idx for idx in peak_idxs if 150 <= idx <= 850]

        if valid_peaks:
            hist, bins = np.histogram(valid_peaks, bins=np.arange(150, 851, 10))
            max_bin_idx = np.argmax(hist)
            target_idx = int((bins[max_bin_idx] + bins[max_bin_idx + 1]) // 2)
        else:
            target_idx = default_target_idx

        zoom_range = slice(target_idx - 50, target_idx + 100)
        p2p_global = np.ptp(np.concatenate(waveforms))

        inset_ax = inset_axes(ax, width="30%", height="30%", loc='upper right', borderpad=2)
        inset_ax.set_facecolor("#111122")
        inset_ax.grid(color='white', alpha=0.2)

        for wf, alpha in zip(waveforms, alpha_decay):
            inset_ax.plot(np.arange(len(wf))[zoom_range], wf[zoom_range], color=base_color, alpha=alpha, linewidth=1)

        if p2p_global <= 10:
            inset_ax.set_ylim(-2, 4)
            inset_ax.set_title("Zoom for low signals", fontsize=10)
        elif p2p_global <= 30:
            inset_ax.set_ylim(-4, 8)
            inset_ax.set_title("Zoom for high signals", fontsize=10)
        else:
            inset_ax.set_ylim(-15, 50)
            inset_ax.set_title("Zoom for very high signals", fontsize=10)
        inset_ax.set_xlim(zoom_range.start, zoom_range.stop)
            
        inset_ax.tick_params(axis='both', which='both', labelsize=8)