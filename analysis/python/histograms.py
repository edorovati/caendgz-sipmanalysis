import numpy as np
import matplotlib.pyplot as plt
from utils import Utils
from waveform_analysis import waveform_analysis
from scipy.optimize import curve_fit
import os

class Histograms:
    @staticmethod
    def get_above_threshold_times(waveforms, threshold):
        baseline = np.mean(waveforms, axis=1, keepdims=True)
        wf_centered = waveforms - baseline
        times = []
        for wf in wf_centered:
            above = np.where(wf > threshold)[0]
            times.extend(above.tolist())
        return np.array(times)

    @staticmethod
    def get_time_differences(times, npz_file):
        times_s = Utils.cell_to_seconds(times, npz_file)
        diffs = np.diff(np.sort(times_s)) * 1e9  # ns
        return diffs

    @staticmethod
    def load_times(npz_file, threshold, unit, mode="times"):
        data = Utils.load_waveforms(npz_file)
        waveforms = data["waveforms"]
        if unit == "mV":
            waveforms = Utils.adc_to_mv(waveforms)
        times = Histograms.get_above_threshold_times(waveforms, threshold)
        if mode == "times":
            return Utils.cell_to_seconds(times, npz_file) * 1e9  # ns
        elif mode == "differences":
            return Histograms.get_time_differences(times, npz_file)
        # let's make an option also for peaks with the "rising_edge" method
        elif mode == "peaks":
            peaks = []
            for wf in waveforms:
                peak_index = waveform_analysis.find_rising_edge(wf, threshold)
                if peak_index != 0:  # se vuoi escludere casi senza fronte
                    peaks.append(peak_index)
            return Utils.cell_to_seconds(np.array(peaks), npz_file) * 1e9
        else:
            raise ValueError(f"Unknown mode: {mode}")

    @staticmethod
    def plot_multiple(files, threshold=5.0, unit='mV', mode='times', subplot=False, save=False):
        """
        Plot multiple distributions from list of .npz files.
        mode: 'times' or 'differences'
        """
        n = len(files)
        fig, axes = plt.subplots(n if subplot else 1, 1, figsize=(8, 3*n if subplot else 4), squeeze=False)
        axes = axes.flatten()

        for idx, file in enumerate(files):
            label = f"{Utils.get_info(file)["bias"]}"
            data = Histograms.load_times(file, threshold, unit, mode)
            ax = axes[idx] if subplot else axes[0]
            ax.hist(data, bins=200, alpha=0.6, label=label)
            ax.set_xlabel("Time [ns]" if mode == "times" else "ΔTime [ns]")
            ax.set_ylabel("Counts")
            ax.grid(True)
            if subplot:
                ax.set_title(f"{label}")
            else:
                ax.legend()

        axes[0].set_title("Sample times above threshold" if mode == "times" else "Time differences above threshold")
        plt.tight_layout()
        if save:
            filename = f"hist_{mode}_{unit}.pdf"
            plt.savefig(filename)
            print(f"Saved to {filename}")
        plt.show()
        
        
    #static method to plot RMS noise histogram where there is the baseline only
    @staticmethod
    def plot_rms_noise_histogram(npz_file, threshold=5.0, unit='mV', save=False):
        """
        Plot RMS noise histogram from .npz file.
        """
        data = Utils.load_waveforms(npz_file)
        waveforms = data["waveforms"]
        if unit == "mV":
            waveforms = Utils.adc_to_mv(waveforms)
        baseline = np.mean(waveforms, axis=1, keepdims=True)
        wf_centered = waveforms - baseline
        rms_noise = np.std(wf_centered, axis=1)
        
        plt.figure(figsize=(8, 4))
        plt.hist(rms_noise, bins=200, alpha=0.6)
        plt.xlabel("RMS Noise [mV]" if unit == "mV" else "RMS Noise [ADC]")
        plt.ylabel("Counts")
        plt.grid(True)
        plt.title("RMS Noise Histogram")
        
        if save:
            filename = f"rms_noise_histogram_{unit}.pdf"
            plt.savefig(filename)
            print(f"Saved to {filename}")
        
        plt.show()
        
    @staticmethod
    def load_times(npz_file, threshold, unit, mode="times"):
        data = Utils.load_waveforms(npz_file)
        waveforms = data["waveforms"]
        if unit == "mV":
            waveforms = Utils.adc_to_mv(waveforms)

        if mode == "times":
            times = Histograms.get_above_threshold_times(waveforms, threshold)
            return Utils.cell_to_seconds(times, npz_file) * 1e9  # ns

        elif mode == "differences":
            times = Histograms.get_above_threshold_times(waveforms, threshold)
            return Histograms.get_time_differences(times, npz_file)

        elif mode == "peaks":
            peaks = []
            for wf in waveforms:
                peak_index = waveform_analysis.find_rising_edge(wf, threshold)
                if peak_index != 0:
                    peaks.append(peak_index)
            peaks = np.array(peaks)
            times_ns = Utils.cell_to_seconds(peaks, npz_file) * 1e9
            # Filtro: escludi picchi prima di 5 ns
            filtered_times = times_ns[times_ns > 5]
            return filtered_times

        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        
    def plot_histograms_from_txt(txt_path, output_dir=".", bins=100):
        # Carica header (prima riga) e dati
        with open(txt_path, "r") as f:
            header = f.readline().strip().split()

        data = np.loadtxt(txt_path, skiprows=1)

        if data.ndim == 1:
            data = data[:, np.newaxis]  # Caso di una sola colonna

        # Verifica o crea cartella output
        os.makedirs(output_dir, exist_ok=True)

        # Istogramma per ogni colonna
        for i, label in enumerate(header):
            plt.figure()
            plt.hist(data[:, i], bins=bins, color="steelblue", edgecolor="black")
            plt.xlabel(label)
            plt.ylabel("Counts")
            plt.title(f"Istogramma di {label}")
            plt.grid(True)
            plt.tight_layout()
            out_path = os.path.join(output_dir, f"hist_{label}.png")
            plt.savefig(out_path, dpi=300)
            plt.close()
            print(f"✅ Salvato: {out_path}")


