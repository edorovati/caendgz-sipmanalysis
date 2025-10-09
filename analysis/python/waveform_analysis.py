import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import time
from sklearn.cluster import DBSCAN
import re
from utils import Utils
from filters import Filters
from plotting import Plotting
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from matplotlib.backends.backend_pdf import PdfPages






class waveform_analysis:

    
    ###########################
    # BASIC FUNCTIONS
    ###########################

    @staticmethod
    def calculate_baseline_with_mask(
        wf, baseline_start, baseline_end,
        window_size=40, std_threshold=1.5,
        min_good_block=30, pre_margin=30, post_margin=150
    ):
        """
        Estimate baseline by excluding regions with signal-like fluctuations,
        based on standard deviation in sliding windows.

        Returns:
        - baseline value
        - indices used for baseline estimation
        """
        region = wf[baseline_start:baseline_end]
        num_points = len(region)
        region_indices = np.arange(baseline_start, baseline_end)

        # Compute std in sliding windows
        stds = np.array([
            np.std(region[i:i+window_size])
            for i in range(num_points - window_size)
        ])

        # Identify noisy windows
        noisy = stds > std_threshold
        noisy_indices = np.where(noisy)[0] + baseline_start

        # Create exclusion mask
        exclude_mask = np.zeros_like(wf, dtype=bool)
        for idx in noisy_indices:
            start = max(0, idx - pre_margin)
            end = min(len(wf), idx + post_margin)
            exclude_mask[start:end] = True

        # Identify good indices
        usable_indices = region_indices[~exclude_mask[baseline_start:baseline_end]]

        # Group into blocks and filter
        if len(usable_indices) > 0:
            diffs = np.diff(usable_indices)
            block_edges = np.where(diffs > 1)[0]
            blocks = np.split(usable_indices, block_edges + 1)
            good_blocks = [b for b in blocks if len(b) >= min_good_block]
            final_indices = np.concatenate(good_blocks) if good_blocks else np.array([])
        else:
            final_indices = np.array([])

        if len(final_indices) > 0:
            baseline = np.median(wf[final_indices])
        else:
            baseline = np.median(region)
            final_indices = np.array([])

        return baseline, final_indices


    @staticmethod
    def group_consecutive(indices):
        """
        Group consecutive indices together to create contiguous spans.
        Useful for visualizing ranges on plots.
        """
        if len(indices) == 0:
            return []
        groups = []
        start = indices[0]
        for i in range(1, len(indices)):
            if indices[i] != indices[i-1] + 1:
                groups.append((start, indices[i-1]))
                start = indices[i]
        groups.append((start, indices[-1]))
        return groups
    
    
    # should this be changed in favor of find_peaks?
    # SASSANDUM EST?
    @staticmethod
    def find_rising_edge(wf, threshold_factor=0.63):
        peak = np.max(wf)
        threshold = peak * threshold_factor
        for i in range(128, 768): #  should be changed dinamically 
            if wf[i] > threshold and wf[i - 1] <= threshold:
                return i
        return 0

    @staticmethod
    def exp_decay(t, A, tau, C):
        return A * np.exp(-t / tau) + C
    

    # SASSANDUM EST! NOT USED ANYMORE!
    # @staticmethod
    # def find_polarity(waveform):
    #     """
    #     Determines the dominant signal polarity of a baseline-subtracted waveform.
    #     Returns 'positive' or 'negative'.
    #     """
    #     max_val = np.max(waveform)
    #     min_val = np.min(waveform)

    #     if abs(min_val) > abs(max_val):
    #         return 'negative'
    #     else:
    #         return 'positive'
       
    
    
    ###########################
    # WAVEFORM PROCESSING POST-BASELINE SUBTRACTION
    ###########################
    
    @staticmethod 
    def align_waveforms(waveforms, rising_edges, mode="True"):
        # consider removing the "auto" mode, not perfectly controllable
        # IMPORTANT: ADD the option to set the offset to 0 and shift the waveforms accordingly?   
           # maybe can be handled in the plotting function 
        """
        Aligns waveforms to their rising edge positions.
        
        Parameters:
            waveforms (list of np.ndarray): The list of waveforms to align.
            rising_edges (list of int): The list of rising edge positions.
            mode (str or bool): 'auto' to detect alignment necessity via clustering,
                                True to force alignment,
                                False to skip alignment.
        Returns:
            list of np.ndarray: Aligned or original waveforms depending on the mode.
        """
        if mode is False:
            return waveforms

        if mode == "auto":
            # Cluster rising edges to detect if waveforms are already aligned
            edge_positions = np.array(rising_edges).reshape(-1, 1)
            db = DBSCAN(eps=4, min_samples=10).fit(edge_positions)
            labels = db.labels_

            # Count number of waveforms in the largest cluster (excluding noise)
            unique_labels, counts = np.unique(labels[labels != -1], return_counts=True)
            max_cluster_size = counts.max() if len(counts) > 0 else 0

            if max_cluster_size >= 10: # this introduces some bias depending on the number of spread waveforms (>100) --> FIX THIS
                print(f"[INFO] Waveforms are already aligned (cluster size = {max_cluster_size}), skipping alignment")
                return waveforms
            else:
                print(f"[INFO] No dominant cluster found (max cluster size = {max_cluster_size}), aligning waveforms")

        # Perform alignment (either auto detected or forced)
        target_idx = 150  # Desired target index for alignment
        aligned = []

        for wf, edge in zip(waveforms, rising_edges):
            shift = target_idx - edge
            if shift > 0:
                aligned_wf = np.pad(wf, (shift, 0), mode='constant')[:len(wf)]
            elif shift < 0:
                aligned_wf = np.roll(wf, shift)
                aligned_wf[shift:] = 0
            else:
                aligned_wf = wf.copy()
            aligned.append(aligned_wf)

        return aligned



    # iper genera-purpose analyzer, included: baseline subtraction, peak detection, alignment, etc.
    @staticmethod
    def analyze_waveforms(wf_data, fs, baseline_range, amplitude_threshold=5.5,
                        max_waveforms=None, unit='mV', get_edges=False,
                        peak_index_range=(128, 768), single_pe_peak=False):
        """
        Filters waveform data based on peak amplitude and peak position after baseline subtraction.
        
        Parameters:
            wf_data (npz): loaded waveform file
            fs (float): sampling rate
            baseline_range (tuple): (start, end) indices for baseline estimation
            amplitude_threshold (float): min peak amplitude to keep waveform
            max_waveforms (int or None): maximum number of waveforms to return (None = no limit)
            unit (str): 'ADC' or 'mV'
            get_edges (bool): if True, return rising edges too
            peak_index_range (tuple): (min_index, max_index) allowed range for the peak position
            single_pe_peak (bool): if True, keep only waveforms with a single peak above threshold

        Returns:
            selected_waveforms: list of filtered waveforms
            rising_edges (optional): list of rising edge indices
        """
        baseline_start, baseline_end = baseline_range
        peak_min_index, peak_max_index = peak_index_range
        selected_waveforms = []
        rising_edges = []
        
        total_checked = 0
        total_kept = 0

        min_peak_distance_samples = int(fs * 5e-9)  # ~5 ns
        min_peak_height = 4.0  # mV, per il filtro single_pe

        for key in wf_data.files:
            wf_array = wf_data[key]
            if wf_array.ndim == 2 and wf_array.shape[0] > 0:
                for i in range(wf_array.shape[0]):
                    wf_adc = wf_array[i]
                    total_checked += 1

                    wf = wf_adc
                    if unit == 'mV':
                        wf = Utils.adc_to_mv(wf)

                    baseline, _ = waveform_analysis.calculate_baseline_with_mask(
                        wf, baseline_start, baseline_end
                    )
                    wf_corrected = wf - baseline
                    amplitude = np.max(wf_corrected)
                    peak_index = np.argmax(wf_corrected)

                    if not (amplitude >= amplitude_threshold and peak_min_index <= peak_index <= peak_max_index):
                        continue

                    if single_pe_peak:
                        peaks, props = find_peaks(wf_corrected, height=min_peak_height, distance=min_peak_distance_samples)
                        if len(peaks) != 1:
                            continue  # scarta se ha più di un picco significativo

                    selected_waveforms.append(wf_corrected)
                    if get_edges:
                        rising_edges.append(waveform_analysis.find_rising_edge(wf_corrected))
                    total_kept += 1

                    if max_waveforms is not None and total_kept >= max_waveforms:
                        break
                if max_waveforms is not None and total_kept >= max_waveforms:
                    break
            else:
                print(f"⚠️ Skipping {key}: not a valid 2D waveform array.")

        print(f"[INFO] Analyzed {total_checked} waveforms, selected {total_kept} (ampl ≥ {amplitude_threshold} {unit}, peak in [{peak_min_index}, {peak_max_index}])")
        
        final_waveforms = []
        for wf in selected_waveforms:
            # wf = Filters.lowpass_filter(wf, fs, cutoff_hz=200e6, order=4)
            final_waveforms.append(wf)

        if get_edges:
            return final_waveforms, rising_edges
        else:
            return final_waveforms

###########################
# Amplitude analysis/distributions
###########################
    def amplitudes(waveforms, sampling=5000, output_filename=None):
        """
        Calcola le ampiezze di picco per una o più waveform.
        
        - waveforms: array 1D (una sola wf) o 2D (più wf)
        - sampling: frequenza di campionamento in MHz (solo per messaggi)
        - output_filename: se specificato, salva i risultati in un file .txt
        
        Restituisce un array di ampiezze di picco.
        """
        waveforms = np.asarray(waveforms)
        if waveforms.ndim == 1:
            waveforms = waveforms[np.newaxis, :]

        peak_vals = np.max(waveforms, axis=1)

        print(f"📊 Trovate {len(peak_vals)} waveform — Sampling: {sampling} MHz")
        print(f"Ampiezze medie: {np.mean(peak_vals):.2f} mV | Mediana: {np.median(peak_vals):.2f} mV")

        if output_filename:
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            np.savetxt(output_filename, peak_vals, fmt="%.6f")
            print(f"✅ Ampiezze salvate in {output_filename}")

        return peak_vals
        

###########################
# Persistence plot with waveform analysis
###########################
    def persistence(self, npz_path, unit='mV', amplitude_threshold=4, num_waveforms=10000, #maybe none?
                color='green', align="False", sampling=None, output_filename=None):
        wf_data = Utils.load_waveforms(npz_path)
        info = Utils.get_info(npz_path)

        selected_waveforms = []
        rising_edges = []
        fs = sampling
        baseline_range = (49, 973)
        selected_waveforms, rising_edges = waveform_analysis.analyze_waveforms(
            wf_data=wf_data,
            fs=fs,
            baseline_range=baseline_range,
            amplitude_threshold=amplitude_threshold,
            max_waveforms=num_waveforms,
            unit=unit,
            get_edges=True
        )

        if not selected_waveforms:
            print("⚠️  No waveforms passed the P2P threshold!")
            return
        # Initialize waveforms as the selected waveforms
        waveforms = selected_waveforms
        
        # === ALIGNMENT ===
        if align:
            waveforms = waveform_analysis.align_waveforms(waveforms, rising_edges, mode=align)
        
        # === Check if waveforms are available before plotting ===
        if not waveforms:
            print("⚠️ No valid waveforms to plot.")
            return

        # === PLOT ===
        plt.style.use('dark_background')
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('#000c1f')

        alpha_decay = np.linspace(0.05, 0.8, len(waveforms))
        Plotting.scope_plot(ax, waveforms, unit=unit, npz_path=npz_path, color=(1.0, 1.0, 0.0), align=align,sampling=sampling, plot_avg=True)

        plt.tight_layout()


        target_idx = 150 if align else None
        Plotting.add_inset_zoom(ax, waveforms, alpha_decay, default_target_idx=target_idx)

        # Add margin to the right for the legend
        fig.subplots_adjust(right=0.81)  # Adjust right margin to leave space for the legend

        # Estrai il valore del bias dal nome del file
        bias_value = info["bias"]
        # Legend text
        filename_info = f'Bias voltage: {bias_value} V\n' 
        filename_info += f"CAEN DGZ DT5742B\n"
        filename_info += f"Sampling rate: {sampling} GS\n"
        filename_info += f"Persistence: {"Auto" if align == "auto" else align}\n"
        # Add legend outside the plot
        fig.text(0.82, 0.5, filename_info, ha='left', va='center', fontsize=8, color='white', fontfamily='monospace')

        # === SAVE ===
        # Remove file extension and get relative path from 'data' --> BE CAREFUL to this!
        output_path = output_filename if output_filename else print(f"⚠️ No output filename provided, skipping save.")
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=1200)
        print(f"✅ Plot saved to: {output_path}")
        
        return waveforms



###########################
# Store median waveform 
###########################
    @staticmethod
    def store_avg(npz_path, unit='mV', amplitude_threshold=4, num_waveforms=10000, align=False, sampling=None):
        wf_data = Utils.load_waveforms(npz_path)
        baseline_range = (49, 973)

        selected_waveforms, rising_edges = waveform_analysis.analyze_waveforms(
            wf_data=wf_data,
            fs=sampling,
            baseline_range=baseline_range,
            amplitude_threshold=amplitude_threshold,
            max_waveforms=num_waveforms,
            unit=unit,
            get_edges=True
        )

        peak_values = []
        peak_indices = []

        for wf in selected_waveforms:
            sub_wave = wf[128:768]
            peak_idx = np.argmax(sub_wave) + 128
            peak_val = wf[peak_idx]
            peak_values.append(peak_val)
            peak_indices.append(peak_idx)

        peak_values = np.array(peak_values).reshape(-1, 1)
        kmeans = KMeans(n_clusters=2, random_state=0)
        kmeans.fit(peak_values)
        labels = kmeans.labels_
        cluster_centers = kmeans.cluster_centers_.flatten()
        lowest_cluster = np.argmin(cluster_centers)

        filtered_waveforms = [
            wf for wf, peak_val, label in zip(selected_waveforms, peak_values.flatten(), labels)
            if label == lowest_cluster
        ]

        if not filtered_waveforms:
            print("⚠️ Nessun waveform selezionato dopo clustering!")
            return [], []

        if align:
            filtered_waveforms = waveform_analysis.align_waveforms(filtered_waveforms, rising_edges, mode=align)

        peak_aligned_waveforms = [
            wf for wf in filtered_waveforms if 145 <= np.argmax(wf) <= 155
        ]

        single_peak_waveforms = []
        for wf in peak_aligned_waveforms:
            peaks, _ = find_peaks(wf, height=amplitude_threshold)
            main_peaks = peaks[(145 <= peaks) & (peaks <= 155)]
            if len(main_peaks) == 1:
                single_peak_waveforms.append(wf)

        if not single_peak_waveforms:
            print("⚠️ Nessuna waveform valida!")
            return [], []

        avg_wf = np.median(single_peak_waveforms, axis=0)
        return single_peak_waveforms, avg_wf



    ###########################
    # LATEST + DEBUG...
    ###########################
    
    def compute_tau(waveforms, sampling=1000, output_filename=None):
        """
        Calcola tau (tempo di decadimento a 1/e dal picco) su una o più waveform.
        - waveforms: array 1D (una sola wf) o 2D (più wf)
        - sampling: frequenza di campionamento in MHz
        - output_filename: se specificato, salva i risultati

        Restituisce una lista di tuple: (tau, tau_err, peak_val)
        """
        waveforms = np.asarray(waveforms)
        dt = 1e3 / sampling  # ns per sample
        print(f"Sampling rate: {sampling} MHz → dt = {dt:.3f} ns")

        results = []

        if waveforms.ndim == 1:
            waveforms = waveforms[np.newaxis, :]

        for idx, wf in enumerate(waveforms):
            peak_idx = np.argmax(wf)
            peak_val = wf[peak_idx]

            if peak_val < 4.0:
                print(f"[{idx+1}] ⏭️  Scartata: picco troppo basso ({peak_val:.2f} mV)")
                continue

            target_val = peak_val / np.e
            tau_ns = None
            tau_err_ns = None

            std = np.std(wf, ddof=1)
            err = std / np.sqrt(len(wf))

            for i in range(peak_idx + 1, len(wf)):
                if wf[i] <= target_val:
                    tau_ns = (i - peak_idx) * dt
                    slope = abs(wf[i] - wf[i - 1])
                    slope = max(slope, 1e-12)
                    tau_err_ns = dt * (err / slope)
                    break

            if tau_ns is not None:
                print(f"[{idx+1}] τ = {tau_ns:.2f} ± {tau_err_ns:.2f} ns (peak = {peak_val:.2f})")
                results.append((tau_ns, tau_err_ns, peak_val))
            else:
                print(f"[{idx+1}] ⚠️ Nessun crossing trovato per 1/e")

        if results and output_filename:
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            with open(output_filename, "w") as f:
                for r in results:
                    f.write(f"{r[0]:.6f}\t{r[1]:.6f}\t{r[2]:.6f}\n")
            print(f"✅ Risultati salvati in {output_filename}")

        return results


    
    @staticmethod
    def persistence_4_tau(npz_path, unit='mV', amplitude_threshold=4, num_waveforms=10000,
                color='green', align="False", sampling=None, output_filename=None):
        wf_data = Utils.load_waveforms(npz_path)
        info = Utils.get_info(npz_path)

        single_peak_waveforms, avg_wf = waveform_analysis.store_avg(
            npz_path=npz_path,
            unit=unit,
            amplitude_threshold=amplitude_threshold,
            num_waveforms=num_waveforms,
            align=align,
            sampling=sampling
        )

        if not single_peak_waveforms:
            print("⚠️ No valid waveforms to plot.")
            return

        waveform_analysis.compute_tau(single_peak_waveforms, avg_wf, sampling, output_filename)

        # Plot
        plt.style.use('dark_background')
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('#000c1f')

        alpha_decay = np.linspace(0.05, 0.8, len(single_peak_waveforms))
        Plotting.scope_plot(ax, single_peak_waveforms, unit=unit, npz_path=npz_path,
                            color=(1.0, 1.0, 0.0), align=align, sampling=sampling, plot_avg=True)

        plt.tight_layout()
        # Add margin to the right for the legend
        fig.subplots_adjust(right=0.81)

        bias_value = info["bias"]
        filename_info = f'Bias voltage: {bias_value} V\n' 
        filename_info += f"CAEN DGZ DT5742B\n"
        filename_info += f"Sampling rate: {sampling} GS\n"
        filename_info += f"Persistence: {"Auto" if align == "auto" else align}\n"
        fig.text(0.82, 0.5, filename_info, ha='left', va='center', fontsize=8, color='white', fontfamily='monospace')

        output_path = output_filename + ".png" if output_filename else info["output_path"]
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=1200)
        print(f"✅ Plot saved to: {output_path}")
        
        return single_peak_waveforms, avg_wf



















###########################
# COMPARE PLOTS - MAYBE SOMEWHERE ELSE?
###########################


    @staticmethod
    def compare_waveforms(npz_paths, unit='mV', amplitude_threshold=4.0, num_waveforms=10000, sampling=None):
        """
        Confronta le waveform da più file .npz: selezione sopra soglia, allineamento e normalizzazione.
        Nessun effetto persistenza: solo overlay comparativo.

        Parameters:
        - npz_paths: lista di percorsi ai file .npz
        - unit: 'mV' o 'ADC'
        - amplitude_threshold: soglia P2P minima
        - num_waveforms: numero di waveform da selezionare per file
        - sampling: frequenza di campionamento in GHz
        """
        try:
            plt.style.use('seaborn-darkgrid')
        except OSError:
            print("⚠️  Style 'seaborn-darkgrid' not found, falling back to 'default'")
            plt.style.use('default')

        all_waveforms = []
        labels = []

        for npz_path in npz_paths:
            wf_data = Utils.load_waveforms(npz_path)
            info = Utils.get_info(npz_path)

            selected_waveforms, rising_edges = waveform_analysis.analyze_waveforms(
                wf_data=wf_data,
                fs=sampling,
                baseline_range=(49, 973),
                amplitude_threshold=amplitude_threshold,
                max_waveforms=num_waveforms,
                unit=unit,
                get_edges=True
            )

            if not selected_waveforms:
                print(f"⚠️ No waveforms selected from: {npz_path}")
                continue

            aligned = waveform_analysis.align_waveforms(selected_waveforms, rising_edges, mode="auto")
            waveforms = [wf / np.max(wf) for wf in aligned]


            all_waveforms.append(waveforms)
            label = os.path.basename(npz_path).replace('.npz', '')
            labels.append(label)

        if not all_waveforms:
            print("❌ No valid waveform sets found.")
            return

        # Plot
        plt.figure(figsize=(10, 6))
        for i, wf_set in enumerate(all_waveforms):
            for wf in wf_set:
                plt.plot(wf, alpha=0.5, label=labels[i] if wf is wf_set[0] else None)

        plt.xlabel("Samples")
        plt.ylabel(f"Amplitude (normalized)")
        plt.title("Overlay of Aligned & Normalized Waveforms")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        outname = "_vs_".join(labels) + "_comparison.png"
        output_dir = os.path.join(Utils.find_data_directory(), "plots")
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, outname)
        plt.savefig(output_path, dpi=1200)
        print(f"✅ Comparison plot saved to: {output_path}")



###########################
# DECAY TIME WITH FIT
###########################
    @staticmethod
    def decay_time(npz_path, unit='mV', amplitude_threshold=4, num_waveforms=10000,
                color='green', align="False", single_pe_filter=False, sampling=None,
                fit_window=50, tail_sigma_cut=2.0, baseline_window=50):

        wf_data = Utils.load_waveforms(npz_path)
        info = Utils.get_info(npz_path)

        fs = sampling
        print(f"[INFO] Sampling rate: {fs} MHz")
        baseline_range = (49, 973)

        selected_waveforms, rising_edges = waveform_analysis.analyze_waveforms(
            wf_data=wf_data,
            fs=fs,
            baseline_range=baseline_range,
            amplitude_threshold=amplitude_threshold,
            max_waveforms=num_waveforms,
            unit=unit,
            get_edges=True
        )

        if not selected_waveforms:
            print("⚠️  No waveforms passed the P2P threshold!")
            return

        txt_output = os.path.join(
            os.path.dirname(npz_path),
            os.path.basename(npz_path).replace(".npz", "_decay_times.txt")
        )

        taus, tau_errors, amps, amp_errors, offsets, offset_errors, _ = \
            waveform_analysis.fit_decay_and_plot_tau_distribution(
                waveforms=selected_waveforms,
                sampling=fs,
                plot=False,
                single_pe_filter=single_pe_filter,
                filename=None,
                txt_output=txt_output,
                fit_window=fit_window,
                tail_sigma_cut=tail_sigma_cut,
                baseline_window=baseline_window
            )

        return taus, tau_errors, amps, amp_errors, offsets, offset_errors
    
    

    @staticmethod
    def fit_waveform_tau(
        waveforms,
        sampling=None,
        filename=None,
        txt_output=None,
        fit_window=50,
        baseline_window=50,     # ignorato per ora
        tail_sigma_cut=2.0      # ignorato per ora
    ):
        import scipy.signal

        if txt_output is None and filename:
            txt_output = os.path.splitext(filename)[0] + ".txt"

        taus, tau_errors, amps, amp_errors, offsets, offset_errors, fitted_waveforms = (
            [], [], [], [], [], [], []
        )

        sampling_rate_hz = sampling * 1e6
        time_step_ns = 1e9 / sampling_rate_hz

        print(f"\n🔍 Inizio fit per {len(waveforms)} waveform(s)...\n")

        for idx, wf in enumerate(waveforms):
            print(f"\n➡️ Waveform {idx+1}/{len(waveforms)}")

            # baseline a zero e soglia minima
            min_height = 4.0  # in mV, direttamente

            # Trovo tutti i picchi
            peaks, props = scipy.signal.find_peaks(wf, height=min_height)

            # Prendo solo quelli con rising-edge tra 145 e 155
            valid_peaks = [p for p in peaks if 145 <= p <= 155]
            if not valid_peaks:
                print(f"⚠️ Nessun picco tra 145 e 155 con altezza ≥{min_height} mV. Skip.")
                continue

            peak_idx = valid_peaks[0]
            peak_val = wf[peak_idx]

            if peak_val < 4.0:
                print(f"⏭️  Scartata: picco troppo basso ({peak_val:.2f} mV)")
                continue

            print(f"✅ Picco valido: idx {peak_idx}, ampiezza = {peak_val:.3f} mV")

            tail = wf[peak_idx : peak_idx + fit_window]
            t = np.arange(len(tail)) * time_step_ns

            if len(tail) < 5:
                print("⚠️ Troppo pochi punti per il fit. Skip.")
                continue

            A0 = tail[0]
            tau0 = 20.0
            C0 = 0.0  # baseline

            try:
                popt, pcov = curve_fit(
                    waveform_analysis.exp_decay, t, tail,
                    p0=[A0, tau0, C0],
                    bounds=([0, 1, -np.inf], [np.inf, 200, np.inf])
                )
                A_fit, tau_fit, C_fit = popt
                err = np.sqrt(np.diag(pcov))

                # ✅ Controllo: fit sospetto se A troppo bassa o tau al limite
                if A_fit <= 0.1 or tau_fit >= 199.0:
                    print(f"⚠️ Fit non attendibile: A = {A_fit:.2f}, τ = {tau_fit:.2f}. Skip.")
                    continue

                print(f"✅ Fit: τ = {tau_fit:.2f}±{err[1]:.2f} ns, A = {A_fit:.2f}, C = {C_fit:.2f}")

                taus.append(tau_fit)
                tau_errors.append(err[1])
                amps.append(A_fit)
                amp_errors.append(err[0])
                offsets.append(C_fit)
                offset_errors.append(err[2])
                fitted_waveforms.append(wf)

            except Exception as e:
                print(f"❌ Fit fallito su picco indice {peak_idx}: {e}")
                continue


        # Scrittura file
        if taus and txt_output:
            os.makedirs(os.path.dirname(txt_output), exist_ok=True)
            with open(txt_output, "w") as f:
                for t0, dt0, A, dA, C, dC in zip(
                    taus, tau_errors, amps, amp_errors, offsets, offset_errors
                ):
                    f.write(f"{t0:.6f}\t{dt0:.6f}\t{A:.6f}\t{dA:.6f}\t{C:.6f}\t{dC:.6f}\n")
            print(f"\n📄 Risultati salvati in '{txt_output}'")
        else:
            print("⚠️ Nessuna waveform fittata con successo.")

        return taus, tau_errors, amps, amp_errors, offsets, offset_errors, fitted_waveforms
 