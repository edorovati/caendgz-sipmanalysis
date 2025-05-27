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







class waveform_analysis:
    def __init__(self, data_dir=None):
        """
        Initialize the waveform analyzer with a default data directory.
        """
        if data_dir is None:
            self.data_dir = Utils.find_data_directory()
        else:
            self.data_dir = os.path.abspath(data_dir)

        print("Using data directory:", self.data_dir)

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

    @staticmethod
    def plot_first_waveform(npz_path, data, unit='ADC'):
        wf_data = Utils.load_waveforms(npz_path) #called for nothing
        info = Utils.get_info(npz_path)
        for key in data.files:
            wf_array = data[key]
            if wf_array.ndim == 2 and wf_array.shape[0] >= 10:
                for i in range(10):
                    wf = wf_array[i+10]
                    if unit == 'mV':
                        wf = Utils.adc_to_mv(wf)
                        baseline, _ = waveform_analysis.calculate_baseline_with_mask(wf, 49, 973)
                        wf -= baseline
                        fs = info["sampling_rate_ghz"]*1e9
                        wf = Filters.lowpass_filter(wf, fs, cutoff_hz=200e6, order=4)
                    plt.plot(wf, alpha=0.4, label=f"{key} #{i}" if i == 0 else "")
            else:
                print(f"⚠️ Skipping {key}: not enough waveforms or not a 2D array.")

        plt.title(f"First 100 waveforms per channel ({unit})")
        plt.xlabel("Sample index")
        plt.ylabel(f"Amplitude [{unit}]")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
        
    @staticmethod
    def find_rising_edge(wf, threshold_factor=0.63):
        peak = np.max(wf)
        threshold = peak * threshold_factor
        for i in range(128, 768):
            if wf[i] > threshold and wf[i - 1] <= threshold:
                return i
        return 0

    @staticmethod
    def find_polarity(waveform):
        """
        Determines the dominant signal polarity of a baseline-subtracted waveform.
        Returns 'positive' or 'negative'.
        """
        max_val = np.max(waveform)
        min_val = np.min(waveform)

        if abs(min_val) > abs(max_val):
            return 'negative'
        else:
            return 'positive'
       
    @staticmethod 
    def align_waveforms(waveforms, rising_edges, mode="auto"):
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




    def persistence(self, npz_path, unit='mV', amplitude_threshold=4, num_waveforms=100, color='green', align="False", sampling = None):
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

        # TO DO! fix alignment 
        # Align waveforms if required --> to be fixed and tested in few cases
        if align:
            waveforms = waveform_analysis.align_waveforms(waveforms, rising_edges, mode=align)

        # === Check if waveforms are available before plotting ===
        if not waveforms:
            print("⚠️ No valid waveforms to plot.")
            return

        # Plot
        plt.style.use('dark_background')
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('#000c1f')

        alpha_decay = np.linspace(0.05, 0.8, len(waveforms))
        Plotting.scope_plot(ax, waveforms, unit=unit, npz_path=npz_path, color=(1.0, 1.0, 0.0), align=align,sampling=sampling)

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

        # Save
        # Remove file extension and get relative path from 'data'
        output_path = info["output_path"]
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=1200)
        print(f"✅ Plot saved to: {output_path}")



    

        # # Plot delle waveform con tau ≈ 5 ns
        # if fits_for_plot:
   
        #         y = wf[peak_idx:]
        #         t = np.arange(len(y)) * time_step_ns
        #         y_fit = waveform_analysis.exp_decay(t, *popt)

        #         plt.figure(figsize=(8, 4))
        #         plt.plot(np.arange(len(wf)), wf, label="Waveform", color='steelblue')
        #         plt.plot(np.arange(peak_idx, peak_idx + len(y_fit)), y_fit, label=f"Fit (τ = {popt[1]:.2f} ns)", color='orange')
        #         plt.axvline(peak_idx, linestyle='--', color='gray', label="Peak index")
        #         plt.xlabel("Sample Index")
        #         plt.ylabel("Amplitude [mV]")
        #         plt.title("Waveform con τ = 5 ns")
        #         plt.legend()
        #         plt.grid(True)
        #         plt.tight_layout()
        #         plt.show()

        # return taus, fitted_waveforms











    #     # THIS HAS TO BE COMPLETELY REWRITTEN IN VIEW OF THE NEW CHANGES --> SEE PERSISTENCE
    #     #waveform_analysis
    #     def recdec_time(self, npz_path, unit='mV', p2p_threshold=4.0, num_waveforms=100, align="auto", custom_filename=None):
    #         info = waveform_analysis.get_info(npz_path)
    #         wf_data = waveform_analysis.load_waveforms(npz_path)
    #         info = Utils.get_info(npz_path)
    #         fs = info["sampling_rate_ghz"]*1e9

    #         baseline_start, baseline_end = 49, 973
    #         selected_waveforms = []
    #         rising_edges = []
    #         peak_heights = []
    #         valid_waveforms_count = 0

    #         for key in wf_data.files:
    #             wf_array = wf_data[key]
    #             if wf_array.ndim == 2 and wf_array.shape[0] > 0:
    #                 for i in range(wf_array.shape[0]):
    #                     wf_adc = wf_array[i]
    #                     wf = wf_adc if unit != 'mV' else wf_adc

    #                     wf = Utils.adc_to_mv(wf_adc) if unit == 'mV' else wf_adc  # conversion enabled
    #                     baseline, _ = waveform_analysis.calculate_baseline_with_mask(wf, baseline_start, baseline_end)
    #                     wf -= baseline
    #                     fs = info["sampling_rate_ghz"]*1e9
    #                     wf = Filters.lowpass_filter(wf, fs, cutoff_hz=200e6, order=4)


    #                     p2p = np.ptp(wf[baseline_start:baseline_end])
    #                     if p2p >= p2p_threshold:
    #                         selected_waveforms.append(wf)
    #                         rising_edges.append(waveform_analysis.find_rising_edge(wf))
    #                         peak_heights.append(np.max(wf))
    #                         valid_waveforms_count += 1

    #                     if valid_waveforms_count >= num_waveforms:
    #                         break
    #                 if valid_waveforms_count >= num_waveforms:
    #                     break

    #         if valid_waveforms_count == 0:
    #             print("⚠️  No waveforms passed the P2P threshold!")
    #             return

    #         if align:
    #             selected_waveforms = waveform_analysis.align_waveforms(selected_waveforms, rising_edges, mode=align)

    #         if not selected_waveforms:
    #             print("⚠️ No valid waveforms to plot.")
    #             return

    #         peak_heights = np.array(peak_heights)
    #         span = 3.0  # mV
    #         n_sigma = 3 # let's say like this but it may change in the future

    #         # Calcola il rumore medio (RMS) nella baseline delle waveform valide
    #         baseline_rms_values = [np.std(wf[baseline_start:baseline_end]) for wf in selected_waveforms]
    #         avg_rms = np.mean(baseline_rms_values)
    #         min_peak_threshold = avg_rms * n_sigma

    #         # Trova il gruppo di waveform con il picco minimo sopra soglia rumore
    #         valid_peak_indices = [i for i, pk in enumerate(peak_heights) if pk > min_peak_threshold]
    #         if not valid_peak_indices:
    #             print("⚠️ I picco sopra la soglia rumore. Riduci n_sigma o controlla i dati.")
    #             return

    #         valid_peaks = np.array([peak_heights[i] for i in valid_peak_indices])
    #         min_peak = np.min(valid_peaks)

    #         # Seleziona waveform attorno al picco più basso valido
    #         close_indices = [i for i in valid_peak_indices if abs(peak_heights[i] - min_peak) < span]

    #         used_waveforms = [selected_waveforms[i] for i in close_indices]
    #         other_waveforms = [selected_waveforms[i] for i in range(len(selected_waveforms)) if i not in close_indices]

    #         avg_waveform = np.mean(used_waveforms, axis=0)

    #         # Plotting
    #         plt.style.use('dark_background')
    #         fig, ax = plt.subplots(figsize=(10, 6))
    #         fig.patch.set_facecolor('black')
    #         ax.set_facecolor('#000c1f')
    #         ax.grid(color='white', alpha=0.2)

    #         alpha_decay = np.linspace(0.05, 0.8, len(selected_waveforms))
    #         x_time = Utils.cell_to_seconds(np.arange(len(wf)), npz_path) * 1e9
    #         for idx, wf in enumerate(other_waveforms):
    #             ax.plot(x_time, wf, color='yellow', alpha=0.3, linewidth=1)
    #         for idx, wf in enumerate(used_waveforms):
    #             ax.plot(x_time, wf, color='red', alpha=0.5, linewidth=1)
    #         ax.set_xlabel("Time [ns]")

    #         ax.plot(x_time, avg_waveform, color='white', linewidth=2.0, label='Average SPE')

    #         ax.set_title(f'SPE Overlay: {len(used_waveforms)} Selected (Red) for Avg + Mean (White) [{unit}]')
    #         ax.set_ylabel(f"Amplitude [{unit}]")
    #         plt.tight_layout()

    #         target_idx = 150 if align else None
    #         # to be fixed with the right time conversion
    #         Plotting.add_inset_zoom(ax, used_waveforms, np.linspace(0.3, 0.8, len(used_waveforms)), default_target_idx=target_idx)
    #         base_filename = custom_filename if custom_filename else os.path.basename(npz_path).replace(".npz", "")
    #         bias_match = re.search(r'bias([\d\.]+)', npz_path)
    #         bias_value = bias_match.group(1) if bias_match else 'N/A'   
    #         # Prepare the legend text
    #         # filename_info = f"File: {os.path.basename(npz_path)}\n"
    #         # filename_info += f"Waveforms: {valid_waveforms_count}\n"
    #         # filename_info += f"Unit: {unit}\n"
    #         # filename_info += f"P2P Threshold: {p2p_threshold} {unit}\n"
    #         filename_info = f'Bias voltage: {bias_value} V\n' 
    #         filename_info += f"CAEN DGZ DT5742B\n"
    #         filename_info += f"Sampling rate: {info["rel_path"] / 1e9} GS\n"
    #         filename_info += f"Persistence: {"Auto" if align == "auto" else align}\n"

    #         os.makedirs("./analysis", exist_ok=True)
            
    #         output_path = os.path.join("./analysis", base_filename + "_SPE_overlay.png")
    #         plt.savefig(output_path, dpi=1200)
    #         print(f"✅ Plot with average waveform saved to: {output_path}")

    # #######################
    @staticmethod
    def exp_decay(t, A, tau, C):
        return A * np.exp(-t / tau) + C



    
    @staticmethod
    def fit_decay_and_plot_tau_distribution(waveforms, sampling=None, plot=True, single_pe_filter=False, filename=None):
        """
        Fit esponenziale del decadimento di ciascuna waveform dal picco in poi,
        raccoglie i tau e opzionalmente mostra la loro distribuzione.

        Parameters:
        - waveforms: lista di array 1D, già baseline-corrected
        - sampling: frequenza di campionamento in MHz
        - plot: se True, mostra la distribuzione dei tau
        - single_pe_filter: se True, accetta solo waveform con un singolo picco ≥ 4 mV

        Returns:
        - lista di tau (in ns) per i fit riusciti
        - lista delle waveform corrispondenti
        """
        from scipy.signal import find_peaks
        from scipy.optimize import curve_fit
        import numpy as np
        import matplotlib.pyplot as plt

        taus = []
        fitted_waveforms = []

        sampling_rate_hz = sampling * 1e6
        time_step_ns = 1e9 / sampling_rate_hz

        for wf in waveforms:
            if single_pe_filter:
                peaks, _ = find_peaks(wf, height=4.0)
                if len(peaks) != 1:
                    continue

            peak_idx = np.argmax(wf)
            y = wf[peak_idx:]
            t = np.arange(len(y)) * time_step_ns

            A0 = y[0]
            tau0 = 50.0  # ns, stima iniziale fissa
            C0 = np.mean(y[-min(10, len(y)):])  # media ultimi N punti

            try:
                popt, _ = curve_fit(
                    waveform_analysis.exp_decay, t, y,
                    p0=[A0, tau0, C0],
                    bounds=([0, 1, -np.inf], [np.inf, 200, np.inf])
                )
                tau_fit = float(popt[1])
                taus.append(tau_fit)
                fitted_waveforms.append(wf)
            except Exception:
                continue
                
        if plot and taus:
            plt.figure(figsize=(7, 4))
            plt.hist(taus, bins=200, color='steelblue', edgecolor='black')
            plt.xlabel("Tau [ns]")
            plt.ylabel("Counts")
            plt.title("Distribuzione costante di decadimento (τ)")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(filename, dpi=1200)
            print(f"✅ Histogram saved to: {filename}")

        return taus






    @staticmethod
    def decay_time(npz_path, unit='mV', amplitude_threshold=4, num_waveforms=100,
                color='green', align="False", single_pe_filter=False, sampling = None):
        wf_data = Utils.load_waveforms(npz_path)
        info = Utils.get_info(npz_path)

        selected_waveforms = []
        rising_edges = []
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

        waveforms = selected_waveforms

        # # === Plot ===
        # plt.style.use('dark_background')
        # fig, ax = plt.subplots(figsize=(10, 6))
        # fig.patch.set_facecolor('black')
        # ax.set_facecolor('#000c1f')

        # alpha_decay = np.linspace(0.05, 0.8, len(waveforms))
        # Plotting.scope_plot(ax, waveforms, unit=unit, npz_path=npz_path, color=(1.0, 1.0, 0.0), align=align)
        # plt.tight_layout()

        # target_idx = 150 if align else None
        # Plotting.add_inset_zoom(ax, waveforms, alpha_decay, default_target_idx=target_idx)

        # # Legend
        # fig.subplots_adjust(right=0.81)
        # bias_value = info["bias"]
        # filename_info = f'Bias voltage: {bias_value} V\n' 
        # filename_info += f"CAEN DGZ DT5742B\n"
        # filename_info += f"Sampling rate: {info["sampling_rate_hz"]/1e9:.2f} GS\n"
        # filename_info += f"Persistence: {"Auto" if align == "auto" else align}\n"
        # fig.text(0.82, 0.5, filename_info, ha='left', va='center', fontsize=8, color='white', fontfamily='monospace')

        # Save figure
        # output_path = info["output_path"]
        # os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # plt.savefig(output_path, dpi=1200)
        # print(f"✅ Plot saved to: {output_path}")

        # === Tau fit ===
        decay_times  = waveform_analysis.fit_decay_and_plot_tau_distribution(
            waveforms=waveforms,
            sampling=sampling,
            plot=True,
            single_pe_filter=single_pe_filter, filename=info["output_path"]+ "_decay_time_distribution.png"
        )
        if decay_times:
            print(f"✅ Decay time fit successful: {len(decay_times)} waveforms fitted")
            print(f"Average tau: {np.mean(decay_times):.2f} ns, Std: {np.std(decay_times):.2f} ns")
        
        # #let's plot the "involved" waveforms, given by the secodn "return" of the fit function
        # if len(decay_times_wf) > 0:
        #     plt.figure(figsize=(10, 6))
        #     plt.plot(decay_times_wf, color='blue', alpha=0.5)
        #     plt.title("Decay Time Waveforms")
        #     plt.xlabel("Sample Index")
        #     plt.ylabel(f"Amplitude [{unit}]")
        #     plt.grid(True)
        #     plt.tight_layout()
        #     plt.show()
        # else:
        #     print("⚠️ No waveforms were fitted successfully.")