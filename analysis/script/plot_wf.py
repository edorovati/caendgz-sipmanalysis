import numpy as np
import matplotlib.pyplot as plt

def visualizza_waveform(
    npz_path,
    sign=1,
    num_waveforms=5,
    x_start=0,
    x_end=1024,
    lowpass_filter=None,
    highpass_filter=None,
    notch_filter=None,
    adc_zero=4096,
    filter_obj=None
):
    def find_good_indices(wf, baseline_start, baseline_end,
                          window_size=40, std_threshold=1.5,
                          min_good_block=30, pre_margin=30, post_margin=150):
        region = wf[baseline_start:baseline_end]
        num_points = len(region)
        region_indices = np.arange(baseline_start, baseline_end)

        stds = np.array([
            np.std(region[i:i+window_size])
            for i in range(num_points - window_size)
        ])
        noisy_indices = np.where(stds > std_threshold)[0] + baseline_start

        exclude_mask = np.zeros_like(wf, dtype=bool)
        for idx in noisy_indices:
            start = max(0, idx - pre_margin)
            end = min(len(wf), idx + post_margin)
            exclude_mask[start:end] = True

        usable_indices = region_indices[~exclude_mask[baseline_start:baseline_end]]

        if len(usable_indices) > 0:
            diffs = np.diff(usable_indices)
            block_edges = np.where(diffs > 1)[0]
            blocks = np.split(usable_indices, block_edges + 1)
            good_blocks = [b for b in blocks if len(b) >= min_good_block]
            final_indices = np.concatenate(good_blocks) if good_blocks else np.array([])
        else:
            final_indices = np.array([])
        return final_indices

    def calculate_baseline(wf, indices):
        return np.median(wf[indices]) if len(indices) > 0 else np.median(wf)

    # Carica dati
    data = np.load(npz_path)
    if "waveforms" not in data:
        raise ValueError(f"'waveforms' non trovato in {npz_path}")
    waveforms = data["waveforms"].squeeze()

    waveforms = waveforms[:num_waveforms]

    all_corr = []

    for i, wf in enumerate(waveforms):
        mv = ((wf - 2048) / adc_zero) * 1000
        good_indices = find_good_indices(mv, 49, 973)
        baseline = calculate_baseline(mv, good_indices)
        #corr = mv - baseline
        mv = mv[x_start:x_end]
        

        # Applica filtri se definiti
        '''
        if filter_obj:
            if lowpass_filter:
                corr = filter_obj.lowpass(mv, lowpass_filter)
            if highpass_filter:
                corr = filter_obj.highpass(mv, highpass_filter)
            if notch_filter:
                corr = filter_obj.notch(mv, notch_filter)
        '''
        all_corr.append(mv)

        dt = 1e9 / 750e6  # 1.333... ns per campione

    # Plot
    for i, corr in enumerate(all_corr):
        time_axis = np.arange(len(corr)) * dt
        plt.plot(time_axis, corr, label=f"WF {i+1}")
    plt.xlabel("Tempo (ns)")
    plt.ylabel("Ampiezza (mV)")
    plt.title("Waveform Sovrapposti")
    plt.grid(True)
    
    plt.show()



# === Qui inizia l'esempio di chiamata ===

file_path = "TEMP-2.5Gs-run3-grosso.npz"  # Metti qui il percorso del tuo file .npz

visualizza_waveform(
    file_path,
    num_waveforms=2000,
    x_start=0,
    x_end=1024
)
