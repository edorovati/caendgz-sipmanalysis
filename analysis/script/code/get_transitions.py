import sys
import os
import glob
import re

base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../python'))
sys.path.append(base_path)

from filters import Filter
from plotter import Plotter
from scanthr import ScanThreshold
from utils import Utils
import numpy as np
import argparse
from colorama import Fore, Style
import uproot
import awkward as ak
'''
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
'''
def find_good_indices(
    wf, baseline_start, baseline_end,
    window_size=40, std_threshold=1.5,
    min_good_block=30, pre_margin=30, post_margin=150
):
    """
    Restituisce gli indici buoni escludendo zone rumorose.
    """
    region = wf[baseline_start:baseline_end]
    num_points = len(region)
    region_indices = np.arange(baseline_start, baseline_end)

    stds = np.array([
        np.std(region[i:i+window_size])
        for i in range(num_points - window_size)
    ])

    noisy = stds > std_threshold
    noisy_indices = np.where(noisy)[0] + baseline_start

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
    """
    Calcola la baseline come mediana dei valori sugli indici forniti.
    """
    if len(indices) > 0:
        return np.median(wf[indices])
    else:
        return np.median(wf)


def calculate_rms(wf, indices):
    """
    Calcola l'RMS nella regione specificata rispetto alla mediana locale,
    non rispetto alla baseline globale.
    """
    if len(indices) == 0:
        return np.float32(0)
    
    local_region = wf[indices]
    local_median = np.median(local_region)
    residuals = local_region - local_median
    rms = np.sqrt(np.mean(residuals ** 2)).astype(np.float32)
    return rms


##################### Waveform Analysis Script ###############################
# This script analyzes waveform data with optional thresholding and filtering.
#
# Features:
# 1. Visualize waveforms with selectable scaling (mV or ADC).
# 2. Apply digital filters (lowpass, highpass, notch).
# 3. Scan waveforms using threshold ranges or fixed thresholds.
# 4. Output hit rate (Hz) as a function of threshold.
#
# Command-line arguments:
# --npz              : Name of the .npz file containing waveform data (default: waveforms_5GS.npz).
# --waveform         : Select scale for waveform display ("mV" or "ADC").
# --sign             : Polarity of transition to detect (+1 or -1).
# --num_waveforms    : Number of waveforms to display (default: 5).
# --lowpass          : Apply a lowpass filter with the given cutoff frequency (Hz).
# --highpass         : Apply a highpass filter with the given cutoff frequency (Hz).
# --notch            : Apply a notch filter at the specified frequency (Hz).
# --scanthr          : Enable threshold scanning mode.
# --range            : Set threshold scanning range. Format: "min-max" or "singlethr value".
# --x_start          : Starting index for waveform slicing (default: 0).
# --x_end            : Ending index for waveform slicing (default: 1024).
#
# Example usages:
# 1. Display waveforms in mV with a lowpass filter:
#    python3 get_transitions.py --npz waveforms_5GS.npz --waveform mV --sign -1 --num_waveforms 10 --lowpass 200e6
#
# 2. Display ADC or mV scaled waveforms:
#    python3 get_transitions.py --npz waveforms_5GS.npz --waveform ADC --sign 1 --num_waveforms 5
#
# 3. Scan thresholds between 0-15 with a lowpass filter:
#    python3 get_transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1 --lowpass 200e6
#
# 4. Use a fixed threshold of 5 and apply a lowpass filter:
#    python3 get_transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6
#
# 5. Same as above, but restrict processing to a subrange of the waveform:
#    python3 get_transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6 --x_start 100 --x_end 900
###############################################################################

########################## === Argument Parsing === ##########################
parser = argparse.ArgumentParser(description="Waveform analysis with thresholds and filters")
parser.add_argument("--waveform", choices=["mV", "ADC"], help="Scale to display waveforms in (mV or ADC)")
parser.add_argument("--sign", type=int, choices=[-1, 1], default=-1, help="Polarity of transition to detect (+1 or -1)")
parser.add_argument("--npz", type=str, default="waveforms_5GS.npz", help="Path to the .npz file containing waveform data")
parser.add_argument("--num_waveforms", type=int, default=5, help="Number of waveforms to display")
parser.add_argument("--lowpass", type=float, help="Lowpass filter cutoff frequency (Hz)")
parser.add_argument("--highpass", type=float, help="Highpass filter cutoff frequency (Hz)")
parser.add_argument("--notch", type=float, help="Notch filter frequency (Hz)")
parser.add_argument("--scanthr", action="store_true", help="Enable threshold scan mode")
parser.add_argument("--range", type=str, help="Threshold range: 'min-max' or 'singlethr value'")
parser.add_argument("--x_start", type=int, default=0, help="Start index of waveform region to process")
parser.add_argument("--x_end", type=int, default=1024, help="End index of waveform region to process")
parser.add_argument("--to_root", action="store_true", help="Esporta tutte le ampiezze in un file ROOT con TTree")

args = parser.parse_args()

########################## === Constants === ##########################
ADC_ZERO = 4096              # Zero level for ADC values
SAMPLING_RATE = 750e6        # Sampling rate in Hz
#######################################################################

########################## === Load Data === ##########################
print("\n" + "="*80)
print(Fore.CYAN + "Data Loading Section" + Style.RESET_ALL)
print("="*80 + "\n")
data = np.load(args.npz)
if "waveforms" not in data:
    print(Fore.RED + f"Error: {args.npz} file not found or missing 'waveforms'." + Style.RESET_ALL)
    exit()
print(Fore.GREEN + "File loaded successfully!" + Style.RESET_ALL)
waveforms = data["waveforms"].squeeze()
##########################################################################

filter_obj = Filter(SAMPLING_RATE)
plotter_obj = Plotter()
scan_obj = ScanThreshold()
num_to_plot = min(args.num_waveforms, len(waveforms))

########################## === Waveform Display === ##########################
'''if args.waveform:
    for i in range(num_to_plot):
        wf = waveforms[i]
        mv = ((wf - 2048) / ADC_ZERO) * 1000  # Primo passaggio: conversione
        baseline, _ = calculate_baseline_with_mask(mv, args.x_start,args.x_end)  # usa 0-300 o valori personalizzati
        corr = mv - baseline
        corr = corr[args.x_start:args.x_end]
        corr = -corr if args.sign == -1 else corr
        if args.lowpass:
            corr = filter_obj.lowpass(corr, args.lowpass)
        if args.highpass:
            corr = filter_obj.highpass(corr, args.highpass)
        if args.notch:
            corr = filter_obj.notch(corr, args.notch)
        plotter_obj.plot_waveform(corr, f"Waveform Displayed {i+1}", "Amplitude", f"Waveform {i+1}", color="blue")
    exit()'''
    '''
if args.waveform:
    all_waveforms = []
    labels = []

    for i in range(num_to_plot):
        wf = waveforms[i]
        mv = ((wf - 2048) / ADC_ZERO) * 1000
        good_indices = find_good_indices(mv, 49, 973)
        baseline = calculate_baseline(mv, good_indices)
        corr = mv - baseline
        corr = corr[args.x_start:args.x_end]
        corr = -corr if args.sign == -1 else corr

        # Applica i filtri
        if args.lowpass:
            corr = filter_obj.lowpass(corr, args.lowpass)
        if args.highpass:
            corr = filter_obj.highpass(corr, args.highpass)
        if args.notch:
            corr = filter_obj.notch(corr, args.notch)

        all_waveforms.append(corr)
        labels.append(f"WF {i+1}")

       

    # Plot all waveforms together
    plotter_obj.plot_multiple_waveforms(
        all_waveforms,
        labels,
        title="Overlayed Waveforms",
        xlabel="Sample",
        ylabel="Amplitude (mV)"
    )
    exit()
'''

###############################################################################

########################## === Threshold Scan Mode === ##########################
if args.scanthr and args.range:
    print(Fore.CYAN + "\nPre-processing the waveforms for scanning..." + Style.RESET_ALL)
    preprocessed_waveforms = []
    waveforms = waveforms[:args.num_waveforms]  # limita il numero di waveform
    for wf in waveforms:
        mv = ((wf - 2048) / ADC_ZERO) * 1000
        good_indices = find_good_indices(mv, 49, 973)
        baseline = calculate_baseline(mv, good_indices)
        corr = mv - baseline
        corr = corr[args.x_start:args.x_end]
        corr = -corr if args.sign == -1 else corr
        if args.lowpass:
            corr = filter_obj.lowpass(corr, args.lowpass)
        if args.highpass:
            corr = filter_obj.highpass(corr, args.highpass)
        if args.notch:
            corr = filter_obj.notch(corr, args.notch)
        preprocessed_waveforms.append(corr)
    print(Fore.GREEN + "Pre-processing completed." + Style.RESET_ALL)

    if "singlethr" in args.range:
        threshold_value = float(args.range.split()[-1])
        #output_file = Utils.generate_output_filename(args.npz)
        output_file = args.npz
        scan_obj.scan_thresholds(
            preprocessed_waveforms,
            args.sign,
            threshold_value,
            threshold_value,
            0.50,
            output_file
        )
        print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)

    elif '-' in args.range:
        #range_min, range_max = map(float, args.range.split('-'))
        #output_file = Utils.generate_output_filename(args.npz)
        output_file = args.npz
        scan_obj.scan_thresholds(
            preprocessed_waveforms,
            args.sign,
            range_min,
            range_max,
            0.5,
            output_file
        )
        print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)
#################################################################################
'''
import matplotlib.pyplot as plt
from scipy.stats import norm
import uproot
import awkward as ak
import numpy as np

if args.to_root:
    input_path = args.npz
    input_dir = os.path.dirname(input_path)
    
    # Prende tutti i file .npz nella directory
    npz_files = sorted(glob.glob(os.path.join(input_dir, "*.npz")))

    # Crea nome cartella output (es. hamamatsu/LASER_ON/ -> HAMA-LASER_ON)
    folder_parts = os.path.normpath(input_dir).split(os.sep)
    output_dir = "-".join(p.upper() for p in folder_parts)
    os.makedirs(output_dir, exist_ok=True)

    print(Fore.CYAN + f"\nEsportazione in ROOT in corso per {len(npz_files)} file..." + Style.RESET_ALL)

    for npz_file in npz_files:
        waveforms = np.load(npz_file)["waveforms"]
        waveforms = waveforms[:args.num_waveforms]  # limita il numero di waveform

        amplitude_max = []
        charge_list = []
        noise_rms_list = []

        for wf in waveforms:
            mv = ((wf - 2048) / ADC_ZERO) * 1000

            good_indices = find_good_indices(mv, 49, 973)
            baseline = calculate_baseline(mv, good_indices)
            corr = mv - baseline
            rms = calculate_rms(corr, good_indices)

            amp_slice = corr[550:650]
            amp_max = np.max(amp_slice).astype(np.float32)
            amplitude_max.append(amp_max)

            dt = (1 / SAMPLING_RATE)
            charge = (np.sum(amp_slice) * dt).astype(np.float32)
            charge_list.append(charge)
            noise_rms_list.append(rms)

        # Estrae nome ROOT file (es. bias42.0V.root)
        basename = os.path.basename(npz_file)
        match = re.search(r"bias(\d+\.\d+)V", basename)
        if match:
            voltage = match.group(1).replace(".", ",")  # usa la virgola come richiesto
            root_filename = f"Bias{voltage}V.root"
        else:
            root_filename = "output.root"  # fallback di sicurezza

        root_path = os.path.join(output_dir, root_filename)

        # Salva su ROOT
        tree_info = {
            "amplitude": ak.Array(amplitude_max),
            "charge": ak.Array(charge_list),
            "noise_rms": ak.Array(noise_rms_list),
        }

        with uproot.recreate(root_path) as root_file:
            root_file["Analysis"] = tree_info

        print(Fore.GREEN + f"Creato: {root_path}" + Style.RESET_ALL)

    print(Fore.YELLOW + "\nElaborazione terminata." + Style.RESET_ALL)
'''
