import numpy as np
import argparse
from utils.filters import Filter
from utils.plotter import Plotter
from utils.scanthr import ScanThreshold
from utils.utils import Utils
from colorama import Fore, Style  # For colored terminal outputs

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
#    python3 transitions.py --npz waveforms_5GS.npz --waveform mV --sign -1 --num_waveforms 10 --lowpass 200e6
#
# 2. Display ADC or mV scaled waveforms:
#    python3 transitions.py --npz waveforms_5GS.npz --waveform ADC --sign 1 --num_waveforms 5
#
# 3. Scan thresholds between 0-15 with a lowpass filter:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1 --lowpass 200e6
#
# 4. Use a fixed threshold of 5 and apply a lowpass filter:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6
#
# 5. Same as above, but restrict processing to a subrange of the waveform:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6 --x_start 100 --x_end 900
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

args = parser.parse_args()

########################## === Constants === ##########################
ADC_ZERO = 2048              # Zero level for ADC values
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
if args.waveform:
    for i in range(num_to_plot):
        wf = waveforms[i]
        corr = wf - np.median(wf)
        mv = (corr / ADC_ZERO) * 1000 if args.waveform == "mV" else corr
        mv = mv[args.x_start:args.x_end]
        mv = -mv if args.sign == -1 else mv
        if args.lowpass:
            mv = filter_obj.lowpass(mv, args.lowpass)
        if args.highpass:
            mv = filter_obj.highpass(mv, args.highpass)
        if args.notch:
            mv = filter_obj.notch(mv, args.notch)
        plotter_obj.plot_waveform(mv, f"Waveform Displayed {i+1}", "Amplitude", f"Waveform {i+1}", color="blue")
    exit()
###############################################################################

########################## === Threshold Scan Mode === ##########################
if args.scanthr and args.range:
    print(Fore.CYAN + "\nPre-processing the waveforms for scanning..." + Style.RESET_ALL)
    preprocessed_waveforms = []
    for wf in waveforms:
        corr = wf - np.median(wf)
        mv = (corr / ADC_ZERO) * 1000
        mv = mv[args.x_start:args.x_end]
        mv = -mv if args.sign == -1 else mv
        if args.lowpass:
            mv = filter_obj.lowpass(mv, args.lowpass)
        if args.highpass:
            mv = filter_obj.highpass(mv, args.highpass)
        if args.notch:
            mv = filter_obj.notch(mv, args.notch)
        preprocessed_waveforms.append(mv)
    print(Fore.GREEN + "Pre-processing completed." + Style.RESET_ALL)

    if "singlethr" in args.range:
        threshold_value = float(args.range.split()[-1])
        output_file = Utils.generate_output_filename(args.npz)
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
        range_min, range_max = map(float, args.range.split('-'))
        output_file = Utils.generate_output_filename(args.npz)
        scan_obj.scan_thresholds(
            preprocessed_waveforms,
            args.sign,
            range_min,
            range_max,
            0.50,
            output_file
        )
        print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)
#################################################################################
