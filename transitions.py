import numpy as np
import argparse
from utils.filters import Filter
from utils.plotter import Plotter
from utils.scanthr import ScanThreshold
from colorama import Fore, Style  # Importing colorama for colored outputs

##################### Basic Analysis of Waveforms ###############################
# This script is designed to analyze waveform data with thresholding and filtering options.
# The script allows for the following functionalities:
# 1. Visualizing waveforms with different scaling (mV or ADC).
# 2. Applying filters (lowpass, highpass, notch) to the waveforms.
# 3. Scanning the waveforms with thresholds and saving the results => Hit rate (Hz) vs Threshold.
# 4. Extract the Hit rate (Hz) at fixed threshold.
#
# The script requires the use of the following arguments:
#
# --npz              : Specify the name of the .npz file to load (default: waveforms_5GS.npz).
# --waveform         : Choose the type of scale for the waveform display (mV or ADC).
# --sign             : Specify the direction of the front to search for (+1 or -1).
# --num_waveforms    : Set the number of waveforms to display (default: 5).
# --lowpass          : Set the cutoff frequency for the lowpass filter in Hz.
# --highpass         : Set the cutoff frequency for the highpass filter in Hz.
# --notch            : Set the frequency for the notch filter in Hz.
# --scanthr          : Enable threshold scanning mode.
# --range            : Define the range for scanning: 'min-max' or 'singlethr value'.
#
# Example usage:
# 1. To visualize waveforms with scaling in mV and apply a lowpass filter:
#    python3 transitions.py --npz waveforms_5GS.npz --waveform mV --sign -1 --num_waveforms 10 --lowpass 200e6
#
# 2. To visualize waveforms with either mV or ADC scaling and change the sign:
#    python3 transitions.py --npz waveforms_5GS.npz --waveform mV/ADC --sign -1/+1 --num_waveforms 10
#
# 3. To scan waveforms within a threshold range (0-15) and apply a lowpass filter:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1 --lowpass 200e6
#
# 4. To scan waveforms within a threshold range (0-15) without any filter:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1
#
# 5. To scan waveforms with a single threshold value (5) and apply a lowpass filter:
#    python3 transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6
#
# These commands give users flexibility to analyze waveforms based on different filter settings and threshold scans.
###############################################################

########################## === Argument Parsing === ##########################
parser = argparse.ArgumentParser(description="Waveform analysis with thresholds and filters")
parser.add_argument("--waveform", choices=["mV", "ADC"], help="Type of scale for displaying waveforms (mV or ADC)")
parser.add_argument("--sign", type=int, choices=[-1, 1], default=-1, help="Direction of the front to search (+1 or -1)")
parser.add_argument("--npz", type=str, default="waveforms_5GS.npz", help="Name of the .npz file to load")
parser.add_argument("--num_waveforms", type=int, default=5, help="Number of waveforms to display")
parser.add_argument("--lowpass", type=float, help="Cutoff frequency for the lowpass filter (in Hz)")
parser.add_argument("--highpass", type=float, help="Cutoff frequency for the highpass filter (in Hz)")
parser.add_argument("--notch", type=float, help="Frequency for the notch filter (in Hz)")
parser.add_argument("--scanthr", action="store_true", help="Enable threshold scanning mode")
parser.add_argument("--range", type=str, help="Range for scanning: 'min-max' or 'singlethr value'")

args = parser.parse_args()

########################## === Constants === ##########################
ADC_ZERO = 2048                                                       #
SAMPLING_RATE = 5e9                                                   #
#######################################################################

########################## === Data Loading === ##########################
print("\n" + "="*80)
print(Fore.CYAN + "Data Loading Section" + Style.RESET_ALL)
print("="*80 + "\n")
data = np.load(args.npz)
if "waveform_ch1" not in data:
    print(Fore.RED + f"Error: {args.npz} file not found or missing 'waveform_ch1'." + Style.RESET_ALL)
    exit()
print(Fore.GREEN + "File loaded successfully!" + Style.RESET_ALL)
waveforms = data["waveform_ch1"].squeeze()
##########################################################################

filter_obj = Filter(SAMPLING_RATE)
plotter_obj = Plotter()
scan_obj = ScanThreshold() 
num_to_plot = min(args.num_waveforms, len(waveforms))


##################################### === Waveform === ##########################################################
if args.waveform:                                                                                               #
    for i in range(num_to_plot):                                                                                #
        wf = waveforms[i]                                                                                       #
        corr = wf - np.median(wf)                                                                               #
        mv = (corr / ADC_ZERO) * 1000 if args.waveform == "mV" else corr                                        #
        mv = -mv if args.sign == -1 else mv                                                                     #
        if args.lowpass:                                                                                        #
            mv = filter_obj.lowpass(mv, args.lowpass)                                                           #
        if args.highpass:                                                                                       #
            mv = filter_obj.highpass(mv, args.highpass)                                                         #
        if args.notch:                                                                                          #
            mv = filter_obj.notch(mv, args.notch)                                                               #
        plotter_obj.plot_waveform(mv, f"Waveform Displayed {i+1}", "Amplitude", f"Waveform {i+1}", color="blue")#
    exit()                                                                                                      #
#################################################################################################################

############################ === Scan/Fixed threshold === ################################
if args.scanthr and args.range:                                                          #
    print(Fore.CYAN + "\nPre-processing the waveforms for scanning..." + Style.RESET_ALL)#
    preprocessed_waveforms = []                                                          #
    for wf in waveforms:                                                                 #
        corr = wf - np.median(wf)                                                        #
        mv = (corr / ADC_ZERO) * 1000                                                    #
        mv = -mv if args.sign == -1 else mv                                              #
        if args.lowpass:                                                                 #
            mv = filter_obj.lowpass(mv, args.lowpass)                                    #
        if args.highpass:                                                                #
            mv = filter_obj.highpass(mv, args.highpass)                                  #
        if args.notch:                                                                   #
            mv = filter_obj.notch(mv, args.notch)                                        #
        preprocessed_waveforms.append(mv)                                                #
    print(Fore.GREEN + "Pre-processing completed." + Style.RESET_ALL)                    #
    if "singlethr" in args.range:                                                        #
        threshold_value = float(args.range.split()[-1])                                  #
        output_file = "single_scan_results.txt"                                          #
        scan_obj.scan_thresholds(                                                        #
            preprocessed_waveforms,                                                      #
            args.sign,                                                                   #
            threshold_value,                                                             #
            threshold_value,                                                             #
            0.1,                                                                         #
            output_file                                                                  #
        )                                                                                #
        print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)      #
    elif '-' in args.range:                                                              #
        range_min, range_max = map(float, args.range.split('-'))                         #
        output_file = "scan_results.txt"                                                 #
        scan_obj.scan_thresholds(                                                        #
            preprocessed_waveforms,                                                      #
            args.sign,                                                                   #
            range_min,                                                                   #
            range_max,                                                                   #
            0.1,                                                                         #
            output_file                                                                  #
        )                               		                                         #
        print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)      #
##########################################################################################
