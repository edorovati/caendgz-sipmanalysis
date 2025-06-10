import sys
import os
import re
import numpy as np
import argparse
from colorama import Fore, Style

'''                            '''                                         '''
================================================================================
Waveform Threshold Scan and Analysis Tool
================================================================================

This script is designed to analyze waveform data (typically from .npz files)
by performing a **threshold scan** or a **single-threshold detection** to 
count discrete transition events such as pulses (e.g., for DCR — Dark Count Rate 
studies in photodetectors). The script supports:

- Single-threshold mode or a multiple-threshold scan (`--scanthr` + `--range`)
- Custom selection of waveform regions (`--x_start`, `--x_end`)
- Optional filtering: lowpass, highpass, or notch filter
- Input as `.npz` file containing waveform arrays
- Output as a `.txt` file containing results from the scan

Usage is configured entirely via command-line arguments.

Requirements:
- The waveform file must contain a key `"waveforms"`
- Threshold scan range: use `'min-max'` format (e.g., `--range 5-25`)
- Single threshold: use `'singlethr value'` (e.g., `--range 'singlethr 12.5'`)
- Output is saved to a `.txt` file specified by the `--output` argument

================================================================================
'''                            '''                                         '''

########################## === Libraries' path === ##########################
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/Users/edoardorovati/Desktop/caendgz-sipmanalysis/analysis/python'))
sys.path.append(base_path)
from filters import Filters
from plotter import Plotter
from scanthr import ScanThreshold
from utils import Utils
from waveform_analysis import waveform_analysis
#################################################################################

########################## === Argument Parsing === ##########################
parser = argparse.ArgumentParser(description="Waveform analysis with thresholds and filters")
parser.add_argument("--sign", type=int, choices=[-1, 1], default=-1, help="Negative or positive signal (-1 or +1)")
parser.add_argument("--npz", type=str, default="waveforms_5GS.npz", help="Path to the .npz file containing waveform data")
parser.add_argument("--num_waveforms", type=int, default=5, help="Number of waveforms to process")
parser.add_argument("--lowpass", type=float, help="Lowpass filter cutoff frequency (Hz)")
parser.add_argument("--highpass", type=float, help="Highpass filter cutoff frequency (Hz)")
parser.add_argument("--notch", type=float, help="Notch filter frequency (Hz)")
parser.add_argument("--scanthr", action="store_true", help="Enable threshold scan mode")
parser.add_argument("--range", type=str, help="Threshold range: 'min-max' or 'singlethr value'")
parser.add_argument("--x_start", type=int, default=0, help="Start index of waveform region to process")
parser.add_argument("--x_end", type=int, default=1024, help="End index of waveform region to process")
parser.add_argument("--sampling_rate", type=float, default=750e6, help="Sampling rate in Hz")
parser.add_argument("--output", type=str, required=True, help="Output file name")

args = parser.parse_args()
#################################################################################

########################## === Load Data === ##########################
print("\n" + "="*80)
print(Fore.CYAN + "Data Loading Section" + Style.RESET_ALL)
print("="*80 + "\n")

data = Utils.load_waveforms(args.npz)
if "waveforms" not in data:
    print(Fore.RED + f"Error: {args.npz} file not found or missing 'waveforms'." + Style.RESET_ALL)
    exit()

print(Fore.GREEN + "File loaded successfully!" + Style.RESET_ALL)
waveforms = data["waveforms"].squeeze()
#################################################################################

########################## === Initialize Tools === ##########################
filter_obj = Filters(args.sampling_rate)
plotter_obj = Plotter()
scan_obj = ScanThreshold()
num_to_plot = min(args.num_waveforms, len(waveforms))
#################################################################################

########################## === Threshold Scan === ##########################
if args.scanthr and args.range:
    print(Fore.CYAN + "\nPre-processing the waveforms for threshold scanning..." + Style.RESET_ALL)
    preprocessed_waveforms = []
    for wf in waveforms:
        mv = Utils.adc_to_mv(wf)
        baseline,_ = waveform_analysis.calculate_baseline_with_mask(mv, 49, 973)
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

    output_file = args.output

    if "singlethr" in args.range:
        threshold_value = float(args.range.split()[-1])
        scan_obj.scan_thresholds(preprocessed_waveforms, args.sign, threshold_value, threshold_value, 0.50, output_file)
    elif '-' in args.range:
        range_min, range_max = map(float, args.range.split('-'))
        scan_obj.scan_thresholds(preprocessed_waveforms, args.sign, range_min, range_max, 0.5, output_file)
    print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)
#################################################################################
