import sys
import os
import re
import numpy as np
import argparse
from colorama import Fore, Style

'''
================================================================================
Waveform Threshold Scan and Analysis Tool
================================================================================

This script analyzes waveform data (from .npz files) by performing either:

- A threshold scan across a given range, or
- A single-threshold detection mode.

It counts discrete transition events such as pulses (e.g., for DCR studies in
photodetectors).

Features:
- Single-threshold mode or multi-threshold scan ("--scanthr" + "--range")
- Custom waveform region selection ("--x_start", "--x_end")
- Optional filtering: lowpass, highpass, notch
- Input: .npz file containing waveform arrays
- Output: .txt file with results

Requirements:
- The .npz file must contain a key "waveforms".
- Threshold scan range: use 'min-max' format (e.g. --range 5-25).
- Single threshold: use 'singlethr value' (e.g. --range 'singlethr 12.5').
- Sampling rate is entered in MHz on the command line.
================================================================================
'''

########################## === Libraries' path === ##########################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))

from filters import Filters
from plotter import Plotter
from scanthr import ScanThreshold
from utils import Utils
from waveform_analysis import waveform_analysis
#################################################################################

########################## === Argument Parsing === ##########################
parser = argparse.ArgumentParser(
    description="Waveform analysis with thresholds and filters.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog="""
Examples of usage:

  1) Run with a single threshold of 12.5 mV:
     python script.py --npz waveforms_5GS.npz --scanthr --range "singlethr 12.5" --output results.txt

  2) Run a threshold scan between 5 and 25 mV:
     python script.py --npz waveforms_5GS.npz --scanthr --range 5-25 --output scan.txt

  3) Apply a lowpass filter at 20 MHz and scan:
     python script.py --npz waveforms_5GS.npz --scanthr --range 10-30 --lowpass 2e7 --output filtered.txt

Sampling rate is given in MHz, e.g. --sampling_rate 750
""",
)

parser.add_argument("--sign", type=int, choices=[-1, 1], default=1,
                    help="Signal polarity: -1 for negative, +1 for positive")
parser.add_argument("--npz", type=str, default="waveforms_5GS.npz",
                    help="Path to the .npz file containing waveform data")
parser.add_argument("--lowpass", type=float,
                    help="Lowpass filter cutoff frequency (Hz)")
parser.add_argument("--highpass", type=float,
                    help="Highpass filter cutoff frequency (Hz)")
parser.add_argument("--notch", type=float,
                    help="Notch filter frequency (Hz)")
parser.add_argument("--scanthr", action="store_true",
                    help="Enable threshold scan mode")
parser.add_argument("--range", type=str,
                    help="Threshold range: 'min-max' or 'singlethr value'")
parser.add_argument("--x_start", type=int, default=0,
                    help="Start index of waveform region to process")
parser.add_argument("--x_end", type=int, default=1024,
                    help="End index of waveform region to process")
parser.add_argument("--sampling_rate", type=float, default=750.0,
                    help="Sampling rate in MHz (will be converted to Hz internally)")
parser.add_argument("--step_voltage_mV", type=float, default=0.20,
                    help="Step size in mV for threshold scanning")
parser.add_argument("--output", type=str, required=True,
                    help="Output file name for results")

args = parser.parse_args()
#################################################################################

########################## === Load Data === ##########################
print("\n" + "="*80)
print(Fore.CYAN + "Data Loading Section" + Style.RESET_ALL)
print("="*80 + "\n")

data = Utils.load_waveforms(args.npz)
if "waveforms" not in data:
    print(Fore.RED + f"Error: {args.npz} file not found or missing 'waveforms'." + Style.RESET_ALL)
    sys.exit(1)

print(Fore.GREEN + "File loaded successfully!" + Style.RESET_ALL)
waveforms = data["waveforms"].squeeze()
#################################################################################

########################## === Initialize Tools === ##########################
sampling_rate_hz = args.sampling_rate * 1e6  # convert MHz → Hz
filter_obj = Filters(sampling_rate_hz)
plotter_obj = Plotter()
scan_obj = ScanThreshold()
#################################################################################

########################## === Threshold Scan === ##########################
if args.scanthr and args.range:
    print(Fore.CYAN + "\nPre-processing the waveforms for threshold scanning..." + Style.RESET_ALL)
    preprocessed_waveforms = []

    for wf in waveforms:
        mv = Utils.adc_to_mv(wf)
        baseline = np.median(mv[49:973])
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
        scan_obj.scan_thresholds(preprocessed_waveforms, args.sign,
                                 threshold_value, threshold_value,
                                 args.step_voltage_mV, output_file)
    elif '-' in args.range:
        range_min, range_max = map(float, args.range.split('-'))
        scan_obj.scan_thresholds(preprocessed_waveforms, args.sign,
                                 range_min, range_max,
                                 args.step_voltage_mV, output_file)

    print(Fore.YELLOW + f"\nResults saved in: {output_file}" + Style.RESET_ALL)
#################################################################################
