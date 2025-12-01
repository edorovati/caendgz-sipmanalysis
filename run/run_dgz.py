import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import uproot
from tqdm import tqdm  # Barra di avanzamento
from datetime import datetime
import os
import subprocess
import time
import threading
import queue

sys.path.append("/eu/caen-dt5742b/python/")
from rwave import rwaveclient

parser = argparse.ArgumentParser(description="Receive and save waveform data from CAEN digitizer.")
parser.add_argument("--filter_ADC", type=float, default=None,
                    help="Apply filter on central peak-to-peak amplitude (ADC units).")
parser.add_argument("--channel", type=int, nargs='+', default=None, 
                    help="Specify channel(s) to acquire data from (e.g., --channel 1 3 4).")
parser.add_argument("--min_events", type=int, default=1000, 
                    help="Minimum number of valid waveforms to accumulate before saving.")
parser.add_argument("--log_file", type=str, default="acquisition_log.txt", 
                    help="Log file to record acquisition details.")
parser.add_argument("--sampling", type=float, default=5000,
                    help="Sampling frequency in MHz (default: 5.0 MHz)")
parser.add_argument("--vbias", type=float, required=True,
    help="Bias voltage applied to the SiPM (in volts, e.g., 35)")
parser.add_argument("--folder", type=str, default=None,
                    help="Optional subfolder in /data where to save output files.")
parser.add_argument('--trg', choices=['NIM', 'sw'], default='sw', help="Trigger type: 'NIM' 'laser' (with NIM) or 'sw'")
parser.add_argument("--shift_waveforms", action="store_true", help="Shift waveforms by first_cell modulo 1024")
parser.add_argument('--filename', type=str, default='output_data.npz',
                    help='Name of the .npz file to save (default: output_data.npz)')
parser.add_argument('--plot', action='store_true', help='Plot the waveforms after acquisition')
args = parser.parse_args()

# safety check for bias voltage
if not (10.0 <= args.vbias <= 80.0):
    print(f"⚠️ Warning: Vbias {args.vbias} V seems out of expected range (10–80 V).")


# Create data directory if it doesn't exist
os.makedirs("../data", exist_ok=True)

# Create base data directory
base_dir = "../data"

# Check if a subfolder was specified
if args.folder:
    save_dir = os.path.join(base_dir, args.folder)
    os.makedirs(save_dir, exist_ok=True)
    print(f"Saving data in subfolder: {save_dir}")
else:
    save_dir = base_dir


HOST = 'localhost'
PORT = 30001

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
NPZ_FILE = os.path.join(save_dir, f"{timestamp}_waveforms_bias{args.vbias}V_{int(args.sampling)/1000}GS.npz")
TREE_NAME = "waveform_tree"


def command_listener(command_queue):
    print("⌨️  Type 'start' to begin, 'stop' to pause, 'exit' to quit.")
    while True:
        cmd = input(">> ").strip().lower()
        command_queue.put(cmd)
        if cmd == 'exit':
            break

def NIM_trg():
    """Set up the NIM TRG with the pulser."""
    subprocess.run("/eu/aimtti/aimtti-cmd.py --address aimtti-tgp3152-00 --list /eu/aimtti-tgp3152/configs/NIM_trg.config", shell=True, check=True)
    return True
    
def laser_trg():
    """Set up the NIM TRG with the pulser."""
    subprocess.run("/eu/aimtti/aimtti-cmd.py --address aimtti-tgp3152-00 --list /home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/configs/laser_trg.conf", shell=True, check=True)
    return True

def acquire_data_sw(chmask, correction):
    """Acquire waveform data from the digitizer."""
    with rwaveclient(HOST, PORT, verbose=True) as rwc:
        if rwc is None:
            return None

        rwc.send_cmd(f'sampling {int(args.sampling)}')
        rwc.send_cmd('grmask 0x1')
        rwc.send_cmd(f'chmask {chmask}')
        rwc.send_cmd('correction on' if correction else 'correction off')
        rwc.send_cmd("start")
        rwc.send_cmd('swtrg 1024')
        rwc.send_cmd('readout')
        rwc.send_cmd('download')
        data = rwc.download()
        rwc.send_cmd('stop')

        return data

def configure_dgz(rwc,chmask, correction):
    """Configure the digitizer settings."""
    rwc.send_cmd(f'sampling {int(args.sampling)}')
    rwc.send_cmd('grmask 0x1')
    rwc.send_cmd(f'chmask {chmask}')
    rwc.send_cmd('correction on' if correction else 'correction off')
    rwc.send_cmd("start")
    return rwc

def acquire_data_NIM(rwc):
    """Acquire waveform data from the digitizer."""
      
    rwc.send_cmd("start")
        
    time.sleep(3)

    rwc.send_cmd('readout')
    rwc.send_cmd('download')
    data = rwc.download()
    
    return data

def close_acquisition(rwc):
    rwc.send_cmd('stop')
    return True
    

def handle_data(data, selected_ch=None, shift_waveforms=False):
    """Process waveform data and return one dictionary per event.
    If shift_waveforms is True, shift waveform array according to first_cell modulo 1024."""
    if data is None:
        print("No data received.")
        return None

    all_events_data = []
    max_channels = len(data[0])

    if selected_ch is not None:
        active_channels = selected_ch
    else:
        active_channels = range(max_channels)

    for event_number, event in enumerate(data):
    #     first_cell = event[0]["first_cell"] % 1024 if shift_waveforms else 0
        
        for ch in active_channels:
            waveform = event[ch]["waveform"]
            # if shift_waveforms:
            #     waveform = np.roll(waveform, -first_cell)
            
            tree_data = {
                "event_number": np.array([event_number]),
                f'waveform_ch{ch}': np.array([waveform]),
                # "trigger_tag": np.array([event[0]["trigger_tag"]]),
                # "first_cell": np.array([event[0]["first_cell"]]),
                # "num_cells": np.array([np.arange(len(waveform))])
            }
            all_events_data.append(tree_data)

    return all_events_data

def apply_filter(event_data, threshold=20, window_size=900):
    """Filter waveform events based on central peak-to-peak amplitude."""
    start_idx = (1024 - window_size) // 2
    end_idx = start_idx + window_size
    filtered_dict = {}

    for event in event_data:
        for key in event:
            if key.startswith("waveform_ch"):
                ch = key.split("waveform_ch")[1]
                wf = event[key].squeeze()
                central_window = wf[start_idx:end_idx]
                peak_to_peak = central_window.max() - central_window.min()
                if peak_to_peak > threshold:
                    if f'waveform_ch{ch}' not in filtered_dict:
                        filtered_dict[f'waveform_ch{ch}'] = []
                    filtered_dict[f'waveform_ch{ch}'].append(wf)

    return filtered_dict


def save_filtered_waveforms_to_npz(filtered_dict, output_file):
    """Save filtered waveforms to compressed .npz file."""
    if filtered_dict:
        for ch_key, wfs in filtered_dict.items():
            print(f"{ch_key}: {len(wfs)} waveforms passed the filter.")
        np.savez_compressed(output_file, **{k: np.array(v) for k, v in filtered_dict.items()})
        print(f"Filtered NPZ file saved as: {output_file}")
    else:
        print("⚠️ No waveform passed the filter.")
    return bool(filtered_dict)


def save_waveforms_to_root(filtered_dict, output_file):
    """Save filtered waveforms into separate ROOT files per channel."""
    for ch_key, waveforms in filtered_dict.items():
        ch = ch_key.replace("waveform_ch", "")
        with uproot.recreate(f'ch{ch}_{output_file}') as file:
            for i, wf in enumerate(waveforms):
                file[f"{TREE_NAME}_event_{i}"] = {
                    "event_number": np.array([i]),
                    f"waveform_ch{ch}": np.array([wf]),
                    # "trigger_tag": np.array([0]),
                    # "first_cell": np.array([0]),
                    # "num_cells": np.array([np.arange(len(wf))])
                }
        print(f"ROOT file saved: ch{ch}_{output_file}")


# def print_trigger_tags(data, selected_ch=None):
#     """Print the trigger tags for each waveform."""
#     if data is None:
#         print("No data received.")
#         return

#     max_channels = len(data[0])

#     if selected_ch is not None:
#         active_channels = selected_ch
#     else:
#         active_channels = range(max_channels)

#     for event_number, event in enumerate(data):
#         for ch in active_channels:
#             trigger_tag = event[0]["trigger_tag"]
#             print(f"Event {event_number}, Channel {ch}: Trigger Tag = {trigger_tag}")

if __name__ == "__main__":
    # if args.trg == 'NIM':
    #     print("NIM pulser should be already ON from the .sh script...")
    #     NIM_trg()
    # elif args.trg == 'laser':
    #     print("Starting laser pulser...")
    #     laser_trg()

    # Print selected channels
    selected_channels = args.channel if args.channel else []
    print(f"Selected channels: {selected_channels}")

    # Warn if no channels specified
    if not selected_channels:
        print("⚠️ Warning: No channels specified, acquisition may fail or use default channels.")

    # Compute the channel mask from selected channels
    chmask = 0
    for ch in selected_channels:
        chmask |= (1 << ch)
    print(f"Using channel mask: 0x{chmask:04x}")

    min_events = args.min_events
    print(f"Accumulating at least {min_events} valid waveforms per channel.")

    # Initialize dictionary to store valid waveforms per channel
    valid_waveforms = {ch: [] for ch in selected_channels}

    run_count = 0

    # Connect to the digitizer client
    with rwaveclient(HOST, PORT, verbose=True) as rwc:
        if rwc is None:
            sys.exit("❌ Failed to connect to the digitizer.")

        print("Configuring digitizer...")
        close_acquisition(rwc)

        # Configure digitizer with dynamic channel mask
        configure_dgz(rwc, chmask, correction=True)

        with open(args.log_file, 'w') as log_file:
            log_file.write(f"Acquisition started. Accumulating at least {min_events} waveforms per channel.\n")
            log_file.write(f"Selected channels: {selected_channels}\n")

            # Progress bar for acquisition
            pbar = tqdm(total=min_events, desc="Min waveforms per channel", ncols=100, unit="waveforms")

            # Continue acquisition until min_events reached on all selected channels
            while min(len(valid_waveforms[ch]) for ch in selected_channels) < min_events:
                run_count += 1
                min_valid = min(len(valid_waveforms[ch]) for ch in selected_channels)
                print(f"\nRun {run_count}: Accumulating... (min={min_valid} / {min_events})")

                # Acquire data depending on trigger type
                if args.trg == 'NIM':
                    print("🔌 Starting acquisition with NIM trigger...")
                    data = acquire_data_NIM(rwc)
                else:
                    print("🔌 Starting acquisition with software trigger...")
                    data = acquire_data_sw(chmask, correction=True)

                if data is None:
                    sys.exit("❌ Acquisition failed.")

                # Process acquired data, optionally shift waveforms
                event_data = handle_data(data, selected_ch=selected_channels, shift_waveforms=args.shift_waveforms)
                if event_data is None:
                    sys.exit("❌ No events to process.")

                # Apply filtering if filter_ADC is set
                if args.filter_ADC is not None:
                    print(f"🔍 Applying filter: peak-to-peak > {args.filter_ADC} ADC")
                    filtered = apply_filter(event_data, threshold=args.filter_ADC)
                    for ch_key, waveforms in filtered.items():
                        ch = int(ch_key.replace("waveform_ch", ""))
                        if ch in valid_waveforms:
                            valid_waveforms[ch].extend(waveforms)
                else:
                    print("No filter applied. Saving all waveforms.")
                    for event in event_data:
                        for ch_key in event:
                            if ch_key.startswith("waveform_ch"):
                                ch = int(ch_key.replace("waveform_ch", ""))
                                if ch in valid_waveforms:
                                    valid_waveforms[ch].append(event[ch_key])

                # Update progress bar
                min_valid = min(len(valid_waveforms[ch]) for ch in selected_channels)
                pbar.update(min_valid - pbar.n)
                log_file.write(f"Run {run_count}: waveforms per channel = { {ch: len(valid_waveforms[ch]) for ch in selected_channels} }\n")

            print(f"\n✅ Done. At least {min_events} waveforms per channel acquired.")

            # Save filtered waveforms per channel as compressed NPZ files
            for ch_key, waves in valid_waveforms.items():
                ch_filename = f"{args.filename}_ch{ch_key}.npz"
                ch_output_path = os.path.join(save_dir, ch_filename)
                np.savez_compressed(ch_output_path, waveforms=np.array(waves))
                print(f"🔸 Saved {len(waves)} waveforms for channel {ch_key} to {ch_output_path}")

            # Optional plotting of waveforms
            if args.plot:
                plt.figure(figsize=(14, 8))
                x_time = (np.arange(1024) / args.sampling * 1e3)  # time in ns
                max_plot_waves = 20000
                for ch in selected_channels:
                    waves = valid_waveforms[ch]
                    for i in range(min(len(waves), max_plot_waves)):
                        plt.plot(x_time, waves[i], alpha=0.5, label=f'ch{ch} Waveform {i+1}')
                plt.title('Sample Waveforms per Channel')
                plt.xlabel('Time [ns]')
                plt.ylabel('Amplitude [ADC]')
                plt.grid()
                plt.tight_layout()
                plt.show()

            print("Waveforms saved in per-channel .npz files.")
            pbar.close()
            close_acquisition(rwc)