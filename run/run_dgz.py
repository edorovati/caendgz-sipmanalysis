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
parser.add_argument('--trg', choices=['NIM', 'sw'], default='sw', help="Trigger type: 'NIM' or 'sw'")
parser.add_argument("--shift_waveforms", action="store_true", help="Shift waveforms by first_cell modulo 1024")
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
OUTPUT_FILE = os.path.join(save_dir, f"{timestamp}_waveforms_bias{args.vbias}V_{int(args.sampling)/1000}GS.root")
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
    

def acquire_data(chmask, correction):
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

    # rwc.send_cmd(f'sampling {int(args.sampling)}')
    # rwc.send_cmd('grmask 0x1')
    # rwc.send_cmd(f'chmask {chmask}')
    # rwc.send_cmd('correction on' if correction else 'correction off')
      
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
        first_cell = event[0]["first_cell"] % 1024 if shift_waveforms else 0
        
        for ch in active_channels:
            waveform = event[ch]["waveform"]
            if shift_waveforms:
                waveform = np.roll(waveform, -first_cell)
            
            tree_data = {
                "event_number": np.array([event_number]),
                f'waveform_ch{ch}': np.array([waveform]),
                "trigger_tag": np.array([event[0]["trigger_tag"]]),
                "first_cell": np.array([event[0]["first_cell"]]),
                "num_cells": np.array([np.arange(len(waveform))])
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
                    "trigger_tag": np.array([0]),
                    "first_cell": np.array([0]),
                    "num_cells": np.array([np.arange(len(wf))])
                }
        print(f"ROOT file saved: ch{ch}_{output_file}")


def print_trigger_tags(data, selected_ch=None):
    """Print the trigger tags for each waveform."""
    if data is None:
        print("No data received.")
        return

    max_channels = len(data[0])

    if selected_ch is not None:
        active_channels = selected_ch
    else:
        active_channels = range(max_channels)

    for event_number, event in enumerate(data):
        for ch in active_channels:
            trigger_tag = event[0]["trigger_tag"]
            print(f"Event {event_number}, Channel {ch}: Trigger Tag = {trigger_tag}")







# if __name__ == "__main__":
#     selected_channels = args.channel if args.channel else None
#     min_events = args.min_events
#     valid_waveforms = []
#     run_count = 0

#     command_queue = queue.Queue()
#     listener_thread = threading.Thread(target=command_listener, args=(command_queue,))
#     listener_thread.daemon = True
#     listener_thread.start()
#     plotted_waveforms = []
#     with rwaveclient(HOST, PORT, verbose=True) as rwc:
#         if rwc is None:
#             sys.exit("❌ Failed to connect to the digitizer.")
        
#         print("🔌 Configuring digitizer...")
#         close_acquisition(rwc)

#         with open(args.log_file, 'w') as log_file:
#             log_file.write(f"Selected channels: {selected_channels}\n")

#             while True:
#                 try:
#                     cmd = command_queue.get(timeout=0.1)
#                     if cmd == 'exit':
#                         print("👋 Exiting...")
#                         plt.figure(figsize=(10, 6))
#                         for i in range(min(100, len(plotted_waveforms))):
#                             plt.plot(plotted_waveforms[i].squeeze(), label=f'Waveform {i+1}')
#                         plt.title('Sample Waveforms')
#                         plt.xlabel('Sample Number')
#                         plt.ylabel('Amplitude (ADC)')
#                         plt.legend()
#                         plt.grid()
#                         plt.show()

#                         break
                    
#                     elif cmd == 'start':
#                         print("▶️  Starting acquisition...")
#                         run_count += 1
#                         pbar = tqdm(total=min_events, desc="Accumulating waveforms", ncols=100, unit="waveforms")
#                         configure_dgz(rwc, 0x0003, correction=True)

#                         while len(valid_waveforms) < min_events:
#                             # check if user typed 'stop'
#                             if not command_queue.empty():
#                                 subcmd = command_queue.get_nowait()
#                                 if subcmd == 'stop':
#                                     print("⏸️  Acquisition paused.")
#                                     break
#                                 elif subcmd == 'exit':
#                                     raise KeyboardInterrupt
                            
#                             print(f"\nRun {run_count}: Accumulating waveforms... ({len(valid_waveforms)} / {min_events})")
#                             if args.trg == 'NIM':
#                                 data = acquire_data_NIM(rwc)
#                             else:
#                                 data = acquire_data(0x0003, correction=True)

#                             if data is None:
#                                 sys.exit("❌ Acquisition failed.")

#                             event_data = handle_data(data, selected_ch=selected_channels, shift_waveforms=args.shift_waveforms)
#                             if event_data is None:
#                                 sys.exit("❌ No events to process.")

#                             if args.filter_ADC is not None:
#                                 filtered = apply_filter(event_data, threshold=args.filter_ADC)
#                                 for ch_key, waveforms in filtered.items():
#                                     valid_waveforms.extend(waveforms)
#                             else:
#                                 for event in event_data:
#                                     for ch_key in event:
#                                         if ch_key.startswith("waveform_ch"):
#                                             valid_waveforms.append(event[ch_key])

#                             pbar.update(len(valid_waveforms) - pbar.n)
#                             log_file.write(f"Run {run_count}: {len(valid_waveforms)} valid waveforms accumulated\n")

#                         pbar.close()

#                         # plot and save only after a run is completed
#                         if len(valid_waveforms) >= min_events:
#                             plt.figure(figsize=(10, 6))
#                             for i in range(min(100, len(valid_waveforms))):
#                                 plt.plot(valid_waveforms[i].squeeze(), label=f'Waveform {i+1}')
#                                 plotted_waveforms.append(valid_waveforms[i])
#                             plt.title('Sample Waveforms')
#                             plt.xlabel('Sample Number')
#                             plt.ylabel('Amplitude (ADC)')
#                             plt.legend()
#                             plt.grid()
#                             plt.show()

#                             print(f"\n💾 {len(valid_waveforms)} valid waveforms accumulated. Saving to NPZ.")
#                             np.savez_compressed(NPZ_FILE, waveforms=np.array(valid_waveforms))
#                             log_file.write(f"Acquisition complete. {len(valid_waveforms)} waveforms saved.\n")
                            
#                             valid_waveforms.clear()

#                 except queue.Empty:
#                     continue
#                 except KeyboardInterrupt:
#                     print("👋 Interrupted. Closing acquisition.")
#                     break

#         close_acquisition(rwc)
#         print(f"Waveforms saved in: {NPZ_FILE}")

#         #plot the waveforms
#         if plotted_waveforms:
#             plt.figure(figsize=(10, 6))
#             for i in range(min(100, len(plotted_waveforms))):
#                 plt.plot(plotted_waveforms[i].squeeze(), label=f'Waveform {i+1}')
#             plt.title('Sample Waveforms')
#             plt.xlabel('Sample Number')
#             plt.ylabel('Amplitude (ADC)')
#             plt.legend()
#             plt.grid()
#             plt.show()

if __name__ == "__main__":
    if args.trg == 'NIM':
        # print("🔌 Starting NIM pulser...")
        print("NIM pulser should be already ON from the .sh script...")
        NIM_trg()

    selected_channels = args.channel if args.channel else None
    print(f"Selected channels: {selected_channels}")

    min_events = args.min_events
    print(f"Accumulating at least {min_events} valid waveforms.")

    valid_waveforms = []
    run_count = 0

    with rwaveclient(HOST, PORT, verbose=True) as rwc:
        if rwc is None:
            sys.exit("❌ Failed to connect to the digitizer.")
        print("🔌 Configuring digitizer...")
        close_acquisition(rwc)
            
        with open(args.log_file, 'w') as log_file:
            log_file.write(f"Acquisition started. Accumulating at least {min_events} waveforms.\n")
            log_file.write(f"Selected channels: {selected_channels}\n")

            pbar = tqdm(total=min_events, desc="Accumulating waveforms", ncols=100, unit="waveforms")

            configure_dgz(rwc, 0x0003, correction=True)
            
            while len(valid_waveforms) < min_events:
                run_count += 1
                print(f"\nRun {run_count}: Accumulating waveforms... ({len(valid_waveforms)} / {min_events})")
                
                # Qui scegli la funzione in base al trigger
                if args.trg == 'NIM':
                    data = acquire_data_NIM(rwc)
                else:
                    data = acquire_data(0x0003, correction=True)
                
                if data is None:
                    sys.exit("❌ Acquisition failed.")

                # print_trigger_tags(data, selected_ch=selected_channels)

                event_data = handle_data(data, selected_ch=selected_channels, shift_waveforms=args.shift_waveforms)
                if event_data is None:
                    sys.exit("❌ No events to process.")

                if args.filter_ADC is not None:
                    print(f"🔍 Applying filter: peak-to-peak > {args.filter_ADC} ADC")
                    filtered = apply_filter(event_data, threshold=args.filter_ADC)

                    for ch_key, waveforms in filtered.items():
                        valid_waveforms.extend(waveforms)
                    print(f"Waveforms after filtering: {len(valid_waveforms)}")
                else:
                    print("No filter applied. Saving all waveforms.")
                    for event in event_data:
                        for ch_key in event:
                            if ch_key.startswith("waveform_ch"):
                                valid_waveforms.append(event[ch_key])

                pbar.update(len(valid_waveforms) - pbar.n)
                log_file.write(f"Run {run_count}: {len(valid_waveforms)} valid waveforms accumulated\n")

            
            print(f"\n{len(valid_waveforms)} valid waveforms accumulated. Saving to NPZ.")
            np.savez_compressed(NPZ_FILE, waveforms=np.array(valid_waveforms))
            log_file.write(f"Acquisition complete. {len(valid_waveforms)} valid waveforms saved.\n")
            pbar.close()
            close_acquisition(rwc)
            

        print(f"Waveforms saved in: {NPZ_FILE}")
        
        
        
