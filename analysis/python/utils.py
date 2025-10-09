import re
import numpy as np
import glob
import os
from scipy.optimize import curve_fit

##################### general utility class ############################################################
# This utility module provides a helper function to generate output filenames for threshold scan results.
# The naming convention is based on a voltage value embedded in the input `.npz` filename.
#
# Main Feature (todo: update if necessary):
#
# 1. generate_output_filename(npz_filename):
#    - Extracts the voltage value from the `.npz` filename using a regular expression.
#    - Computes the difference (`delta`) between the extracted voltage and a reference voltage (51.5V).
#    - Formats the difference to use a comma for decimals (e.g., 0,5OV.txt instead of 0.5OV.txt).
#    - Returns a string to be used as the name of the output file, which is intended to store threshold scan results.
#
# Use Cases:
# - Automatically generate standardized filenames when saving threshold scan outputs.
# - Easily identify the voltage setting for a given scan result based on the filename.
# - Avoid manual renaming or mislabeling of analysis outputs.
########################################################################################################


class Utils:
    @staticmethod
    def waveform_to_npz(waveforms, filename, sampling_rate):
        """
        Save waveforms to a .npz file with metadata.
        
        Parameters:
            waveforms (list of np.ndarray): List of waveforms to save.
            filename (str): Name of the output .npz file.
            sampling_rate (float): Sampling rate in MHz? let's see.
        """
        np.savez(filename, waveforms=waveforms, sampling_rate=sampling_rate)
        print(f"Waveforms saved to {filename}") 
    
    @staticmethod
    def generate_output_filename(npz_filename: str) -> str:
        """
        Generates an output filename in the format '[delta]OV.txt' based on the voltage value
        embedded in the input `.npz` filename. The voltage is compared against a fixed reference of 51.5V.

        Args:
            npz_filename (str): Name of the input .npz file (e.g., 'waveforms_750.0_52V.npz').

        Returns:
            str: Formatted output filename (e.g., '0,5OV.txt' or '3OV.txt' for whole numbers).
        """
        match = re.search(r'(\d+(?:\.\d+)?)V', npz_filename)
        if not match:
            raise ValueError(f"Unable to find a voltage value in the filename: {npz_filename}")
        
        voltage = float(match.group(1))
        delta = round(voltage - 51.5, 2)
        
        if delta.is_integer():
            formatted_delta = str(int(delta))
        else:
            formatted_delta = f"{delta}".replace('.', ',')

        return f"{formatted_delta}OV.txt"
    
    
    @staticmethod
    def load_waveforms(npz_path):
        """
        Load waveform data from a given .npz file.
        """
        data = np.load(npz_path)
        return data


    @staticmethod
    def get_info(npz_path):
        """
        Extract sampling rate and other metadata from the filename or path.
        """
        info = {}

        # Extract sampling rate in GHz
        sampling_match = re.search(r'(\d+\.?\d*)GS', npz_path)
        if sampling_match:
            info["sampling_rate_ghz"] = float(sampling_match.group(1))
            info["sampling_rate_hz"] = info["sampling_rate_ghz"] * 1e9
        else:
            info["sampling_rate_ghz"] = None
            info["sampling_rate_hz"] = None

        # Extract bias voltage if present
        bias_match = re.search(r'bias([\d\.]+)', npz_path)
        info["bias"] = bias_match.group(1) if bias_match else 'N/A'

        # Create relative path for saving output
        rel_path = os.path.splitext(npz_path)[0]
        try:
            rel_path = rel_path.split("data/")[1]  # keep only path after 'data/'
        except IndexError:
            rel_path = os.path.basename(rel_path)  # fallback to flat
        info["rel_path"] = rel_path

        # Output path for saving figures or processed data. Let's create a 'plots' directory if it doesn't exist.
        if not os.path.exists("../plots"):
            os.makedirs("../plots")
        info["output_path"] = os.path.join("../plots", rel_path)

        return info
    
    @staticmethod
    def cell_to_seconds(num_cells, sampling = None): 
        """
        Convert a number of DRS4 cells (samples) into time in seconds,
        based on the sampling rate extracted from the filename.
        Parameters:
            num_cells (int or float or np.ndarray): Number of cells to convert.
            filename (str): Filename containing sampling rate in the format '*<rate>GS*'.
        Returns:
            float or np.ndarray: Time duration in seconds.
        """
        seconds = num_cells / sampling
        return seconds

    @staticmethod
    def adc_to_mv(adc_array):
        """
        Convert raw ADC values to millivolts assuming a 12-bit ADC centered at 2048.
        """
        return ((adc_array - 2048) / 4096.0) * 1000  # mV
    
    @staticmethod
    def exponential_decay(t, A, tau, y0):
        return A * np.exp(-t / tau) + y0

    @staticmethod
    def parse_range(col_range_str):
        """Parses a string like '2', '2-4', '2-3,5' into a list of 0-based column indices."""
        indices = set()
        parts = col_range_str.split(',')
        for part in parts:
            if '-' in part:
                start, end = map(int, part.split('-'))
                indices.update(range(start - 1, end))
            else:
                indices.add(int(part) - 1)
        return sorted(indices)

    @staticmethod
    def read_data_file(filepath, num_columns):
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]  # esclude righe vuote

        if not lines:
            raise ValueError(f"File {filepath} è vuoto.")

        header = []
        data_start_idx = 0

        # Cerca la prima riga che inizia con un numero (float)
        for i, line in enumerate(lines):
            first_token = line.split()[0]
            try:
                float(first_token)
                data_start_idx = i
                break
            except ValueError:
                header.append(line)

        data_lines = lines[data_start_idx:]
        if len(data_lines) == 0:
            raise ValueError(f"Nessuna riga dati trovata in {filepath}")

        parsed_data = []
        for line in data_lines:
            parts = line.split()
            if len(parts) != num_columns:
                raise ValueError(f"Incompatible number of columns in {filepath}: expected {num_columns}, got {len(parts)}")
            parsed_data.append([float(x) for x in parts])

        return header, parsed_data


  
    @staticmethod
    def merge_files(files_or_dir, num_columns, sum_columns, output_file):
        if len(files_or_dir) == 1 and os.path.isdir(files_or_dir[0]):
            files = glob.glob(os.path.join(files_or_dir[0], '*'))
        else:
            files = files_or_dir

        if not files:
            raise ValueError("No input files found.")

        all_data = []
        header_saved = None

        for idx, file in enumerate(files):
            header, data = Utils.read_data_file(file, num_columns)

            # Forza data a lista di liste anche se è una sola riga
            if isinstance(data, np.ndarray):
                if data.ndim == 1:
                    data = [data.tolist()]
                else:
                    data = data.tolist()

            if header_saved is None:
                header_saved = header

            if sum_columns:
                if idx == 0:
                    all_data = data
                else:
                    for i, row in enumerate(data):
                        for col in sum_columns:
                            all_data[i][col] += row[col]
            else:
                all_data.extend(data)

        with open(output_file, 'w') as f_out:
            if header_saved:
                for h in header_saved:
                    f_out.write(h + '\n')
            for row in all_data:
                f_out.write('\t'.join(f'{x:.2f}' if isinstance(x, float) else str(x) for x in row) + '\n')
