import re
import numpy as np
import os

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
    def find_data_directory(start_path=None, target_folder="data", max_depth=5):
        """
        Search upward from the given path for a folder named `data`.

        :param start_path: starting directory (default: directory of this file)
        :param target_folder: name of the folder to find
        :param max_depth: how many levels to go up
        :return: absolute path to the found folder or raises FileNotFoundError
        """
        if start_path is None:
            start_path = os.path.dirname(__file__)
        current_path = os.path.abspath(start_path)

        for _ in range(max_depth):
            candidate = os.path.join(current_path, target_folder)
            if os.path.isdir(candidate):
                return candidate
            # Go one level up
            current_path = os.path.dirname(current_path)

        raise FileNotFoundError(f"Folder '{target_folder}' not found within {max_depth} levels from {start_path}")

    
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

        # Output path for saving figures or processed data
        info["output_path"] = os.path.join("analysis", rel_path + ".png")

        return info
    
    @staticmethod
    def cell_to_seconds(num_cells, filename): 
        """
        Convert a number of DRS4 cells (samples) into time in seconds,
        based on the sampling rate extracted from the filename.
        Parameters:
            num_cells (int or float or np.ndarray): Number of cells to convert.
            filename (str): Filename containing sampling rate in the format '*<rate>GS*'.
        Returns:
            float or np.ndarray: Time duration in seconds.
        """
        sampling_rate_hz = Utils.get_info(filename)["sampling_rate_hz"]
        seconds = num_cells / sampling_rate_hz
        return seconds

    @staticmethod
    def adc_to_mv(adc_array):
        """
        Convert raw ADC values to millivolts assuming a 12-bit ADC centered at 2048.
        """
        return ((adc_array - 2048) / 4096.0) * 1000  # mV