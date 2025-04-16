import re

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
    def generate_output_filename(npz_filename: str) -> str:
        """
        Args:
            npz_filename (str): Name of the input .npz file (e.g., 'waveforms_750.0_52V.npz').

        Returns:
            str: Formatted output filename (e.g., '0,5OV.txt' or '3OV.txt' for whole numbers).
        """
        # Extract voltage value using regex (e.g., matches '52V' as 52)
        match = re.search(r'(\d+(?:\.\d+)?)V', npz_filename)
        if not match:
            raise ValueError(f"Unable to find a voltage value in the filename: {npz_filename}")
        
        voltage = float(match.group(1))
        delta = round(voltage - 51.5, 2)  # Compute delta from reference voltage

        # Format delta: use comma instead of dot for decimals; remove decimal for integers
        if delta.is_integer():
            formatted_delta = str(int(delta))
        else:
            formatted_delta = f"{delta}".replace('.', ',')

        return f"{formatted_delta}OV.txt"
