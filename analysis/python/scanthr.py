import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
# Aggiungi queste due righe subito dopo gli altri import
import uproot


##################### Threshold Scanning for Waveform Analysis ###############################
# This module provides functionality to analyze waveforms by applying threshold-based scans.
# It is designed to detect transitions (i.e., threshold crossings) in waveform signals and
# to evaluate how many such transitions occur for different threshold levels.
#
# Main Features:
#
# 1. get_transitions(waveform_mv, threshold, sign):
#    - Detects the indices where the waveform crosses the threshold in the direction specified by sign.
#    - This function implements a basic threshold-crossing detection mechanism with an "arming" logic
#        to ensure that it only detects one crossing per pulse or excursion above the threshold.
#    - The logic is as follows:
#        * Wait for the signal to go below (or above, depending on sign) 50% of the threshold (arming phase).
#        * Once "armed", detect the first point where the signal crosses the threshold.
#        * Record the index of this crossing, then disarm until the waveform resets below the 50% level again
#    - It returns the indices (sample numbers) where the signal crosses a defined threshold.
#    - The `sign` parameter determines the direction of the crossing:
#        * sign = +1 => looks for upward crossings (e.g., negative pulse going up),
#        * sign = -1 => looks for downward crossings (e.g., positive pulse going down).
#    - It uses a simple "arming" mechanism to avoid detecting multiple transitions too close together.
#
# 2. scan_thresholds(waveforms, sign, range_min, range_max, step, output_file):
#    - Scans through a range of thresholds and counts how many transitions are detected at each level.
#    - Saves the results to a text file (threshold vs. count), useful for hit rate analysis or noise studies.
#
# Use Cases:
# - Determine the optimal detection threshold in noisy signals.
# - Analyze how signal detectability changes with different filtering.
# - Produce threshold vs. event count plots for calibration.
###############################################################################################


class ScanThreshold:
    def __init__(self):
        pass

    @staticmethod
    def get_transitions(waveform_mv: np.ndarray, threshold: float, sign: int) -> list[int]:
        """
        :param waveform_mv: The waveform data as a numpy array (in millivolts).
        :param threshold: The threshold level (in millivolts) to detect crossings.
        :param sign: The direction of threshold crossing:
                     +1 to detect upward crossings (e.g., negative pulse going up),
                     -1 to detect downward crossings (e.g., positive pulse going down).
        :return: A list of indices (samples) where threshold crossings occurred.
        """
        transitions = []
        armed = False
        for i, y in enumerate(waveform_mv):
            if not armed and sign * y > sign * threshold * 0.5:
                continue
            armed = True
            if sign * y < sign * threshold:
                continue
            transitions.append(i)
            armed = False
        return transitions

    def scan_thresholds(self, waveforms_mv: list[np.ndarray], sign: int, range_min: float, range_max: float, step: float, output_file: str):
        """
        :param waveforms_mv: Array of waveforms already in mV.
        :param sign: Direction of threshold (+1 or -1).
        :param range_min: Minimum threshold (in mV).
        :param range_max: Maximum threshold (in mV).
        :param step: Step between thresholds (in mV).
        :param output_file: File .txt to save the results.
        """
        with open(output_file, 'w') as f:
            f.write("Threshold (mV)   Counts\n")
            for threshold_mv in np.arange(range_min, range_max + step, step):
                count_above_threshold = 0
                for wf_mv in waveforms_mv:
                    transitions = self.get_transitions(wf_mv, threshold_mv, sign)
                    count_above_threshold += len(transitions)
                f.write(f"{threshold_mv:.2f}    {count_above_threshold}\n")
                print(f"Threshold: {threshold_mv:.2f}, Counts: {count_above_threshold}")
                
    
