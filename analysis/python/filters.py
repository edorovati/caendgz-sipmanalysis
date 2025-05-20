import numpy as np
from scipy.signal import butter, filtfilt, iirnotch

##################### Filter Class for Waveform Analysis ###############################
# This class provides various filter functions to process waveform data.
# Filters are useful in many signal processing applications to clean up the data,
# remove noise, or focus on specific frequency ranges of interest.
#
# 1. Lowpass Filter:
#    - Passes frequencies below a given cutoff and attenuates higher ones.
#    - Useful to remove high-frequency noise or smooth out fast transitions.
#    - Use when: you want to reduce fast noise or spikes.
#
# 2. Highpass Filter:
#    - Passes frequencies above a given cutoff and attenuates lower ones.
#    - Useful to eliminate low-frequency drift, baseline wander, or DC offset.
#    - Use when: you want to remove slow-varying trends or baseline drift.
#
# 3. Notch Filter:
#    - Suppresses a narrow band of frequencies centered around a target frequency.
#    - Useful to eliminate known interferences (e.g., power line noise at 50/60 Hz).
#    - Use when: you want to remove a single frequency from the waveform.
#
# All filter methods rely on SciPy's signal processing library and are intended for
# 1D numpy arrays representing waveform data.
########################################################################################

class Filters:
    def __init__(self, sampling_rate: float):
        """
        :param sampling_rate: Sampling rate of the waveform in Hz (e.g., 5e9 for 5 GS/s).
        """
        self.sampling_rate = sampling_rate
        self.nyquist = 0.5 * self.sampling_rate

    # def lowpass(self, waveform: np.ndarray, cutoff_freq: float, order: int = 4) -> np.ndarray:
    #     """
    #     :param waveform: The waveform data to filter.
    #     :param cutoff_freq: The cutoff frequency (in Hz) for the lowpass filter.
    #     :param order: The order of the filter. Higher orders result in a steeper roll-off.
    #     :return: The filtered waveform.
    #     """
    #     normal_cutoff = cutoff_freq / self.nyquist
    #     b, a = butter(order, normal_cutoff, btype="low")
    #     return filtfilt(b, a, waveform)
    
    @staticmethod
    def lowpass_filter(wf, fs, cutoff_hz=200e6, order=4): # to be merged with the other lowpass function
        """
        Apply a low-pass Butterworth filter to reduce high-frequency noise.
        
        Parameters:
        - wf: waveform array
        - fs: sampling frequency in Hz
        - cutoff_hz: cutoff frequency (default: 200 MHz)
        - order: filter order (default: 4)
        
        Returns:
        - filtered waveform array
        """
        nyq = fs / 2.0
        normal_cutoff = cutoff_hz / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        return filtfilt(b, a, wf)

    def highpass(self, waveform: np.ndarray, cutoff_freq: float, order: int = 4) -> np.ndarray:
        """
        :param waveform: The waveform data to filter.
        :param cutoff_freq: The cutoff frequency (in Hz) for the highpass filter.
        :param order: The order of the filter. Higher orders result in a steeper roll-off.
        :return: The filtered waveform.
        """
        normal_cutoff = cutoff_freq / self.nyquist
        b, a = butter(order, normal_cutoff, btype="high")
        return filtfilt(b, a, waveform)

    def notch(self, waveform: np.ndarray, notch_freq: float, q_factor: float = 30) -> np.ndarray:
        """
        :param waveform: The waveform data to filter.
        :param notch_freq: The frequency (in Hz) to be attenuated.
        :param q_factor: The quality factor (Q) determines the width of the notch. Higher values result in a narrower notch.
        :return: The filtered waveform.
        """
        w0 = notch_freq / self.nyquist
        b, a = iirnotch(w0, q_factor)
        return filtfilt(b, a, waveform)
