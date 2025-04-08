# Waveform Analysis and Threshold Scanning Tools

This repository provides two key classes to process and analyze waveform data:

1. **Filter Class:**  
   Offers a set of filtering functions (lowpass, highpass, and notch) designed to clean up and smooth out your waveform data by removing noise or unwanted frequency components.

2. **ScanThreshold Class:**  
   Contains functionality to perform threshold-based scanning on waveforms. It detects threshold-crossings in a signal and scans through a range of thresholds to evaluate how many crossings occur at each level. This is particularly useful for optimizing detection thresholds and studying noise characteristics.

---

## Filter Class

The `Filter` class is meant for processing 1D numpy arrays representing waveform data. It uses SciPy's signal processing functions to apply three types of filters:

### 1. Lowpass Filter
- **Purpose:**  
  Passes frequencies below a specified cutoff while attenuating higher frequencies.
- **Usage:**  
  Ideal when you want to reduce high-frequency noise or smooth out fast signal transitions.
- **Method:**  
  - `lowpass(waveform, cutoff_freq, order=4)`  
    *Parameters:*  
    - `waveform`: Input waveform data.
    - `cutoff_freq`: Cutoff frequency in Hz.
    - `order`: Filter order (default is 4), with higher orders resulting in a steeper roll-off.
  - **Operation:**  
    It normalizes the cutoff frequency with respect to the Nyquist frequency (half the sampling rate) and applies a Butterworth lowpass filter using forward-backward filtering (`filtfilt`) to avoid phase distortion.

### 2. Highpass Filter
- **Purpose:**  
  Passes frequencies above a specified cutoff, attenuating lower frequencies.
- **Usage:**  
  Useful for eliminating low-frequency drift, baseline wander, or DC offsets.
- **Method:**  
  - `highpass(waveform, cutoff_freq, order=4)`  
    *Parameters:*  
    - `waveform`: Input waveform data.
    - `cutoff_freq`: Cutoff frequency in Hz.
    - `order`: Filter order (default is 4).
  - **Operation:**  
    Similar to the lowpass filter, it normalizes the cutoff frequency and applies a Butterworth highpass filter with forward-backward filtering.

### 3. Notch Filter
- **Purpose:**  
  Attenuates a narrow band of frequencies centered on a target frequency.
- **Usage:**  
  Ideal for removing known interference, such as power line noise (50/60 Hz), from the signal.
- **Method:**  
  - `notch(waveform, notch_freq, q_factor=30)`  
    *Parameters:*  
    - `waveform`: Input waveform data.
    - `notch_freq`: The frequency to be attenuated (in Hz).
    - `q_factor`: Quality factor determining the notch's width; higher values mean a narrower notch.
  - **Operation:**  
    It calculates the normalized frequency using the Nyquist frequency and applies an IIR notch filter using `iirnotch` from SciPy followed by `filtfilt` for zero-phase filtering.

---

## ScanThreshold Class

The `ScanThreshold` class is designed to analyze waveforms by detecting threshold crossings and scanning a range of thresholds. It is particularly useful in setups where you wish to evaluate the hit rate at various detection thresholds.

### Key Features

#### 1. Transition Detection
- **Method:**  
  - `get_transitions(waveform_mv, threshold, sign)`
- **Purpose:**  
  Detects the indices (sample numbers) where the waveform crosses a specified threshold.
- **How It Works:**  
  - **Threshold Crossing Logic:**  
    The method monitors the waveform (given in millivolts) and waits for the signal to first drop (or rise) to 50% of the threshold value. This "arming" phase avoids false or multiple detections from a single pulse.
  - **Sign Parameter:**  
    - `sign = +1`: Detects upward crossings (e.g., a negative pulse rising above the threshold).
    - `sign = -1`: Detects downward crossings (e.g., a positive pulse falling below the threshold).
  - **Operation:**  
    Once the signal is "armed" (i.e., after reaching 50% of the threshold), the first crossing over the full threshold is detected, recorded, and the mechanism resets until the signal drops back below the 50% level. This ensures that only one transition per pulse is recorded.

#### 2. Threshold Scanning
- **Method:**  
  - `scan_thresholds(waveforms_mv, sign, range_min, range_max, step, output_file)`
- **Purpose:**  
  Iterates over a range of threshold values and counts the number of threshold crossings at each level across multiple waveforms.
- **Usage:**  
  - **Inputs:**
    - `waveforms_mv`: A list or array of waveform data in mV.
    - `sign`: Direction for threshold crossing (either +1 or -1).
    - `range_min` and `range_max`: Define the minimum and maximum threshold values to scan (in mV).
    - `step`: The step increment between successive threshold levels.
    - `output_file`: The file in which to save the resulting threshold counts.
  - **Operation:**  
    For each threshold value in the defined range, the method iterates through all provided waveforms, counts the number of detected transitions using `get_transitions`, and then writes the threshold and corresponding counts to a text file. This allows the user to analyze how the hit rate varies with threshold level, which is important for calibration and noise evaluation.

---

## How These Classes Are Useful

- **Filtering:**  
  By applying lowpass, highpass, or notch filters, you can pre-process your waveform data to remove noise, drift, or specific interference frequencies. This improves the quality of your signal analysis and helps in obtaining cleaner, more reliable waveforms.

- **Threshold Scanning:**  
  The ability to scan thresholds and count transitions is crucial in applications such as SiPM signal analysis, where determining the optimal threshold for event detection is necessary. It helps in understanding the relationship between the chosen threshold (in mV) and the resulting hit rate (Hz), leading to improved detector calibration and performance analysis.

---

This README provides a detailed explanation of the two main classes included in the repository. They can be integrated into larger analysis pipelines or used as standalone tools for specific waveform processing tasks.
