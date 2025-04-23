## General Classes

This repository provides three key classes:

1. **Filter Class**
   - Offers a set of filtering functions (lowpass, highpass, and notch) designed to clean up and smooth out your waveform data by removing noise or unwanted frequency components.

2. **ScanThreshold Class**
   - Contains functionality to perform threshold-based scanning on waveforms. It detects threshold-crossings in a signal and scans through a range of thresholds to evaluate how many crossings occur at each level. This is particularly useful for optimizing detection thresholds and studying noise characteristics.

3. **Utils Class**
   - Provides utility functions, such as generating standardized output filenames for threshold scan results based on input waveform filenames. This ensures consistent naming conventions across analyses.

---

## Filter Class

The `Filter` class processes 1D numpy arrays representing waveform data using SciPy's signal processing functions. It offers three filter types:

### 1. Lowpass Filter

- **Purpose**: Passes frequencies below a specified cutoff while attenuating higher frequencies.
- **Usage**: Ideal for reducing high-frequency noise or smoothing fast signal transitions.
- **Method**:
  ```python
  lowpass(waveform, cutoff_freq, order=4)
  ```  
  **Parameters**:
  - `waveform`: Input waveform data (1D numpy array).
  - `cutoff_freq`: Cutoff frequency in Hz.
  - `order`: Filter order (default: 4). Higher orders yield a steeper roll-off.

  **Operation**:
  - Normalizes the cutoff frequency with respect to the Nyquist frequency (half the sampling rate).
  - Applies a Butterworth lowpass filter using forward-backward filtering (`filtfilt`) for zero phase distortion.

### 2. Highpass Filter

- **Purpose**: Passes frequencies above a specified cutoff while attenuating lower frequencies.
- **Usage**: Useful for eliminating low-frequency drift, baseline wander, or DC offsets.
- **Method**:
  ```python
  highpass(waveform, cutoff_freq, order=4)
  ```
  **Parameters**:
  - `waveform`: Input waveform data.
  - `cutoff_freq`: Cutoff frequency in Hz.
  - `order`: Filter order (default: 4).

  **Operation**:
  - Normalizes the cutoff frequency.
  - Applies a Butterworth highpass filter with forward-backward filtering.

### 3. Notch Filter

- **Purpose**: Attenuates a narrow band of frequencies centered on a target frequency.
- **Usage**: Ideal for removing known interference (e.g., 50/60 Hz power line noise).
- **Method**:
  ```python
  notch(waveform, notch_freq, q_factor=30)
  ```
  **Parameters**:
  - `waveform`: Input waveform data.
  - `notch_freq`: Frequency to attenuate (Hz).
  - `q_factor`: Quality factor determining notch width (default: 30). Higher values yield narrower notches.

  **Operation**:
  - Computes normalized frequency relative to the Nyquist frequency.
  - Designs an IIR notch filter using SciPy's `iirnotch`.
  - Applies zero-phase filtering via `filtfilt`.

---

## ScanThreshold Class

The `ScanThreshold` class analyzes waveforms by detecting threshold crossings and scanning a range of thresholds to count detections.

### Key Features

1. **Transition Detection**
   - **Method**:
     ```python
     get_transitions(waveform_mv, threshold, sign)
     ```
   - **Purpose**: Detects sample indices where the waveform crosses a specified threshold.
   - **Parameters**:
     - `waveform_mv`: Waveform data in millivolts.
     - `threshold`: Threshold level in millivolts.
     - `sign`: +1 for upward crossings, -1 for downward crossings.
   - **Operation**:
     - Waits for the signal to reach 50% of the threshold to "arm" detection and avoid false triggers.
     - Records the first full crossing of the threshold post-arming.
     - Resets when the signal returns below 50% of the threshold to detect subsequent pulses.

2. **Threshold Scanning**
   - **Method**:
     ```python
     scan_thresholds(
         waveforms_mv, sign,
         range_min, range_max,
         step, output_file
     )
     ```
   - **Purpose**: Iterates over a range of thresholds to count crossings across multiple waveforms.
   - **Parameters**:
     - `waveforms_mv`: List or array of waveform data arrays in millivolts.
     - `sign`: Direction for threshold crossing (+1 or -1).
     - `range_min` / `range_max`: Minimum and maximum thresholds to scan (mV).
     - `step`: Increment between thresholds (mV).
     - `output_file`: Path to save a text file with threshold vs. counts.
   - **Operation**:
     - For each threshold value, counts transitions in all waveforms using `get_transitions`.
     - Writes a two-column text file: threshold value and total count.

---

## Utils Class

The `Utils` class provides helper functions for consistent naming and file management:

### Filename Generation

- **Method**:
  ```python
  generate_output_filename(npz_filename: str) -> str
  ```
- **Purpose**: Extracts the voltage value from an input `.npz` filename, computes its difference from the reference voltage (51.5V), and returns a standardized output filename.
- **Usage**:
  - Input: `waveforms_750.0_52V.npz`
  - Reference voltage: 51.5V
  - Voltage difference: 52.0V - 51.5V = 0.5V
  - Output filename: `0,5OV.txt`
- **Operation**:
  - Uses a regular expression to extract the voltage from the filename.
  - Calculates the difference from 51.5V.
  - Formats the difference: no decimals for integers; replace decimal point with a comma otherwise.

