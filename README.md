# Waveform Analysis  (TODO update read me)

This repository contains a study of SiPM signals acquired from Caen digitizers. The tool not only applies various filters (lowpass, highpass, and notch) to clean the signals but also visualizes the waveforms and performs a threshold scan. The threshold scan is carried out versus the threshold value in mV, allowing you to analyze the hit rate (Hz) at different threshold levels.

---

## Table of Contents

- [Requirements](#requirements)
- [Usage](#usage)
  - [Terminal Commands](#terminal-commands)
  - [Examples](#examples)
- [Project Structure](#project-structure)

---

## Requirements

- **Python 3.12**
- **Python Libraries:**
  - `numpy`
  - `argparse`
  - `colorama`

---

## Requirements

- **Python 3.12+**
- **Python Libraries:**
  - `numpy`
  - `argparse`
  - `colorama`
  - *(Make sure the utility files in `utils` are present: `filters.py`, `plotter.py`, `scanthr.py`)*

---

## Usage

Below are the terminal commands needed to use the script along with examples. This tool is controlled via command line arguments.

### Terminal Commands

**Visualize Waveforms:**  
Run the script to display waveforms with the specified scale (mV or ADC) and filtering options. Replace the argument values as needed.

```bash
python3 transitions.py --npz <filename.npz> --waveform <mV/ADC> --sign <1/-1> --num_waveforms <number> [--lowpass <freq_Hz>] [--highpass <freq_Hz>] [--notch <freq_Hz>]
```
**Perform Threshold Scanning:**
Run the script to execute a threshold scan over a specified range or using a single threshold value. 

```bash
python3 transitions.py --npz <filename.npz> --scanthr --range <min-max or "singlethr <value>"> --sign <1/-1> [--lowpass <freq_Hz>] [--highpass <freq_Hz>] [--notch <freq_Hz>]
```

### Examples

1. Visualize Waveforms in mV Scale with a Lowpass Filter:

```bash
  python3 transitions.py --npz waveforms_5GS.npz --waveform mV --sign -1 --num_waveforms 10 --lowpass 200e6
```

2. Visualize Waveforms in ADC Scale 
```bash 
  python3 transitions.py --npz waveforms_5GS.npz --waveform ADC --sign 1 --num_waveforms 10
```

3. Perform Threshold Scanning on a Continuous Range (0-15) with a Lowpass Filter:
```bash 
  python3 transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1 --lowpass 200e6
```

4. Perform Threshold Scanning on a Continuous Range (0-15) Without Filtering:

```bash 
  python3 transitions.py --npz waveforms_5GS.npz --scanthr --range 0-15 --sign 1
```

5. Perform Threshold Scanning with a Single Threshold (5 mV) and a Lowpass Filter:

```bash 
  python3 transitions.py --npz waveforms_5GS.npz --scanthr --range "singlethr 5" --sign -1 --lowpass 200e6
```

## Project Structure
Here, depicted a brief schematics of folder: 

```bash 
caendgz-sipmanalysis
├── transitions.py          # Main script for waveform analysis
├── README.md               # Project documentation and guide
└── utils/
    ├── filters.py          # Implementation of filters (lowpass, highpass, notch)
    ├── plotter.py          # Functions for waveform visualization
    └── scanthr.py          # Functions for threshold scanning
```






