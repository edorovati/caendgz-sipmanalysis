# Dark Count Rate (DCR) and Waveform Analysis Toolkit

## Overview
This repository provides a suite of tools for analyzing single-photon detector data and waveform transitions:

- **dcr.cpp** – C++ application that:
  - Reads threshold scan files and a 0 V baseline,
  - Subtracts baseline counts,
  - Computes derivatives to identify peaks,
  - Fits Gaussian curves to obtain single p.e. parameters,
  - Calculates dark count rate (DCR) and cross-talk (CT) ratios,
  - Writes results and graphs into `output.root`.

- **dcr_plot.cpp** – C++ application that:
  - Opens `output.root`, reads the summary TTree,
  - Normalizes DCR to acquisition window,
  - Constructs TGraphErrors for p.e. mean, DCR, and CT,
  - Performs a linear fit on p.e. mean vs. overvoltage,
  - Draws plots with error bands and fit info,
  - Saves updated canvases back to `output.root`.

- **transition.py** – Python script for waveform handling:
  - Loads `.npz` waveform datasets,
  - Applies optional digital filters (lowpass, highpass, notch),
  - Displays waveforms in mV or ADC units,
  - Scans thresholds across a range or single value,
  - Outputs hit-rate vs. threshold scans.


## Usage Examples

### 1. Run Dark Count Rate Analysis
```bash
# Generate output.root with DCR & CT results
./dcr
```

Results saved in `output.root` (folder `Summary/` holds the TTree and graphs).

### 2. Plot Results and Fit
```bash
# Add fit & produce canvases back into output.root
./dcr_plot
```

### 3. Waveform Visualization & Filtering
```bash
# Display first 10 waveforms in mV, negative polarity, lowpass @200 MHz
python3 transition.py \
  --npz data/waveforms_5GS.npz \
  --waveform mV --sign -1 --num_waveforms 10 --lowpass 2e8
```

### 4. Threshold Scanning Mode
```bash
# Scan thresholds from 0 mV to 15 mV and save hit-rate curve
python3 transition.py \
  --npz data/waveforms_5GS.npz \
  --scanthr --range 0-15 --sign 1
```

Output: `[Number]OV.txt` (or similar).

---

