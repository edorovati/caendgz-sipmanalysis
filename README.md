# Waveform Analysis with Threshold and Filters

This project is a waveform analysis tool that allows you to visualize, filter, and analyze signal data in NumPy format. By using different filtering options (lowpass, highpass, notch) and threshold scanning, the tool produces graphical outputs and saves the results into text files (e.g., "scan_results.txt").

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Waveform Visualization](#waveform-visualization)
  - [Filter Application](#filter-application)
  - [Threshold Scanning](#threshold-scanning)
- [Examples](#examples)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- **Flexible Visualization:** Displays waveforms in either mV or ADC scale.
- **Custom Filtering:** Applies lowpass, highpass, and notch filters to clean up signals.
- **Threshold Scanning:** Automatically analyzes the hit rate (Hz) versus different thresholds.
- **Argument-based Configuration:** Customize the analysis by specifying the input file, scale type, number of waveforms to display, and filtering settings.

---

## Requirements

- **Python 3.6+**
- **Python Libraries:**
  - `numpy`
  - `argparse`
  - `colorama`
  - *(Make sure the utility files in `utils` are present: `filters.py`, `plotter.py`, `scanthr.py`)*

---

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/your-username/waveform-analysis.git
   cd waveform-analysis
