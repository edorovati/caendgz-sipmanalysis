# DCR – Dark count rate of Silicon PhotonMultiplier

Code for DCR analysis based on counts of peaks in waveforms and subsequently theshold scan for determining primarly the SPAD amplitude, and then DCR and CT probability.

## Directory Structure

```bash
analysis/
├── dcr.sh                  # Main orchestrator script
├── code/                   # Folder which contains codes
│   ├── get_transitions.py  # Step 1: extract transition thresholds
│   ├── dcr.cpp             # Step 2: ROOT macro to build TTree
│   └── dcr_plot.cpp        # Step 3: ROOT macro to plot transitions
├── Data/                   # Input data -> npz files
└── Scan/                   # Temporary scan output directory
```

>  Folders `Data/` must exist before running `dcr.sh`. This folder should contains npz files from digitizer.

---

## Usage

```bash
bash dcr.sh <scan?> <produce_root?> <plot?>
```

### Command-Line Arguments

| Arg    | Values       | Description                                        |
|--------|--------------|----------------------------------------------------|
| `ans1` | `yes`/`no`   | Perform threshold scan (executes `get_transitions.py`) |
| `ans2` | `yes`/`no`   | Produce ROOT file (runs `produce_tree.C`)          |
| `ans3` | `yes`/`no`   | Generate plots (runs `plot.C`)                     |

### Examples

```bash
# Run all steps
bash dcr.sh yes yes yes

# Only scan and build tree
bash dcr.sh yes yes no

# Skip scan, build tree & plot
bash dcr.sh no yes yes
```

---
