# LowFieldQuantitativeMT
Code files for low-field semisolid magnetization transfer spin parameter quantification via NMR spectroscopy

> **Note:** The contents of this repository were originally designed for use on a 6.5 mT low-field MR system.

## Citations
Please cite the following papers if you use the pulse sequences and/or code provided:

## Repository Contents

### `TNMR_sequenceFiles/`
Tecmag TNMR data files (`.tnt`) and associated pulse sequence files (`.tps`) for acquiring Z-spectroscopy data at multiple saturation B₁ amplitudes. Both single-frequency and two-frequency Z-spectroscopy sequences are provided. See the [README](TNMR_sequenceFiles/README.md) in that subfolder for details.

### `MATLAB_processingScripts/`
MATLAB functions and scripts for processing Z-spectroscopy data and fitting semisolid MT models to quantify spin parameters, including the MT model functions and semisolid lineshape functions. See the [README](MATLAB_processingScripts/README.md) in that subfolder for details.
