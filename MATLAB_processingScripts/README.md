# MATLAB_processingScripts

MATLAB functions and scripts for processing Z-spectroscopy data and fitting semisolid MT models to quantify spin parameters from data acquired with the sequences in `TNMR_sequenceFiles/`.

| Script | Description |
|---|---|
| `ZspecProc.m` | Takes an array of spectral integrals computed by the NMR_UI program, reshapes the array according to the offset and B₁ amplitude dimensions, and calculates the Z-spectrum for each B₁ value. |
| `MTfit_multi_prepInputs.m` | Helper script (run section by section) that assembles datasets to be fitted into input structures, pairing each dataset with the appropriate quantitative MT model function and fitting parameters, for use with `MTfit_multi.m`. |
| `MTfit_multi.m` | Takes the inputs prepared by `MTfit_multi_prepInputs.m` and performs simultaneous nonlinear least-squares fitting across all included datasets. Returns the fitted spin parameters, confidence intervals, and plots of the data alongside the fitted models. Supports Lorentzian, Gaussian, super-Lorentzian, and Kubo-Tomita semisolid lineshapes. |

## `fittingModels/`
MT model functions and semisolid lineshape functions called by `MTfit_multi.m`. See the [README](fittingModels/README.md) in that subfolder for details.

## `TNMRfileloadScripts/`
Helper functions for parsing Tecmag `.tnt` data files. Used internally by the processing scripts above; these do not need to be called directly. See the [README](TNMRfileloadScripts/README.md) in that subfolder for details.
