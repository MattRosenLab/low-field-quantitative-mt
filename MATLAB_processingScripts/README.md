# MATLAB_processingScripts

MATLAB functions and scripts for processing Z-spectroscopy data and fitting semisolid MT models to quantify spin parameters from data acquired with the sequences in `TNMR_sequenceFiles/`.

| Script | Input sequence | Description |
|---|---|---|
| `ZspecProc.m` | `ZSpecMultB1vals` | Loads Z-spectroscopy data, computes Z-spectra normalized to the far-offset reference, and displays MTR asymmetry curves as a function of offset frequency and B₁ amplitude. |
| `MTfit_multi_prepInputs.m` | `ZSpecMultB1vals`, `2freqZSpecMultB1vals_*` | Script (run section by section) that prepares the `dataset_list` and `param_defs` input variables for `MTfit_multi()` from workspace data. Supports 1-sided Z-spec, 2-sided (alternating frequency) Z-spec, and selective inversion-recovery datasets. |
| `MTfit_multi.m` | — | Performs simultaneous nonlinear least-squares fitting of one or more Z-spectroscopy and/or inversion-recovery datasets to a semisolid MT model. Supports Lorentzian, Gaussian, super-Lorentzian, and Kubo-Tomita semisolid lineshapes. |

## `fittingModels/`
MT model functions and semisolid lineshape functions called by `MTfit_multi.m`. See the [README](fittingModels/README.md) in that subfolder for details.

## `TNMRfileloadScripts/`
Helper functions for parsing Tecmag `.tnt` data files. Used internally by the processing scripts above; these do not need to be called directly. See the [README](TNMRfileloadScripts/README.md) in that subfolder for details.
