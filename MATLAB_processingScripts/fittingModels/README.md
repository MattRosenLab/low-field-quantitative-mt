# fittingModels

MT model functions called by `MTfit_multi.m` in the parent folder, and semisolid pool lineshape functions used by the MT models.

| Script | Description |
|---|---|
| `mt_model.m` | Standard two-pool semisolid MT model (Henkelman et al., MRM 1993). Accepts a single saturation offset per evaluation. Supports Lorentzian, Gaussian, super-Lorentzian, and Kubo-Tomita semisolid lineshapes. |
| `mt_model_w_dipolar.m` | Extended two-pool MT model that includes a dipolar relaxation term *T*<sub>d</sub> (Morrison et al., J Magn Reson B 1995). Accepts a single saturation offset per evaluation. Supports the same lineshape options as `mt_model.m`. |
| `mt_model_2freq_offFromWater.m` | Two-pool MT model for two-frequency Z-spectroscopy, where two offsets symmetric about a specified carrier offset `Txoff` are applied simultaneously. Based on Henkelman et al., MRM 1993. Supports the same lineshape options as `mt_model.m`. |
| `biexpMTfitv2.m` | Biexponential recovery model for selective inversion-recovery data containing a semisolid pool (Gochberg & Gore, MRM 2003; Gochberg et al., MRM 1999). Estimates the initial semisolid magnetization from the inversion pulse duration. |

## `lineshapeFunctions/`
Semisolid pool lineshape functions used by the MT models above. See the [README](lineshapeFunctions/README.md) in that subfolder for details.
