# lineshapeFunctions

Semisolid pool lineshape functions used by the MT model functions in the parent folder. Each function returns the RF saturation rate *R*<sub>rfb</sub> for the semisolid pool given the pool's linewidth parameter and the applied RF field.

| Script | Lineshape | Description |
|---|---|---|
| `RF_superlorentzian.m` | Super-Lorentzian | Super-Lorentzian lineshape computed by numerical integration over orientation angles (see, for example, Bieri & Scheffler, MRM 2006). Appropriate for modeling the semisolid pool in tissue. |
| `RF_superlorentzian_roundTop.m` | Super-Lorentzian (smoothed) | Same as above, with the singularity near the peak (±200 Hz from the semisolid resonance) replaced by a cubic spline interpolation to avoid numerical instability at low offsets. |
| `RF_KuboTomita.m` | Kubo-Tomita | Kubo-Tomita lineshape based on a truncated infinite series (Iino, MRM 1994). Interpolates between Lorentzian (*σ*·*τ* ≪ 1) and Gaussian (*σ*·*τ* ≫ 1) lineshapes depending on the parameter `sigtau`. |
