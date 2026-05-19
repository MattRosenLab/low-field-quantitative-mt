# TNMR_sequenceFiles

Tecmag TNMR data files (`.tnt`) and associated pulse sequence files (`.tps`) for Z-spectroscopy measurements at multiple saturation B₁ amplitudes on a low-field MR system.

Two sequence types are provided:

1. **Single-frequency Z-spectroscopy** — a single saturation offset is applied per scan, sweeping across a range of offset frequencies
2. **Two-frequency Z-spectroscopy** — two saturation offsets symmetric about a specified carrier offset are applied simultaneously per scan (cosine-modulated RF), processed with `mt_model_2freq_offFromWater.m`

| File(s) | Type | Description | Processing script |
|---|---|---|---|
| `ZSpecMultB1vals.tps` | Pulse sequence | Single-frequency Z-spectroscopy at multiple B₁ values | — |
| `ZSpecMultB1vals.tnt` | Data | Single-frequency Z-spectroscopy data | `ZspecProc.m` |
| `2freqZSpecMultB1vals_9offsets_ActiveTR.tps` | Pulse sequence | Two-frequency Z-spectroscopy (9 carrier offset values, active TR) | — |
| `2freqZSpecMultB1vals_4-0kHz.tnt` | Data | Two-frequency Z-spectroscopy data, carrier offset range 0–4 kHz | `MTfit_multi_prepInputs.m` / `MTfit_multi.m` |
| `2freqZSpecMultB1vals_20-5kHz.tnt` | Data | Two-frequency Z-spectroscopy data, carrier offset range 5–20 kHz | `MTfit_multi_prepInputs.m` / `MTfit_multi.m` |
| `2freqZSpecMultB1vals_70-22kHz.tnt` | Data | Two-frequency Z-spectroscopy data, carrier offset range 22–70 kHz | `MTfit_multi_prepInputs.m` / `MTfit_multi.m` |
