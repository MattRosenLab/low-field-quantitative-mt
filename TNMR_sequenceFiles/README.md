# TNMR_sequenceFiles

Tecmag TNMR data files (`.tnt`) and associated pulse sequence files (`.tps`) for Z-spectroscopy measurements at multiple saturation B₁ amplitudes on a low-field MR system.

Two sequence types are provided:

1. **Single-frequency Z-spectroscopy** — a single saturation offset is applied per scan, sweeping across a range of offset frequencies
2. **Two-frequency Z-spectroscopy** — the carrier remains on-resonance with water, and cosine-modulated RF pulses produce simultaneous saturation at ±offset frequencies relative to the carrier. The filename numbers indicate the range of saturation offsets covered by each data file.

| File(s) | Type | Description |
|---|---|---|
| `ZSpecMultB1vals.tps` | Pulse sequence | Single-frequency Z-spectroscopy at multiple B₁ values |
| `ZSpecMultB1vals.tnt` | Data | Single-frequency Z-spectroscopy data |
| `2freqZSpecMultB1vals_9offsets_ActiveTR.tps` | Pulse sequence | Two-frequency Z-spectroscopy (9 saturation offset values, active TR) |
| `2freqZSpecMultB1vals_4-0kHz.tnt` | Data | Two-frequency Z-spectroscopy data, saturation offset range 0–4 kHz |
| `2freqZSpecMultB1vals_20-5kHz.tnt` | Data | Two-frequency Z-spectroscopy data, saturation offset range 5–20 kHz |
| `2freqZSpecMultB1vals_70-22kHz.tnt` | Data | Two-frequency Z-spectroscopy data, saturation offset range 22–70 kHz |
