# Online Dataset Acquisition Log v0

Status: narrow log for externally acquired empirical datasets.

## Acquired And Used

| Dataset | Source | Local Path | Route | Status |
| --- | --- | --- | --- | --- |
| Rhodamine 6G UV-Vis spectra, Unit 09 Section 3 | `https://github.com/Edinburgh-Chemistry-Teaching/Data-driven-chemistry/tree/main/Unit_09/data/Section3` | `datasets/edinburgh_data_driven_chemistry_unit09_section3/` | `observation_equation_fisher` | admitted as `beer_lambert_uvvis_probe.v0` |
| Q7 capacitor voltage decay | `https://github.com/Amey-Thakur/COMPUTATIONAL-METHODS-AND-MODELING-FOR-ENGINEERING-APPLICATIONS` via local file `Q7_CapacitorVoltageDecay.m` | `datasets/Q7_CapacitorVoltageDecay.m` | `observation_equation_fisher` | admitted as `rc_decay_q7_probe.v0` |

Why this dataset is useful:

- it contains spectra for known Rhodamine 6G concentrations and an unknown sample;
- it gives the Beer-Lambert route a real sensor surface rather than a toy scalar formula;
- it forces calibration, wavelength choice, residual noise, and analysis-function inversion to be declared;
- it refuses static storage Hessian, hidden load, thermodynamic concentration modelling, and molar absorptivity claims without a declared path length.

Generated artifacts:

- `outputs/empirical_datasets/beer_lambert_uvvis_probe.json`
- `outputs/empirical_datasets/beer_lambert_uvvis_probe.md`
- `outputs/empirical_datasets/rc_decay_q7_probe.json`
- `outputs/empirical_datasets/rc_decay_q7_probe.md`

## Q7 RC Decay Note

The Q7 data are useful as a small textbook exponential-discharge probe, not as a high-quality lab dataset. It provides nine voltage samples, declares `V0 = 100 V`, and supports a residual-estimated local Fisher precision for `tau`.

Boundaries:

- the late point `3.5 s, 2 V -> 4.0 s, 3 V` is nonmonotone and is recorded as a warning;
- no independent voltage-instrument noise is declared;
- no resistance is declared, so capacitance is not inferred;
- storage-Hessian and hidden-load routes remain refused.

## Searched But Not Used

The local supercapacitor discharge dataset under `datasets/Hnandreas-Supercapacitor-Discharge-Measurements-25-F-50-F-DUT-Sets--1d88756/` is not used for the current pass. It is a constant-current supercapacitor characterization set, not the clean exponential RC decay dataset needed for the current `rc_decay_tau_observation_equation` route.

A small public RC CSV candidate was attempted but resolved to a `404: Not Found` body and was removed.
