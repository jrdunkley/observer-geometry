# Beer-Lambert UV-Vis Probe v0

- dataset: `datasets\edinburgh_data_driven_chemistry_unit09_section3`
- source: https://github.com/Edinburgh-Chemistry-Teaching/Data-driven-chemistry
- route: `observation_equation_fisher`
- status: `admissible_empirical_probe`
- lambda_star_nm: `527.0`

## Admission

Spectra and concentrations compile into a Beer-Lambert calibration/analysis-function probe after declaring a common wavelength, known standards, residual noise model, and unknown-sample inversion.

Refused routes:

- `static storage Hessian`
- `hidden_load`
- `thermodynamic concentration model`
- `molar absorptivity claim without declared path length`

## Calibration

- known calibration observations: 12
- alpha: 0.0053648535
- beta_absorbance_per_mM: 0.96667229
- residual sigma absorbance: 0.0043149359
- r_squared: 0.98551515

## Unknown Sample C

| replicate | absorbance | estimated concentration mM | local precision inverse mM2 |
| ---: | ---: | ---: | ---: |
| 1 | 0.075838707 | 0.072903563 | 44558.378 |
| 2 | 0.077083878 | 0.074191663 | 44372.732 |
| 3 | 0.07730259 | 0.074417916 | 44339.293 |

Mean estimated concentration: `0.073837714 mM`.

## Boundary

This is a sensor-law information object and analysis-function probe. It is not a thermodynamic concentration model, not a storage Hessian, and not a hidden-load object. Molar absorptivity is not claimed because path length is not declared in the downloaded data.
