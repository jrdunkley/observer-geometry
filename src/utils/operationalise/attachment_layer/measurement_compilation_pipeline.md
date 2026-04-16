# Measurement Compilation Pipeline v0

Status: design note, not a promoted schema.

This file records the next disciplined compiler shape suggested by the new metrology, process-modeling, and eigenvalue documents.

## Principle

A physical hook is compiled only after both questions are answered:

```text
What finite object can the module see?
What measurement process exposes it?
```

The existing attachment cards mainly answer the first question. The new sources supply the vocabulary for the second.

## Compiler Stages

### 1. Formula Surface

Input: AP/textbook/lab formula.

Allowed outcomes:

- continue to physical attachment test;
- refuse as bare formula;
- route to measurement/inference layer.

Required declarations:

- named measurand;
- physical variables;
- units;
- regime;
- validity range.

### 2. Physical Attachment Test

Question: does the formula already declare a finite physical module object?

Routes:

- exact storage Hessian;
- local quadratic approximation;
- coupled storage with declared ceiling;
- finite modal family;
- not yet attachable.

This is the current card layer.

### 3. Measurement Model Test

Question: if the formula is not a storage object, does it define a measurement model?

Routes:

- `measurement_equation_local_linearisation`: output is a function of uncertain inputs.
- `observation_equation_fisher`: observations arise from parameters plus error.
- `calibration_analysis_function`: calibrated sensor indications are inverted to estimate a measurand.
- `process_control_gate`: check standards, drift, bias, repeatability, reproducibility, gauge R&R.

This layer must preserve refusal states.

### 4. Local Curvature Extraction

Allowed extractions:

- Jacobian-pushed covariance or precision for a measurement equation;
- Fisher information for an observation equation;
- covariance of fitted calibration/analysis function;
- finite modal mass/stiffness matrices after basis and truncation.

Required declarations:

- operating point;
- uncertainty or noise model;
- covariance/correlation assumptions;
- sensitivity/Jacobian or likelihood derivative;
- support/validity range;
- adequacy check.

### 5. Adequacy Gate

Refuse or downgrade when:

- linearisation fails its local check;
- residuals show structure the model does not capture;
- denominator or positive-only inputs are given inappropriate Gaussian treatment;
- extrapolation is attempted outside a calibrated range without theoretical support;
- process control indicates drift, bias, or unstable variability;
- finite-mode truncation is not declared;
- empirical polynomial basis is ill-conditioned or used outside range.

### 6. Module Call Licensing

Initial rule:

```text
storage Hessian cards may call static visible precision;
ceiling cards may call hidden-load only after ceiling dominance is checked;
measurement/inference cards may initially emit evidence/local-curvature objects only;
process-control gates may refuse, but do not call theorem-local kernel surfaces.
```

## Route Matrix

| Source Hook | Best Route | Why |
| --- | --- | --- |
| Beer-Lambert absorbance | `measurement_equation_local_linearisation` or `observation_equation_fisher` | Sensor law becomes concentration precision only after noise and calibration declarations. |
| RC voltage decay | `observation_equation_fisher` | Dynamic exponential is not storage; sampled voltage trace gives local information about `tau`. |
| Thermistor calibration | `calibration_analysis_function` | Instrument indication needs fitted calibration and inverse analysis function with range limits. |
| Load-cell calibration | `calibration_analysis_function` | Force-to-voltage model maps standards to indications; use requires inverse or fitted measurement function. |
| Check standard monitoring | `process_control_gate` | Validates bias/variability stability of the measurement process. |
| Standing string waves | `finite_modal_eigenbasis` | Infinite wave surface becomes finite only after boundary, basis, truncation, and modal metrics. |
| Acoustic resonator | `finite_modal_eigenbasis` plus sensor model | Pressure/velocity field needs finite modal coordinates and microphone/speaker sensor surface. |
| Empirical polynomial calibration | `calibration_analysis_function` with adequacy gate | Useful only inside calibrated range; basis and conditioning matter. |

## Next Build Decision

The next build should not be a generic evidence schema. It should be one narrow proof packet:

```text
measurement_equation_local_linearisation_proofs
```

Minimum contents:

- Beer-Lambert concentration from absorbance;
- cadmium calibration standard concentration by ratio model;
- a negative-control case where first-order linearisation fails or is inadequate.

Current status: implemented as `measurement_equation_linearisation_proofs.py`.

After that, build:

```text
observation_equation_fisher_proofs
```

Minimum contents:

- RC decay `tau`;
- thermistor or load-cell calibration residual gate.

Current status: implemented as `observation_equation_fisher_proofs.py` with RC decay, repeated-current mean, load-cell linear calibration, and two refusal controls.

The finite-mode route is also now implemented as `finite_modal_eigenbasis_proofs.py` with a fixed-end three-mode string, a point-sensor observer-map record, and two refusal controls.

The minimal calibration-analysis route is now implemented as `calibration_analysis_function_proofs.py` with a linear load-cell analysis function, local uncertainty propagation, and two refusal controls.

Only after these proof packets stay stable should a JSON schema be promoted.
