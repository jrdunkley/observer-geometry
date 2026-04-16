# Measurement And Inference Source Digest v0

Status: first pass over the newly supplied operationalise documents.

This is not a new attachment-card schema. It records what the new documents contribute to the operationalisation pipeline. The central result is that they supply the missing layer between a physical hook and a live use: measurand declaration, measurement model, uncertainty model, adequacy checks, calibration/inversion, and process control.

## Sources Read

| Source | Pages | Operational Role |
| --- | ---: | --- |
| `JCGM_GUM_6_2020.pdf` | 103 | Primary metrology source for developing and using measurement models, including explicit/implicit models, statistical observation equations, local linearisation, Monte Carlo adequacy checks, random variation, shared effects, and cause-and-effect analysis. |
| `NIST.TN.1900.pdf` | 105 | Practical uncertainty guide distinguishing measurement equations from observation equations, with examples for thermistor calibration, molecular weight, cadmium calibration standard, random effects, Monte Carlo propagation, and model adequacy. |
| `mpc.pdf` | 426 | NIST measurement process characterization: check standards, bias, short-term and long-term variability, calibration designs, gauge R&R, uncertainty budgets, propagation of error, and process control. |
| `pmd.pdf` | 60 | NIST process modeling: deterministic function plus random component, model selection/fitting/validation, least squares, weighted least squares, experimental design, residual diagnostics, calibration, and load-cell/ultrasonic/thermal-expansion case studies. |
| `5-Eigenvalue and eigenvector.ipynb` | n/a | Linear algebra route for finite-mode thinking: diagonalisation, matrix powers, matrix exponentials, stability by eigenvalues, Markov steady states, Hermitian/symmetric structure, and similarity transformations. |

## Broad Finding

The director's suggestion is right, but should be sharpened:

```text
AP hook -> physical attachment test -> measurement model -> local/inferential curvature -> adequacy/control gate -> runnable object or refusal
```

The new documents do not mainly add physical formulas. They add the missing compilation discipline for using formulas with data.

## Concepts That Integrate Strongly

### Measurement Equation

A measurement equation declares a measurand as a function of input quantities:

```text
Y = f(X1, ..., Xn)
```

This integrates cleanly with the current attachment discipline. It says: the formula is still not an observer object until the measurand, inputs, units, uncertainty model, and validity range are declared.

Module route:

```text
input covariance/precision + Jacobian of f -> local output covariance/precision
```

This is a local positive curvature route, not a static storage Hessian route.

### Observation Equation

An observation equation declares how observed data arise from a measurand or model parameters:

```text
observation = h(parameter, controls) + error
```

This directly supports the existing RC decay and Beer-Lambert local-curvature proofs. It also explains why those hooks were not failing: they were measurement-inference objects, not storage objects.

Module route:

```text
sensor law + noise model + sensitivities -> Fisher / local precision
```

Allowed route should initially be evidence/local-curvature only. It should not license hidden-load calls without a separately declared ceiling.

### Local Linearisation

GUM-6 gives a clear operational rule:

```text
nonlinear model -> local linear model around declared input estimates -> sensitivity coefficients -> uncertainty propagation
```

This is almost exactly the missing phrase for "declared local positive curvature". The attachment layer should record:

- operating point;
- sensitivity/Jacobian;
- input uncertainty model;
- adequacy check for the linear approximation;
- fallback to Monte Carlo when linearisation is not adequate.

The important refusal is also clear: if the derivative vanishes or the first-order approximation is misleading, a first-order local-curvature object is not admissible without a higher-order or Monte Carlo route.

### Calibration And Analysis Functions

The NIST and GUM examples separate:

- calibration function: known standards -> instrument indication;
- analysis or measurement function: instrument indication -> estimated measurand.

This is a major vocabulary improvement for the attachment layer. It stops us from treating "sensor surface" as a passive list. A sensor surface often needs a calibration/analysis function before it is a valid observer map.

Strong future cases:

- thermistor calibration;
- load-cell calibration;
- Beer-Lambert concentration from absorbance;
- RC time constant from voltage trace;
- resistivity probe process control.

### Process Control

`mpc.pdf` adds the "is this measurement process still valid?" layer:

- check standards;
- bias monitoring;
- short-term repeatability;
- long-term variability;
- drift;
- gauge R&R;
- uncertainty budgets;
- calibration designs.

This should not become a theorem call. It is a gate around the attachment object. A card can be physically and inferentially admissible while its measurement process is out of control.

### Finite Mode And Eigenbasis Thinking

The eigenvalue notebook directly supports the pending finite-mode route:

```text
operator or recurrence -> finite matrix -> eigenbasis -> modal coordinates -> truncated quadratic family
```

This is the right grammar for standing waves, acoustic modes, coupled resonators, and stability of linear dynamic systems. It should be used before any wave card is admitted.

## What Does Not Yet Integrate

### Raw Model Selection

Model selection, AIC/BIC, residual plots, and empirical function choice are operationally essential, but they do not themselves produce a module object. They are adequacy gates and route-selection tools.

### Gauge R&R And Check Standards

These are not physical Hessians or observer objects by themselves. They determine whether a sensor process is stable enough to trust the object it claims to expose.

### Monte Carlo

Monte Carlo is an evaluation route, not an attachment object. It can carry an uncertainty distribution through a model when linearisation is not adequate, but it does not license a theorem-local call unless the output object is converted to an admissible finite covariance/precision or distributional object.

### Empirical Polynomials

Empirical polynomial and rational models are useful calibration functions, but they require a declared validity range and adequacy evidence. They should be refused outside the calibrated range unless there is a theoretical model supporting extrapolation.

## Immediate Recommendation

Do not add a broad evidence schema yet.

Add a narrow `measurement_compilation` workstream with four admissibility routes:

1. `measurement_equation_local_linearisation`
2. `observation_equation_fisher`
3. `calibration_analysis_function`
4. `finite_modal_eigenbasis`

Each route should have explicit refusal conditions.

Current implementation note: the first two proof packets now exist as generated artifacts:

- `measurement_equation_linearisation_proofs.py`
- `observation_equation_fisher_proofs.py`
- `finite_modal_eigenbasis_proofs.py`
- `calibration_analysis_function_proofs.py`

The route-level admissibility contract is now recorded in:

- `admissibility_contract_v0.md`
- `admissibility_contract_v0.json`

The first static empirical intake contract is now recorded in:

- `empirical_dataset_protocol_v0.md`
- `dataset_templates/rc_decay_trace_template.csv`
- `dataset_templates/rc_decay_metadata_template.json`
- `dataset_templates/beer_lambert_calibration_template.csv`
- `dataset_templates/beer_lambert_metadata_template.json`

These templates are deliberately not a live pipeline. They require metadata before data can become a measurement object.

The first acquired online empirical dataset is now recorded in:

- `online_dataset_acquisition_log.md`
- `datasets/edinburgh_data_driven_chemistry_unit09_section3/metadata.json`
- `beer_lambert_uvvis_probe.py`
- `outputs/empirical_datasets/beer_lambert_uvvis_probe.json`
- `outputs/empirical_datasets/beer_lambert_uvvis_probe.md`
- `rc_decay_q7_probe.py`
- `outputs/empirical_datasets/rc_decay_q7_probe.json`
- `outputs/empirical_datasets/rc_decay_q7_probe.md`

This promotes the Beer-Lambert route and the small Q7 RC decay route from proof packets to empirical probes, while still refusing storage-Hessian and hidden-load interpretations.

## Strong Next Proof Objects

1. `beer_lambert_concentration_measurement_equation`
   - `c = A / (epsilon l)`;
   - local sensitivity and uncertainty propagation;
   - refusal if denominator uncertainty is mishandled or if absorbance is outside calibrated range.

2. `rc_decay_tau_observation_equation`
   - voltage observations over time;
   - declared noise model;
   - Fisher precision for `tau`;
   - residual adequacy gate.

3. `thermistor_calibration_analysis_function`
   - polynomial calibration function plus inverse analysis function;
   - validity range;
   - residual/Monte Carlo adequacy gate.

4. `acoustic_modal_observer_with_declared_boundary`
   - now partially implemented as `acoustic_tube_three_mode_pressure_observer`;
   - finite-mode analogue after the now-carded string basis, point-sensor observer, and reference-ceiling hidden-load card;
   - must declare acoustic variables, boundary conditions, finite truncation, and sensor convention;
   - raw impedance, attenuation, reflection, and ultrasound imaging remain refused without further dynamic/sensor declarations.

5. `check_standard_process_gate`
   - not a module object;
   - a process-control wrapper that can refuse otherwise admissible measurement cards if bias/drift/variability are out of control.

## Ontological Update

The earlier card grammar attaches physical formulas to finite module objects. The new documents show that live scientific use requires a second object:

```text
measurement model object
```

It is not a Hessian by itself. It is a declared route by which a sensor process produces a local positive curvature object, uncertainty distribution, or refusal.

That is the correct bridge from AP hooks to a compilation pipeline.
