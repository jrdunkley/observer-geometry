# Empirical Dataset Protocol v0

Status: narrow intake contract for the first manual datasets.

This is not a live-feed layer and not a broad evidence schema. It is a disciplined way to accept small static datasets and decide whether they compile into a declared local positive curvature object or a refusal.

## Intake Rule

Every dataset must arrive with two objects:

1. a table with measured values;
2. metadata declaring the measurement model, units, controls, uncertainty/noise convention, and validity range.

Without the metadata, the data are not an attachment object. They are only observations.

## Current Targets

| Dataset Kind | Route | Table Columns | Emitted Object If Admissible | Mandatory Refusals |
| --- | --- | --- | --- | --- |
| `rc_decay_trace` | `observation_equation_fisher` | `time_s`, `voltage_V` | local Fisher precision for `tau` | no time samples, no voltage noise model, nonpositive `tau` convention, zero-sensitivity time grid, residual gate failure, hidden-load request |
| `beer_lambert_calibration` | `measurement_equation_local_linearisation` or `observation_equation_fisher` | `concentration`, `absorbance` | local precision for concentration or calibration slope | no absorbance uncertainty, no path length, missing calibration range, denominator not positive, saturated/nonlinear absorbance regime, hidden-load request |

## RC Decay Admission

Model:

```text
V_i = V0 exp(-t_i / tau) + E_i
E_i ~ N(0, sigma_V^2)
```

The dataset is admissible only after declaring:

- time unit and voltage unit;
- whether the trace is charge or discharge;
- start-time convention;
- known or fitted `V0`;
- noise model, at minimum `sigma_V`;
- whether offset voltage is negligible or fitted;
- residual adequacy gate.

Allowed output:

```text
scalar Fisher precision for tau
```

Refused output:

```text
static storage Hessian
hidden_load
```

## Beer-Lambert Admission

Model:

```text
A = epsilon path_length concentration + E
E ~ N(0, sigma_A^2)
```

or, when `epsilon` and `path_length` are supplied:

```text
c = A / (epsilon path_length)
```

The dataset is admissible only after declaring:

- concentration unit;
- absorbance unit or convention;
- wavelength;
- path length;
- blank correction convention;
- absorbance uncertainty or replicate model;
- calibrated concentration and absorbance range.

Allowed output:

```text
local concentration precision or slope/calibration Fisher precision
```

Refused output:

```text
thermodynamic concentration model
static storage Hessian
hidden_load
```

## Placement

Use the templates in `dataset_templates/`. When actual datasets arrive, place copies under:

```text
operationalise/attachment_layer/incoming_data/
```

The first pass should be static and boring: inspect columns, inspect metadata, fit or refuse, record the route, and only then promote any reusable runner.
