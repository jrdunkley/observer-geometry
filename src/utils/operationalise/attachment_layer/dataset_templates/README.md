# Dataset Templates

These templates are for small manual datasets that will stress the measurement and inference routes.

Use them as shapes, not as data. Keep the CSV columns exactly as named unless a later intake note declares a change.

## RC Decay

Files:

- `rc_decay_trace_template.csv`
- `rc_decay_metadata_template.json`

Minimum useful dataset:

```text
time_s, voltage_V
```

A good first trace has at least 8 to 12 samples after `t = 0`, spans a substantial voltage drop, and includes a declared voltage noise or sensor resolution.

## Beer-Lambert

Files:

- `beer_lambert_calibration_template.csv`
- `beer_lambert_metadata_template.json`

Minimum useful dataset:

```text
concentration, absorbance
```

A good first calibration has several concentration levels, includes a blank-correction convention, and stays in the linear absorbance range.

## Refusal Rule

If a dataset arrives without metadata, the right action is refusal with a request for the missing declarations. Do not infer units, noise, wavelength, path length, or start-time convention from the numbers alone.
