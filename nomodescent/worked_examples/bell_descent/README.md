# Bell As A Descent / Compatibility Problem

Run:

```bash
python -m worked_examples.bell_descent.run_main
```

This worked example expresses the Bell-square family as a `ProblemSpec` and
asks whether the four pairwise visible laws admit a common Gaussian descent.

What it shows:

- the compatible sample is admitted by the engine through the audited affine
  PSD search
- the variance-only sample fails by an exact repeated-marginal / linear
  inconsistency certificate
- the correlator-only sample remains compatible at the linear level but fails
  the audited PSD search
- the compatible and variance-only samples have identical normalized correlator
  matrices
