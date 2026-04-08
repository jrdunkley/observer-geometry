# Bell / Common-Gluing Failure Beyond Correlators

Reproduce from a clean installed environment:

```bash
pip install -e .
python -m examples.bell_common_gluing.run_all
```

Then inspect:

- [summary.json](outputs/summary.json)
- [audit.json](outputs/audit.json)

## What Is Shown

There exists an explicit Gaussian family where:

- the normalized correlator summary is unchanged
- the arcsine-CHSH correlator obstruction is absent
- common Gaussian gluing still fails at the full law level

The second obstruction in this release object is repeated-marginal variance
inconsistency. The example therefore demonstrates a strict gap between a
correlator summary and the full law-level Gaussian compatibility problem.

The open variance-only neighborhood certified in the validator is the `5 x 5`
grid on:

- `delta in [0.18, 0.26]`
- `rho in [0.56, 0.66]`

## What Is Not Claimed

- this is not a complete Bell theory
- this is not a non-Gaussian locality result
- this does not claim more than the explicit Gaussian common-gluing separation
  certified by the scripts

## Why This Is Difficult To Dismiss

- the compatible and variance-only samples have exactly the same normalized
  correlator matrix
- the compatible and variance-only samples therefore cannot be separated at
  correlator level
- incompatibility in the variance-only sample is certified by a direct linear
  completion residual
- compatibility in the variance-consistent slice is cross-checked by a PSD
  feasibility search and by re-observation through `nomogeo`

## Why This Is Not Just A Restatement Of CHSH

- the compatible sample `(delta, rho) = (0.0, 0.6)` and the variance-only
  sample `(delta, rho) = (0.22, 0.6)` have literally the same correlator matrix
- the arcsine margin is the same positive value for both
- the incompatibility appears only after keeping the full law-level marginal
  data instead of collapsing to correlators

## Public Outputs

- headline figure:
  [phase_diagram.svg](outputs/phase_diagram.svg)
- residual figure:
  [variance_slice.svg](outputs/variance_slice.svg)
- summary table:
  [comparison_table.csv](outputs/comparison_table.csv)

See also:

- [CLAIMS.md](CLAIMS.md)
- [LLM_AUDIT_BRIEF.md](LLM_AUDIT_BRIEF.md)
- [RESULT.md](RESULT.md)

