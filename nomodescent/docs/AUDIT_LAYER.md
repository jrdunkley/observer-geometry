# Audit Layer

Every `DescentResult` carries an `AuditReport`.

## What It Records

- `exact_assumptions`
- `approximate_assumptions`
- `authoritative_inputs`
- `residuals`
- `theorem_layer`
- `falsification_route`

## Why This Matters

The scientific question is not only:

`what conclusion did the engine return?`

It is also:

`which part of that conclusion is theorem-grade, which part is modelling, and what would most quickly break it?`

## Reading Results Correctly

- if the theorem layer says `exact`, the conclusion lives inside the exact
  linear / Gaussian regime implemented by the engine
- if the theorem layer mentions audited approximate PSD search, the result must
  not be reported as theorem-grade
- if the falsification route says supply a positive definite affine completion,
  the problem is still open to direct constructive disproof
