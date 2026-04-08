# Problem Contract

`nomodescent` starts from a structured `ProblemSpec`.

## Core Objects

- `ObserverSpec`
  - explicit linear observer matrix
  - provenance and notes

- `VisibleEvidenceSpec`
  - observer reference
  - evidence kind:
    - `covariance`
    - `precision`
  - exact flag
  - optional tolerance
  - authoritative flag

- `ConstraintSpec`
  - extra structured constraints not already captured by visible evidence
  - exact vs approximate is explicit

- `CeilingSpec`
  - optional fixed ceiling/reference object when a problem needs one

- `AssumptionLedger`
  - keeps exact statements separate from approximate modelling assumptions

- `GoalSpec`
  - declares what the caller is asking:
    - factorisation
    - relation classification
    - common completion
    - tower check
    - minimal refinement

## Exact Versus Approximate Separation

This is the central discipline:

- exact observed matrices go into `VisibleEvidenceSpec`
- modelling choices go into the `AssumptionLedger`
- audited approximate search is never silent; it appears in the result
  classification and the audit report

## Current Expressive Boundary

The `v0.2` problem contract is built for:

- linear observer matrices
- Gaussian / quadratic visible evidence
- finite-dimensional latent spaces

It is not yet a generic schema for arbitrary nonlinear science problems.
