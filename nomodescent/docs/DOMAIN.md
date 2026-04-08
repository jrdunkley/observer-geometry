# Supported Domain

`nomodescent v0.2` is intentionally narrow.

## Exact Domain

- linear observers
- Gaussian / quadratic visible evidence
- finite-dimensional matrix problems
- exact factorisation and coarsening logic
- exact quotient tower checks
- exact compatibility classification whenever linear constraints already decide
  the problem
- exact obstruction certificates for:
  - factorisation failure
  - repeated marginal inconsistency
  - linear constraint inconsistency
  - unique affine PSD failure
  - explicit non-nestedness

## Audited Approximate Domain

The only approximate layer included in `v0.2` is:

- deterministic low-dimensional PSD-feasibility search on an exact affine
  covariance family

This layer is always marked as approximate in the result classification and
audit object.

## Deliberately Deferred

- nonlinear observers
- non-Gaussian latent completion
- generic heuristic search over observer spaces
- automated document parsing
- unconstrained optimisation over modelling choices
