# QD Core Scope

This note fixes the narrow operational Quotient Descent core for the current
release pass.

## In Scope

The operational core contains five items:

1. `factorisation_test(...)`
2. `common_refinement_test(...)`
3. `common_descent_test(...)`
4. `false_collapse_diagnostic(...)`
5. `minimal_refinement_search(...)`

Plus one convenience classifier:

- `classify_qd_relation(...)`

## Exact Input Domain

- finite-dimensional linear observers
- Gaussian / quadratic visible objects
- matrix-based compatibility and refinement questions
- finite candidate observer families for refinement search

## Exact Output Classes

- `exact_factorisation`
- `equivalent_observers`
- `exact_common_refinement`
- `exact_common_descent`
- `incompatible_by_linear_inconsistency`
- `incompatible_by_psd_obstruction`
- `underdetermined_affine_family`
- `non_nested_observers`
- `false_collapse_detected` when the richer incompatibility is exact

## Audited Approximate Boundary

Approximation may enter only through the already-existing deterministic affine
PSD search in `common_descent_test(...)`.

When that happens, the QD layer must preserve the audited status of outcomes
such as:

- `approximate_common_descent`
- `incompatible_by_approximate_psd_search`

Any diagnostic built on top of those results inherits that audited approximate
status.

## Certificate Types In Scope

- `factorisation_failure`
- `candidate_not_common_refinement`
- `repeated_marginal_inconsistency`
- `linear_constraint_inconsistency`
- `psd_obstruction`
- `affine_family_underdetermined`
- `approximate_psd_search_failure`

## Deferred

- broad search over observer spaces
- nonlinear observer models
- arbitrary symbolic closure machinery
- generic closure or spectral-shadow diagnostics beyond the current exact layer
- any widening beyond the declared linear / Gaussian / finite-family regime
