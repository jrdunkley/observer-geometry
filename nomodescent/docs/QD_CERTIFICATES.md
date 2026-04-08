# QD Certificates

This note explains the machine-readable obstruction certificates used by the
narrow Quotient Descent operational core.

## Exact Certificates

- `factorisation_failure`
  - the claimed factor map does not exist within tolerance
- `candidate_not_common_refinement`
  - the proposed refinement does not span at least one observer row space
- `repeated_marginal_inconsistency`
  - the same observed linear functional carries incompatible visible marginal
    variance data
- `linear_constraint_inconsistency`
  - the visible linear constraints admit no common symmetric covariance
- `psd_obstruction`
  - the unique affine completion exists but is not positive definite
- `affine_family_underdetermined`
  - the exact linear constraints define an affine family and the current exact
    engine deliberately refuses to guess inside it

## Audited Approximate Certificates

- `approximate_psd_search_failure`
  - the deterministic low-dimensional affine PSD search found no positive
    definite completion

These are not theorem-grade incompatibility certificates. They are audited
approximate outcomes and must be read together with the audit layer.

## Reading Rule

When a result is negative, read in this order:

1. classification
2. certificate kinds
3. exact vs approximate flag
4. theorem layer
5. falsification route
