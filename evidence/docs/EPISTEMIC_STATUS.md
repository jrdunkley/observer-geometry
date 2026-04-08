# Epistemic Status

Every evidence item carries two separate labels.

## Extraction Mode

- `exact_extraction`
  - the source states the quantity directly
- `encoded_inference`
  - the quantity was encoded from a documented modelling choice
- `open_ambiguity`
  - the source does not determine the quantity cleanly

## Epistemic Status

- `exact`
- `parsed`
- `inferred`
- `conjectural`
- `ambiguous`

Rules enforced in code:

- `exact_extraction` must use epistemic status `exact`
- `open_ambiguity` must use epistemic status `ambiguous`
- `open_ambiguity` items cannot be authoritative for downstream assembly

The package is intentionally hostile to mixing these labels up.
