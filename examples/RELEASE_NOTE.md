# Demonstration Release Note

The `nomogeo` kernel is treated as fixed for this release. These examples sit
on top of that validated kernel and are intended as the first public
demonstrations.

## Included

- `entanglement_hidden_load`
  - establishes an exact executable Gaussian identity:
    `tau = 2 I(A:B)` for the stated families
  - does not claim a universal entanglement theorem beyond the stated Gaussian
    setup

- `bell_common_gluing`
  - establishes a strict gap between correlator-level compatibility and full-law
    Gaussian compatibility
  - does not claim a complete Bell theory outside the Gaussian common-gluing
    scope used here

- `arrow_rank_deficiency`
  - establishes a synthetic observer-geometry worked example outside physics
  - does not claim an empirical theory of voting systems or real elections

- `graph_frontier_declared_certificate`
  - establishes the 0.3.2 local quadratic breakpoint where a declared observer
    is stationary without being an exact branch
  - does not claim global observer optimization, probability-support events, or
    full-law branch probabilities

## Release Standard

Each demonstration now contains:

- disciplined claim language
- exact reproduction commands
- main computation and validation scripts
- frozen outputs
- machine-readable `summary.json`
- machine-readable `audit.json`
- an LLM audit brief for external model inspection

## Kernel Amendment Included

This release also promotes the fixed-ceiling inverse theorem already latent in
the hidden-load layer:

- no global inverse of `(H, C) -> Phi_C(H)` is claimed
- the exact ceiling-conditioned inverse theorem is now explicit through
  `inverse_visible_class(T, Lambda, ...)`

See:

- [docs/inverse_theorem.md](../docs/inverse_theorem.md)
- [docs/theorem_map.md](../docs/theorem_map.md)
- [docs/validation_note.md](../docs/validation_note.md)

