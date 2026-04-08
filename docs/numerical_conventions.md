# Numerical Conventions

## Matrix Class Assumptions

- `visible_precision`, `canonical_lift`, `hidden_projector`, and `local_visible_calculus` require:
  - `H` symmetric positive definite
  - `C` surjective
- `hidden_load` accepts positive semidefinite ceilings `T`, but restricts to `S = Ran(T)`.
- `visible_from_hidden_load` requires an explicit reduced-vs-ambient declaration whenever full-rank ceilings make both matrix shapes coincide.
- `clock` and `transport_hidden_load` enforce the theorem domain `Lambda >= 0` on active support.

## Solver Policy

- SPD solves use Cholesky (`scipy.linalg.cho_factor` / `cho_solve`).
- Support restriction and hidden-load square roots use eigendecomposition.
- The package does not form unstable inverses directly in the core visible geometry.

## Symmetry Policy

Derived symmetric quantities are explicitly symmetrised by

```text
A <- (A + A^T) / 2
```

This is used after numerically stable algebra, not as a substitute for it.

## Tolerances

Central tolerance type:

- [`src/nomogeo/validation.py`](../src/nomogeo/validation.py)

Default values:

- absolute tolerance: `1e-10`
- relative tolerance: `1e-8`

Rank and PSD cutoffs are derived from the current spectrum:

```text
cutoff = max(atol, rtol * max(1, max_abs_spectrum))
```

The active support, hidden rank, and rank-sensitive checks all use that central rule.

## Batch Behaviour

`batch_map` preserves task order exactly.

- `backend="process"` uses `ProcessPoolExecutor` when available.
- In restricted environments where process workers are disallowed, it falls back to deterministic threads rather than silently dropping work.

