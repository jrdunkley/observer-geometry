# Release Scope 0.25

This note freezes the exact public claim for `nomogeo 0.25.0`.

## Exact Public Scope

The `0.25.0` release covers:

- exact visible precision `Phi_C(H)`
- canonical lift and hidden projector
- local visible calculus `(V, Q)` and determinant-curvature split
- support-aware hidden-load coordinates beneath a fixed ceiling
- fixed-ceiling inverse theorem surfaces
- hidden-load transport, clock, and contraction-factor composition
- exact closure-adapted whitening, leakage / visibility scores, leakage channels, same-rank observer comparison, and commuting-family observer synthesis
- exact fixed-observer `(Phi, R, K)` chart coordinates
- exact reconstruction from fixed-observer chart coordinates
- exact observer-to-observer transition law for fixed observers of the same visible rank
- exact fixed-observer current / forcing identity `Q = Phi J R^{-1} J^T Phi`

## Explicitly Out Of Scope

The `0.25.0` release does not claim:

- noncommuting frontier optimisation as a public observer-design API
- `sum_a D_a^2` as a final approximate observer rule
- tower-refinement APIs beyond the documented invariant-flag criterion
- varying-observer connection geometry
- field-theoretic or spatially varying precision language
- empirical pipelines, plotting layers, notebooks, or application-specific wrappers

## Numerical Boundary

The exact theorems are not weakened by numerical caution, but the implementation has a practical conditioning boundary.

- validated comfortably through the current hardening range
- conditioning warnings are surfaced in metadata rather than hidden
- transition / reconstruction numerics should be treated cautiously once `cond(H)` approaches `1e8`

## Operational Boundary

- the install surface is part of the release boundary
- `tools/install_surface_smoke.py` remains a required release check
- environment contamination that imports a different `nomogeo` copy is treated as a release failure, not a user mistake

## Release Reading Rule

Read `0.25.0` as:

- exact fixed-observer and fixed-observer-derived geometry in the declared linear / Gaussian regime
- explicit refusal to overclaim beyond that regime

Do not read `0.25.0` as:

- the full observer-transition programme
- the varying-observer programme
- the full approximate observer-design programme
