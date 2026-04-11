# nomogeo

Agent-facing observer-geometry workspace.

Public docs: https://docs.nomogenetics.com/

Current kernel release: `nomogeo 0.30.0`.

Use this repo as a three-layer stack:

- `nomogeo`
  - exact linear / Gaussian kernel
- `nomodescent`
  - exact observer relation and common-descent layer, plus explicitly audited approximation where enabled
- `evidence`
  - evidence encoding and problem assembly with explicit exact / inferred / ambiguous status

Supporting surfaces:

- https://docs.nomogenetics.com/
  - public documentation for install, API overview, examples, and validation
- [examples/README.md](examples/README.md)
  - public demonstrations and runnable examples
- [tests/micro_case_studies/README.md](tests/micro_case_studies/README.md)
  - tiny "we took paper X and did Y" case studies
- [docs/claim_hierarchy.md](docs/claim_hierarchy.md)
  - epistemic split across exact, audited approximate, synthetic, and micro-real outputs
- [docs/release_scope_0_30.md](docs/release_scope_0_30.md)
  - exact release boundary for the support-stratified observation-field layer
- [LICENSE](LICENSE)
- [CITATIONS.md](CITATIONS.md)

## Exact Domain

This workspace is disciplined around:

- linear observers
- finite-dimensional Gaussian / quadratic visible objects
- exact matrix identities where theorems apply
- explicitly audited approximation only where the exact engine deliberately stops

It is not a generic scientific assistant, generic PDF reader, or unconstrained search system.

## Kernel Surface

`nomogeo` keeps scope narrow:

- exact visible precision `Phi_C(H) = (C H^{-1} C^T)^{-1}`
- canonical lift and hidden projector
- local visible calculus `(V, Q)` and determinant-curvature split
- exact closure-adapted whitening, leakage / visibility scores, leakage-channel
  diagnostics, same-rank observer comparison, and commuting-family observer
  synthesis
- exact fixed-observer chart coordinates `(Phi, R, K)`, chart reconstruction,
  observer-transition law, and fixed-observer current / forcing diagnostics
- support-aware hidden-load parametrisation beneath a ceiling
- fixed-ceiling inverse theorem
- hidden-load transport and determinant clock
- contraction factors for associative hidden composition
- observation-field coordinates `Pi <-> Lambda`, support-stratum transport,
  finite birth/death restarts, kernel Schur-jet event classification,
  local coupled birth extraction, and sampled interval-family diagnostics
- thin Donsker-Varadhan and quotient-side Gaussian layers

Runtime deps stay minimal: `numpy`, `scipy`.

## Working Directories

The workspace has three install roots. Run commands from the correct root.

- repo root
  - `python -m pytest -q`
  - `python -m examples.entanglement_hidden_load.run_all`
  - `python -m examples.bell_common_gluing.run_all`
  - `python -m examples.arrow_rank_deficiency.run_all`
  - `python -m tools.stack_soak`
- `nomodescent/`
  - `python -m pytest`
  - `python -m worked_examples.bell_descent.run_main`
  - `python -m worked_examples.free_gaussian_rg.run_main`
  - `python -m worked_examples.replication_fragility.run_main`
- `evidence/`
  - `python -m pytest`
  - `python -m worked_examples.bell_evidence_encoding.run_main`
  - `python -m worked_examples.replication_protocol_encoding.run_main`
  - `python -m worked_examples.benchmark_blindness_encoding.run_main`
  - `python -m micro_real_bundles.bell_counts_bundle.run_main`
  - `python -m micro_real_bundles.iris_protocol_mismatch.run_main`
  - `python -m micro_real_bundles.leaderboard_benchmark_slice.run_main`

## Quick Kernel Use

```python
import numpy as np
from nomogeo import (
    canonical_lift,
    hidden_load,
    inverse_visible_class,
    kernel_schur_jet_from_coefficients,
    pi_from_hidden_load,
    support_stratum_transport,
    visible_precision,
)

H = np.array([[3.0, 1.0], [1.0, 2.0]])
C = np.array([[1.0, 0.0]])
phi = visible_precision(H, C)
lift = canonical_lift(H, C)

T = np.diag([2.0, 1.0, 0.0])
Lambda = np.diag([0.3, 0.8])
X = inverse_visible_class(T, Lambda, lambda_representation="reduced")
load = hidden_load(T, X)

Pi = pi_from_hidden_load(load.reduced_lambda)
transport = support_stratum_transport(load.reduced_lambda, np.diag([0.2, 0.4]))
jet = kernel_schur_jet_from_coefficients([np.diag([0.0, 1.0]), np.diag([1.0, 0.0])])
```

## Important Boundaries

- The fixed-ceiling inverse theorem is exact only after choosing the ceiling `T`. It does not invert the global map `(H, C) -> Phi_C(H)`.
- If `rank(T) = n`, reduced and ambient hidden-load coordinates can have the same shape. In that case you must pass `lambda_representation="reduced"` or `"ambient"`.
- For long hidden composition, use `hidden_contraction(...)` and `load_from_hidden_contraction(...)`. Raw load coordinates are not the associative object.
- The `0.30.0` observation-field layer is exact but narrow: support-stable transport is reduced-coordinate diagnostics, restart maps require explicit nested support bases, kernel jets control leading small-eigenvalue behaviour only, sampled interval diagnostics certify samples only, and no global field simulator or noncommuting optimiser is exposed.
- outside exact Gaussian law mode, treat nomogeo as exact for supplied local quadratic Hessian or Fisher geometry, not as an exact engine for full non Gaussian laws.

## Verification

```bash
python -m pytest -q
python tools/install_surface_smoke.py
python tools/validation_sweep.py
python tools/stack_soak.py
```

For theorem and validation maps, start with:

- [docs/theorem_map.md](docs/theorem_map.md)
- [docs/release_scope_0_30.md](docs/release_scope_0_30.md)
- [docs/release_scope_0_25.md](docs/release_scope_0_25.md)
- [docs/inverse_theorem.md](docs/inverse_theorem.md)
- [docs/validation_note.md](docs/validation_note.md)
- [docs/stack_validation_note.md](docs/stack_validation_note.md)

