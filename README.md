# nomogeo

Agent-facing observer-geometry workspace.

Public docs: https://docs.nomogenetics.com/

Current kernel release: `nomogeo 0.30.0`.

Use this repo as a three-layer stack:

- `nomogeo`
  - exact linear, quadratic, Gaussian, and declared special-sector kernel
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
- [docs/use_case_map.md](docs/use_case_map.md)
  - practical routing from research goals to exact APIs, examples, and boundaries
- [docs/agent_requirements.md](docs/agent_requirements.md)
  - expected agent quality bar, workflow patterns, and external-proof gates
- [docs/release_scope_0_30.md](docs/release_scope_0_30.md)
  - exact release boundary for the support-stratified observation-field layer
- [LICENSE](LICENSE)
- [CITATIONS.md](CITATIONS.md)

## Exact Domain

This workspace is disciplined around:

- linear observers
- finite-dimensional Gaussian / quadratic visible objects
- exact matrix identities where theorems apply
- exact special law sectors where the law structure is explicitly supplied
- explicitly audited approximation only where the exact engine deliberately stops

It is not a generic scientific assistant, generic PDF reader, or unconstrained search system.

Outside exact Gaussian law mode, read `nomogeo` as exact for supplied local
quadratic Hessian/Fisher geometry, plus separately proved special sectors such
as affine-hidden Gaussian fibres. It is not an exact engine for arbitrary full
non-Gaussian laws.

## What You Can Do With It

Use the workspace to make hard observer-geometry questions concrete:

- understand a theory, paper, or scientific idea more clearly by encoding its
  claimed observers, hidden variables, local Hessians, event strata, or
  evidence bundles as explicit finite objects
- check whether a claim is really supported by separating exact theorem output,
  audited approximation, and unsupported full-law extrapolation
- find the right approach to a hard problem by comparing quotient geometry,
  closure scores, hidden-load coordinates, branch diagnostics, and residual
  margins before committing to a model
- push research forward by stress-testing exact quadratic claims against
  non-Gaussian pathologies, affine-hidden exact sectors, weighted-family
  frontiers, and support-event boundaries
- build apps and live tools around the small stable kernel: visible precision,
  hidden load, local calculus, interval diagnostics, evidence encoders, and
  deterministic batch wrappers
- apply it to datasets by first turning the dataset into an explicit covariance,
  Hessian/Fisher estimate, weighted symmetric family, evidence bundle, or
  problem assembly with declared provenance
- compare different observers with same-rank score comparisons, leakage and
  visibility scores, intrinsic local-geometry ensembles, declared-ladder
  dimension-cost intervals, exact-branch Hessian diagnostics where the branch
  hypothesis is already satisfied, and declared local graph-frontier
  certificates where the stationarity and margin conditions are supplied
- diagnose and understand failures, including noncommuting closure failure,
  missing residual margins, ill-conditioning, support-stratum transitions,
  probability-support mismatch, and log-determinant branch flips
- build simpler models that still keep the important structure through quotient
  precision, fixed-ceiling hidden-load coordinates, minimal hidden
  realisations, and contraction-factor composition
- track changes, branches, and critical events through local visible calculus,
  support-stratum transport, kernel Schur jets, semisimple event charges,
  affine-hidden staged elimination, affine-hidden branch-reversal diagnostics,
  and weighted-family branch Hessians

The important restriction is that every use case starts from declared finite
objects. If the object is only a prose claim, raw paper, or raw dataset, the
first step is to encode the relevant observer, Hessian/Fisher/covariance,
weighted family, residual bound, or law-sector data explicitly.

For a task-to-surface routing table, see
[docs/use_case_map.md](docs/use_case_map.md).

## Kernel Surface

`nomogeo` keeps scope narrow:

- exact visible precision `Phi_C(H) = (C H^{-1} C^T)^{-1}`
- canonical lift and hidden projector
- local visible calculus `(V, Q)` and determinant-curvature split
- exact closure-adapted whitening, leakage / visibility scores, leakage-channel
  diagnostics, same-rank observer comparison, and commuting-family observer
  synthesis
- simple-spectrum closure certificates for exact common-closure obstruction
- exact fixed-observer chart coordinates `(Phi, R, K)`, chart reconstruction,
  observer-transition law, and fixed-observer current / forcing diagnostics
- intrinsic, ceiling-mediated, and coordinate-split local quadratic ensemble diagnostics
- support-aware hidden-load parametrisation beneath a ceiling
- fixed-ceiling inverse theorem
- hidden-load transport and determinant clock
- contraction factors for associative hidden composition
- observation-field coordinates `Pi <-> Lambda`, support-stratum transport,
  finite birth/death restarts, kernel Schur-jet event classification,
  local coupled birth extraction, and sampled interval-family diagnostics
- thin Donsker-Varadhan and quotient-side Gaussian layers
- theorem-local rank-one and rank-k covariance/Fisher perturbation diagnostics
- residual-margin certificates for bounded branch/observer score residuals
- exact affine-hidden Gaussian-fibre reduction, including the variable-precision fibre-volume term, guarded fibre-dominance diagnostics, and finite branch-reversal checks
- finite weighted-family frontier evaluators, declared-ladder cost intervals, exact-branch Hessian diagnostics, general graph-frontier Hessians, and sufficient declared local certificates, without a noncommuting optimiser

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
    declared_frontier_local_certificate,
    declared_ladder_dimension_cost_intervals,
    exact_branch_hessian,
    general_graph_frontier_hessian,
    hidden_load,
    inverse_visible_class,
    kernel_schur_jet_from_coefficients,
    pi_from_hidden_load,
    variable_precision_affine_hidden_reduction,
    support_stratum_transport,
    visible_precision,
    weighted_family_frontier_scores,
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

affine = variable_precision_affine_hidden_reduction(
    np.array([0.0, 0.0]),
    np.zeros((2, 1)),
    np.array([[[0.5]], [[2.0]]]),
)

family = [np.diag([1.0, 0.0, 3.0])]
B = np.array([[1.0], [0.0], [0.0]])
frontier = weighted_family_frontier_scores(family, B)
branch = exact_branch_hessian(family, B)
graph = general_graph_frontier_hessian(family, B)
certificate = declared_frontier_local_certificate(family, B, mode="max")
ladder = declared_ladder_dimension_cost_intervals(np.array([4.0, 5.8, 6.4]), np.array([1.0, 2.0, 4.0]))
```

## Important Boundaries

- The fixed-ceiling inverse theorem is exact only after choosing the ceiling `T`. It does not invert the global map `(H, C) -> Phi_C(H)`.
- If `rank(T) = n`, reduced and ambient hidden-load coordinates can have the same shape. In that case you must pass `lambda_representation="reduced"` or `"ambient"`.
- For long hidden composition, use `hidden_contraction(...)` and `load_from_hidden_contraction(...)`. Raw load coordinates are not the associative object.
- The `0.30.0` observation-field layer is exact but narrow: support-stable transport is reduced-coordinate diagnostics, restart maps require explicit nested support bases, kernel jets control leading small-eigenvalue behaviour only, sampled interval diagnostics certify samples only, and no global field simulator or noncommuting optimiser is exposed.
- The affine-hidden reducer is an exact special full-law sector with supplied `A`, `J`, and hidden precision `D`; it is not arbitrary non-Gaussian marginalisation.
- Weighted-family frontier APIs evaluate supplied finite quadratic families, declared graph-chart Hessians, and sufficient local certificates; they do not choose observers globally or certify full-law branch probabilities.
- `exact_branch_hessian` remains strict and requires an already-exact branch. Use `general_graph_frontier_hessian` for declared-observer local quadratic Hessians outside the exact-branch sector.
- `declared_frontier_local_certificate` is sufficient and can be vacuous; when the stationarity residual is nonzero it certifies a nearby local optimizer with a displacement bound, not global optimality of the supplied observer.
- Declared-ladder dimension-cost intervals rank only the supplied finite ladder. They are not a Grassmannian optimizer or observer discovery routine.
- Local quadratic ensembles summarize exact samplewise Hessian/Fisher geometry. They do not estimate mixture masses, cumulants, remote wells, or probability-support events.
- Matrix support strata are ranks/kernels of supplied PSD matrix paths, not hard probability supports or atoms.

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
