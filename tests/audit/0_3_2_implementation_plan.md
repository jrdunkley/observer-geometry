# 0.3.2 Implementation Plan

Date: 2026-04-11

Status: planning artifact. No module source changes are made here. The theory
anchor is `papers/0.3.2_Technical_Note_1.tex`.

## Executive Readiness

0.3.2 is ready for phased implementation, provided the implementation keeps the
same scope discipline as the technical note:

- exact Gaussian law mode stays exact Gaussian law mode;
- general graph-frontier Hessians and local certificates are exact local
  quadratic geometry;
- affine-hidden branch reversal is an exact special full-law sector;
- residual, ensemble, and fibre-dominance helpers are finite declared
  diagnostics, not generic non-Gaussian law engines;
- global noncommuting observer optimization and generic non-Gaussian branch
  probabilities remain blocked.

The strongest implementation win is the general graph-frontier Hessian. It
strictly contains the existing exact-branch Hessian as the off-block-zero sector
and correctly handles stationary non-exact observers where the exact-branch
proxy fails.

## Final Sense Checks

The note `papers/0.3.2_Technical_Note_1.tex` was created and compiled with
`pdflatex -interaction=nonstopmode -halt-on-error`. The generated build
artifacts were removed after the compile check. The log had no real document
errors after cleanup.

Overclaim scan:

- the TeX explicitly says the graph Hessian is not a generic non-Gaussian
  full-law selector;
- exact-branch Hessian is not softened into a tolerance-based approximate
  branch Hessian;
- the local certificate is sufficient and can be vacuous;
- fibre dominance is a guarded diagnostic tuple, not an invariant naked ratio;
- blocked items remain marked blocked.

## Final Numerics

The following stress scripts were rerun after the TeX note was created:

- `python audit\0_3_2_general_graph_hessian_check.py`
- `python audit\0_3_2_graph_hessian_invariance_check.py`
- `python audit\0_3_2_declared_frontier_certificate_check.py`
- `python audit\0_3_2_shelved_feature_recovery_check.py`
- `python audit\0_3_2_certificate_brutality_check.py`
- `python audit\0_3_2_non_exact_stationary_frontier_check.py`
- `python audit\0_3_2_near_branch_variation_check.py`
- `python -m pytest tests\test_frontier.py tests\test_perturbation.py tests\test_affine.py -q`
- `python -m pytest -q`

Key outcomes:

- exact-branch reduction matches current module Hessian: `3.55e-15`;
- non-exact graph-Hessian finite-difference relative Frobenius error:
  `3.35e-07`;
- arbitrary-frame gradient finite-difference error: `9.52e-09`;
- bilinear basis-change invariance error: `5.68e-13`;
- exact-branch arbitrary-frame reduction in the module complement frame:
  `1.42e-14`;
- symmetric non-exact breakpoint gives Hessians `-16, 0, 8` at
  `e=1, sqrt(3), 2`, while the invalid exact-branch proxy stays `-24`;
- certificate brutality check found `0` projector-bound violations with
  2 percent slack over the tested radius grid;
- sharpened certificate constant was about `2.06x` less pessimistic than the
  earlier coarse constant in the tested synthetic cases;
- rank-k perturbation residuals for `k=1,2,3` stayed below `2.4e-16`;
- declared-ladder and affine-hidden branch-reversal finite probes found no
  mismatches.
- targeted module tests for the planned frontier/perturbation/affine surfaces
  passed: `18 passed`;
- the full module suite passed: `100 passed`.

## Implementation Phases

### Phase 1 - Low-Risk Finite Algebra

Files likely touched:

- `src/nomogeo/perturbation.py`
- `src/nomogeo/types.py`
- `src/nomogeo/__init__.py`
- `tests/test_perturbation.py`
- `docs/theorem_map.md`
- `docs/validation_note.md`
- `README.md`

Features:

- add `rank_k_covariance_perturbation(...)`;
- keep `rank_one_covariance_perturbation(...)` as a compatibility wrapper or
  one-vector special case;
- add a result dataclass reporting the full and visible rank-k update terms,
  formula residual, singular values, update rank, and `rank_bound=2*k`;
- keep the coordinate visible split explicit.

Readiness:

- theorem-backed by Woodbury;
- low blast radius;
- strong tests already exist for `k=1`; add `k=2,3` and degenerate-rank
  cases.

Non-claims:

- not a non-Gaussian law selector;
- not a coordinate-free covariance perturbation theorem unless the observer
  contract is separately generalized.

### Phase 2 - Declared Finite Comparisons

Files likely touched:

- `src/nomogeo/frontier.py` or a new finite-comparison helper module;
- `src/nomogeo/types.py`;
- `src/nomogeo/__init__.py`;
- `tests/test_frontier.py`;
- `docs/theorem_map.md`;
- `README.md`.

Features:

- add `declared_ladder_dimension_cost_intervals(scores, dimensions, cost_domain=None)`;
- return all pairwise crossings, candidate winning intervals, never-winning
  candidates, and tie flags;
- add `finite_candidate_residual_margin(...)` only if it genuinely generalizes
  the current `residual_margin_ordering(...)` without confusing the existing
  two-candidate API.

Readiness:

- exact finite affine inequality theorem;
- safe because it only ranks a supplied declared ladder.

Non-claims:

- not global Grassmannian optimization;
- not an observer discovery algorithm.

### Phase 3 - General Graph-Frontier Hessian

Files likely touched:

- `src/nomogeo/frontier.py`;
- `src/nomogeo/types.py`;
- `src/nomogeo/__init__.py`;
- `tests/test_frontier.py`;
- `docs/theorem_map.md`;
- `docs/validation_note.md`;
- `README.md`.

Feature:

- add `general_graph_frontier_hessian(family, B, weights=None, mu=0.0, complement_basis=None, tolerances=None)`.

Required result fields:

- `gradient`: graph-chart first variation in the chosen complement frame;
- `hessian_operator`: polarized Hessian on `Hom(U,W)`;
- `second_variation_operator`: same object with explicit sign convention;
- `eigenvalues`;
- `status`: `strict_local_max`, `strict_local_min`, `saddle`,
  `stationary_degenerate`, or `nonstationary`;
- `stationarity_residual`;
- `off_block_norm`;
- `basis`;
- `complement_basis`;
- `weights`;
- `metadata`.

Implementation contract:

- if no complement basis is supplied, compute the module's default null-space
  complement as `exact_branch_hessian` does;
- if a complement basis is supplied, validate that it is orthonormal and
  orthogonal to `B`;
- document that raw Hessian matrices depend on the complement basis, while the
  bilinear form and spectrum are invariant under orthogonal frame changes;
- do not weaken `exact_branch_hessian`; keep it strict.

Tests:

- exact-branch reduction against `exact_branch_hessian`;
- non-exact finite-difference graph Hessian;
- arbitrary-frame gradient finite differences;
- basis-change bilinear invariance;
- symmetric-pair sign flip at `sqrt(3)`;
- bad basis/complement validation.

Non-claims:

- not a branch optimizer;
- not a probability-support classifier;
- not a full non-Gaussian law result.

### Phase 4 - Declared Local Frontier Certificate

Files likely touched:

- `src/nomogeo/frontier.py`;
- `src/nomogeo/types.py`;
- `src/nomogeo/__init__.py`;
- `tests/test_frontier.py`;
- `docs/theorem_map.md`;
- `docs/validation_note.md`;
- `README.md`.

Feature:

- add `declared_frontier_local_certificate(...)` as a wrapper around the
  general graph Hessian plus a conservative `L_cert(rho)` bound.

Required inputs:

- family;
- declared observer basis `B`;
- weights and `mu`;
- `mode="max"` or `mode="min"`;
- chart radius `rho`, with theorem default requiring `0 < rho <= 1`;
- optional complement basis;
- tolerances.

Required result fields:

- graph-Hessian result or enough embedded fields for inspection;
- `eps`;
- `lambda_margin`;
- `L_cert`;
- `r0`;
- `left_4eps_over_lambda`;
- `displacement_bound_2eps_over_lambda`;
- `certificate_passes`;
- `certificate_kind`: sufficient local max/min certificate or vacuous;
- metadata notes saying the certificate is sufficient and can be vacuous.

Tests:

- symmetric non-exact cases pass/fail exactly as in the probe:
  max passes at `e=1` and `e=1.7`, fails at `sqrt(3)`, min passes at
  `e=1.75` and `e=2`;
- imbalance cases show residual-controlled vacuity;
- random near-branch cases include both nonvacuous and vacuous certificates;
- certificate constants reject invalid `rho` and invalid `mu`.

Non-claims:

- does not certify the original observer is optimal when `eps > 0`; it
  certifies a nearby local optimizer with a displacement bound;
- not global over the Grassmannian;
- not a law-level branch probability.

### Phase 5 - Affine-Hidden Diagnostics

Files likely touched:

- `src/nomogeo/affine.py`;
- `src/nomogeo/types.py`;
- `src/nomogeo/__init__.py`;
- `tests/test_affine.py`;
- `docs/theorem_map.md`;
- `README.md`.

Features:

- add `affine_hidden_branch_reversal(...)` for finite branch arrays of
  variational and fibre terms, or implement it as a small helper around
  `variable_precision_affine_hidden_reduction(...)`;
- add guarded fibre-dominance diagnostics only as a tuple:
  `fibre_centered_norm`, `variational_centered_norm`, optional ratio, ratio
  defined flag, denominator floor, and declared sample measure label.

Readiness:

- exact in the affine-hidden Gaussian-fibre sector;
- diagnostics are useful but should be labelled as diagnostics.

Non-claims:

- no generic marginalization;
- no unguarded fibre-dominance scalar.

### Phase 6 - Documentation, Examples, and Release Gate

Files likely touched:

- `README.md`;
- `docs/theorem_map.md`;
- `docs/validation_note.md`;
- possibly `examples/`.

Documentation updates:

- add `0.3.2_Technical_Note_1.tex` as the source for graph Hessian and local
  certificates;
- update README use-case language to include:
  - compare declared observer ladders under dimension cost;
  - diagnose local frontier stability at declared observers;
  - classify rank-k covariance/Fisher perturbation channels;
  - detect affine-hidden fibre-volume branch reversal in the exact sector;
  - certify finite candidate gaps with supplied residual or sampling bounds.

Validation gate after implementation:

- `python -m pytest -q`;
- targeted tests:
  - `python -m pytest tests/test_frontier.py -q`;
  - `python -m pytest tests/test_perturbation.py -q`;
  - `python -m pytest tests/test_affine.py -q`;
- rerun the audit probes and compare the key numeric tolerances above;
- run README/theorem-map overclaim scans for:
  `generic non-Gaussian`, `global optimizer`, `probability support`,
  `approximate exact-branch`, and `full-law selector`.

## Deferred / Blocked

Do not implement these in 0.3.2:

- generic non-Gaussian branch probabilities;
- automatic higher-order cumulant branch selector;
- global noncommuting Grassmannian optimizer;
- probability support or restart event engine from Hessians alone;
- unbounded-score ensemble concentration for log-det clocks near support
  boundaries;
- coordinate-free affine-hidden stage-sign interpretation beyond a declared
  hidden measure convention.

## Order Of Work

Recommended order:

1. implement rank-k perturbation and tests;
2. implement declared-ladder interval helper and tests;
3. implement general graph-frontier Hessian and tests;
4. implement declared local certificate and tests;
5. implement affine-hidden branch reversal and guarded fibre-dominance
   diagnostics;
6. update public exports, docs, theorem map, README, and examples;
7. run the full validation gate.

This order keeps the low-risk exact finite algebra separate from the broader
frontier Hessian work, and it keeps the certificate layer dependent on the
graph-Hessian API rather than duplicating formulas.

## Final Readiness Judgement

The plan is ready for implementation authorization. The one caution is that
the declared local frontier certificate should be positioned as a sufficient
certificate and expected to be vacuous in many practical near-branch cases. Its
value is not universal acceptance; its value is that when it passes, it passes
for theorem-level reasons.
