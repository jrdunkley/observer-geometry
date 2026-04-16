# 0.3.2 Implementation Record

Date: 2026-04-11

Status: implementation completed; ready for user-directed QA.

## Implemented

- `rank_k_covariance_perturbation`
  - exact coordinate-split covariance/Fisher Woodbury diagnostic;
  - reports the difference of two rank-k PSD terms and the rank bound `2k`;
  - keeps `rank_one_covariance_perturbation` backward-compatible as the
    one-vector special case.

- `declared_ladder_dimension_cost_intervals`
  - exact finite declared-ladder phase diagram for scores `s_j - c d_j`;
  - reports pairwise crossings and winning intervals;
  - does not claim global observer optimization.

- `general_graph_frontier_hessian`
  - exact declared-observer graph-chart gradient and Hessian;
  - recovers `exact_branch_hessian` in the off-block-zero sector;
  - handles stationary non-exact observers without weakening
    `exact_branch_hessian`.

- `declared_frontier_local_certificate`
  - sufficient local max/min certificate using stationarity residual, Hessian
    margin, chart radius, and conservative `L_cert`;
  - reports vacuity explicitly.

- `affine_hidden_branch_reversal`
  - exact finite branch flip diagnostic inside the affine-hidden Gaussian-fibre
    sector.

- `guarded_fibre_dominance`
  - diagnostic tuple with centered fibre norm, centered variational norm, and
    a ratio only when the denominator floor is met.

## Files Touched

- `src/nomogeo/types.py`
- `src/nomogeo/perturbation.py`
- `src/nomogeo/frontier.py`
- `src/nomogeo/affine.py`
- `src/nomogeo/__init__.py`
- `tests/test_perturbation.py`
- `tests/test_frontier.py`
- `tests/test_affine.py`
- `README.md`
- `docs/theorem_map.md`
- `docs/validation_note.md`

## Validation

Commands run:

- `python -m pytest tests\test_perturbation.py tests\test_frontier.py tests\test_affine.py -q`
  - result: `27 passed`
- `python -m pytest -q`
  - result: `109 passed`
- `python tools\install_surface_smoke.py`
  - result: imported from workspace and completed the fixed-observer smoke
- `python tools\validation_sweep.py`
  - result: completed; max energy-split residual `2.84e-13`, max tower-law
    residual `3.28e-14`, max hidden clock residual `3.55e-15`
- `python audit\0_3_2_general_graph_hessian_check.py`
- `python audit\0_3_2_graph_hessian_invariance_check.py`
- `python audit\0_3_2_declared_frontier_certificate_check.py`
- `python audit\0_3_2_shelved_feature_recovery_check.py`
- `python audit\0_3_2_certificate_brutality_check.py`
- `python audit\0_3_2_non_exact_stationary_frontier_check.py`
- `python audit\0_3_2_near_branch_variation_check.py`

Key validation values remained aligned with the research record:

- exact-branch reduction: `3.55e-15`;
- non-exact graph-Hessian finite-difference relative Frobenius error:
  `3.35e-07`;
- arbitrary-frame gradient finite-difference error: `9.52e-09`;
- bilinear basis-change invariance error: `5.68e-13`;
- symmetric-pair Hessians at `e = 1, sqrt(3), 2`: `-16, 0, 8`;
- projector-bound brutality violations with 2 percent slack: `0`;
- rank-k recovery residuals in the shelved-feature probe stayed below
  `2.4e-16`.

## Scope Boundary

The implementation keeps the 0.3.2 theory boundary:

- exact Gaussian and special affine-hidden law sectors remain separate from
  local quadratic graph-frontier geometry;
- `exact_branch_hessian` remains strict;
- local certificates are sufficient and may be vacuous;
- declared-ladder intervals rank only supplied candidates;
- no generic non-Gaussian branch probabilities, global Grassmannian optimizer,
  or probability-support/restart engine was added.
