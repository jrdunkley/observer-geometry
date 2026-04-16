# 0.3.1 Development Plan

Date: 2026-04-11

Status: implemented for Phases 0-4 after user authorization. Phases 5+ remain
deferred unless separately authorized.

## Release Principle

Every new public surface must be typed by layer:

- `Q`: exact supplied quadratic/Fisher/SPD geometry;
- `E`: exact special full-law sector with explicit law structure;
- `C`: certified correction from supplied residual/higher-order data;
- `R`: research-only, not a public exact engine.

The implementation must never turn a `Q` or `C` result into a generic
non-Gaussian full-law claim.

## Phase 0 - Stabilize Provisional Surfaces

Already present and kept within their existing scope:

- `intrinsic_local_quadratic_ensemble`;
- `ceiling_mediated_local_quadratic_ensemble`;
- `coordinate_local_quadratic_ensemble`;
- `rank_one_covariance_perturbation`;
- `residual_margin_ordering`;
- `simple_spectrum_closure_certificate`.

Action taken:

1. Reviewed result dataclasses for field names and metadata notes.
2. Reviewed public exports.
3. Keep ensemble split exactly as now: intrinsic mode has no hidden load/clock;
   mediated mode requires ceilings.
4. Keep rank-one perturbation interpretation limited to white/aligned
   one-channel and coloured rank-two diagnostics.
5. Keep simple-spectrum certificate rejecting degenerate anchors.
6. Run full tests and stress scripts again.

No new mathematical scope is intended in Phase 0.

## Phase 1 - Variable-Precision Affine-Hidden Reducer

Layer: `E`.

Proposed module: `src/nomogeo/affine.py`.

Implemented public function:

```python
variable_precision_affine_hidden_reduction(A, J, D, tolerances=None)
```

Input contract:

- `A`: scalar or array of shape `sample_shape`;
- `J`: array of shape `sample_shape + (hidden_dim,)`;
- `D`: SPD array of shape `sample_shape + (hidden_dim, hidden_dim)`;
- all samples independent; no interpolation or residual estimation.

Output dataclass:

```python
AffineHiddenReductionResult(
    action=A,
    coupling=J,
    hidden_precision=D,
    hidden_mean=-D^{-1}J,
    variational_action=A - 0.5 J^T D^{-1}J,
    fibre_volume=0.5 log det D,
    visible_action=variational_action + fibre_volume,
    sample_shape=...,
    metadata=...
)
```

Non-claims:

- additive constants independent of visible coordinates are not determined;
- this is not arbitrary non-Gaussian marginalization;
- visible law may be non-Gaussian because `A`, `J`, and `D` may vary with `v`;
- branch selection must include the fibre-volume term.

Release tests:

- scalar hidden fibre numerical integration vs formula;
- multi-hidden random SPD batch formula;
- hostile branch flip: flat variational action, exact log-det selection;
- malformed `A/J/D` shape rejection;
- non-SPD hidden precision rejection;
- fixed-precision sector has constant fibre-volume after centering.

## Phase 2 - Staged Affine-Hidden Elimination

Layer: `E`.

Implemented as a narrow public helper plus tests.

Implemented public function:

```python
staged_affine_hidden_elimination(A, J, D, eliminate, tolerances=None)
```

Input contract:

- `eliminate`: explicit hidden indices to integrate out;
- remaining hidden precision is the Schur complement;
- shifted action includes `0.5 log det D_ee - 0.5 J_e^T D_ee^{-1} J_e`.

Release tests:

- one-step vs all single-coordinate elimination orders;
- block elimination vs repeated scalar eliminations;
- determinant factor composition;
- near-singular but SPD hidden block warnings/rejections under tolerance.

## Phase 3 - Weighted-Family Frontier Evaluators

Layer: `Q`.

Proposed module: extend `src/nomogeo/adapted.py` or add
`src/nomogeo/frontier.py`. Prefer `frontier.py` if Phase 3 and Phase 4 are
implemented together; prefer `adapted.py` if only a small evaluator is added.

Implemented public function:

```python
weighted_family_frontier_scores(family, B, weights=None, mu=0.0, tolerances=None)
```

Input contract:

- `family`: already-whitened symmetric operators, not raw precision
  perturbations unless a separate `H`-aware wrapper is added;
- `B`: orthonormal visible basis;
- `weights`: nonnegative finite weights, normalized internally only for
  reporting if requested;
- `mu >= 0`.

Output dataclass:

```python
WeightedFamilyFrontierResult(
    leakage=L_nu,
    visible_score=S_nu,
    captured_curvature=Tr(P M_nu),
    energy_split_residual=captured_curvature - visible_score - leakage,
    penalized_score=S_nu - mu * L_nu,
    moment_operator=M_nu,
    projector=P,
    weights=...,
    metadata=...
)
```

Non-claims:

- this is not a selector policy;
- top eigenspaces of `M_nu` optimize captured curvature, not necessarily the
  leakage-visibility frontier;
- finite samples certify finite samples unless a quadrature/measure theorem is
  supplied.

Release tests:

- random weighted energy split;
- unweighted agreement with current `closure_scores` convention;
- zero-weight rejection or handling policy;
- data-shaped Iris/leaderboard/Bell energy split;
- malformed weights and non-symmetric family rejection.

## Phase 4 - Exact-Branch Hessian Diagnostic

Layer: `Q`.

Implemented public function:

```python
exact_branch_hessian(family, B, weights=None, mu=0.0, tolerances=None)
```

Input contract:

- `B` is an orthonormal exact branch basis;
- every family member must have off-block norm below tolerance relative to
  `B`; otherwise return/reject as `not_exact_branch`;
- no optimization or continuation is performed.

Output dataclass:

```python
ExactBranchHessianResult(
    hessian_contract=H_mu_U,
    second_variation_operator=-2 * H_mu_U,
    eigenvalues=eig(H_mu_U),
    min_eigenvalue=...,
    nullity=...,
    status="strict_max" | "degenerate" | "saddle" | "local_min",
    off_block_norm=...,
    basis=...,
    metadata=...
)
```

Sign convention:

The returned contract operator is
`H_mu_U=(1+mu) C_U^* C_U - G_U`. Its quadratic form is negative one-half of
the second variation of `F_mu`. Therefore:

- `H_mu_U > 0`: strict local maximum;
- `H_mu_U >= 0` with kernel: degenerate;
- mixed sign: saddle;
- `H_mu_U < 0`: local minimum.

Release tests:

- diagonal strict maximum, local minimum, and saddle witnesses;
- noncommuting internal blocks with exact off-block zero;
- finite-difference graph-chart Hessian comparison;
- rejection of non-exact branch off-blocks;
- equality-threshold degeneracy.

## Phase 5 - Residual Certification Extensions

Layer: `C`.

Existing public function: `residual_margin_ordering`.

Candidate additions:

```python
hessian_margin_certificate(model_hessian_gap, hessian_residual_bound, ...)
branch_drift_certificate(lambda_min, gradient_residual_bound, radius, ...)
```

Input contract:

- all residual bounds are supplied, not estimated;
- strict inequalities must remain strict;
- equality means not certified.

Release tests:

- Weyl-threshold strictness;
- equality boundary;
- constructive score reversal at `gap < 2R`;
- quadratic-linear drift example with exact shift inside theorem bound;
- rejection of negative bounds or nonpositive curvature gaps.

## Phase 6 - Support-Event Charge

Layer: `Q`.

Current state: no new API is needed for the core 0.3.1 charge. The existing
`semisimple_event_block` already exposes:

- `order`;
- `dimension`;
- `pole_coefficient`;
- `clock_log_coefficient = order * dimension / 2` on death-like blocks.

Action after authorization:

1. Keep docs aligned with the 0.3.1 wording `q_sup=m_sup*r/2`.
2. Add at most a small alias or doc example if needed.
3. Do not add probability-support claims.

## Deferred Items

Positive packet lift:

- exact and stress-tested;
- defer until complex/Hermitian dtype policy is explicit;
- likely belongs in a separate `packet.py` bridge module if accepted.

Fibre-cumulant branch diagnostics:

- keep in audit/research scripts;
- no public generic branch selector until law data structures and residual
  bounds exist.

Probability-support events:

- no implementation from Hessian or kernel jets alone;
- require explicit atoms, truncation boundaries, tail masses, or law-family
  state.

Laplace entropic quotient:

- useful asymptotic research helper;
- not an exact public reducer unless the API is explicitly asymptotic and
  small-noise ordered.

## Implemented Authorization Batch

The smallest high-value authorized batch has been implemented:

1. Phase 0 stabilization;
2. Phase 1 variable-precision affine-hidden reducer;
3. Phase 2 tower tests/helper;
4. Phase 3 weighted-family frontier evaluator;
5. Phase 4 exact-branch Hessian diagnostic.

Phase 5 can follow if the user wants residual certification expanded in the
same development cycle. Phase 6 is mostly a documentation alignment task.
