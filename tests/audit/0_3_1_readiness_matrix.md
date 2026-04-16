# 0.3.1 Readiness Matrix

Date: 2026-04-11

This began as a planning artifact and now records the implemented narrow
0.3.1 batch. Further module expansion beyond the items listed here still
requires a new checkpoint.

## Executive Judgement

The 0.3.1 technical note is coherent and materially sharpens the programme. It
does not license a generic non-Gaussian branch engine. It does license a typed
extension plan:

1. exact intrinsic quotient geometry on `(H,C)` through `Phi_C(H)`;
2. exact ceiling-mediated hidden-load/clock geometry only after explicit
   ceiling data are supplied;
3. exact special full-law sectors, especially fixed- and variable-precision
   affine-hidden fibres;
4. certified higher-order/residual layers that carry explicit residual,
   cumulant, support, or tail data.

The planning stress tests support the note's central split. The most important
new implementation opportunity is the variable-precision affine-hidden reducer,
because it is exact, small, and exposes an actual branch-verdict correction:
the `0.5 log det D(v)` fibre-volume term. The highest-risk opportunity is the
weighted-family noncommuting frontier, because objective semantics, Hessian
sign conventions, and optimizer boundaries can easily be mislabelled.

## Validation Snapshot

Current non-invasive checks:

- `python -m pytest --collect-only -q`: 100 tests collected.
- `python -m pytest -q`: 100 passed.
- `python audit/0_3_1_planning_stress_tests.py`: overall passed.
- `python audit/weakpoint_research_numerics.py`: passed.
- `python audit/branch_non_gaussian_followup_numerics.py`: passed.
- `python audit/claude_notes_followup_numerics.py`: passed.
- `python tools/adapted_hardening.py`: passed.
- `python tools/connection_flatness_hardening.py`: passed.
- `python tools/validation_sweep.py`: passed.

The planning stress script writes
`audit/outputs/0_3_1_planning_stress_tests.json`.

The concrete implementation plan is in `audit/0_3_1_dev_plan.md`.

## Stress-Test Conclusions

`no_canonical_ceiling`: passed. A hidden shear preserving `C` leaves
`Phi_C(H)` fixed to numerical precision while the coordinate visible block
varies over a wide range. This confirms that hidden load/clock cannot be
advertised as intrinsic to `(H,C)` alone.

`intrinsic_covariance`: passed. Random latent basis changes preserve `Phi`;
visible basis changes follow the expected congruence law. This supports the
intrinsic ensemble API boundary.

`variable_precision_elimination`: passed. Numerical integration over hidden
fibres matches
`A(v)+0.5 log det D(v)-0.5 J(v)^T D(v)^{-1}J(v)` up to centering. A hostile
branch case has flat variational score but exact full-law selection by the
log-det term, so branch verdicts can genuinely change.

`variable_precision_tower`: passed. Random staged eliminations agree with
one-step elimination to machine precision.

`weighted_family_energy`: passed. The integrated identity
`Tr(P M_nu)=S_nu(P)+L_nu(P)` and exact-branch stationarity hold on random
finite families.

`branch_hessian_status`: passed. A minimal diagonal example produces mixed
Hessian status and the graph-chart finite-difference Hessian matches the
contract sign convention: the returned contract operator is the negative
one-half of the second variation.

`residual_certification`: passed. The strict score-margin rule is sharp:
`gap > 2R` is safe, equality is not strict, and `gap < 2R` admits constructive
reversal. The branch-drift bound contains the exact quadratic-linear shift.

`support_event_charge`: passed. Synthetic semisimple support death gives the
expected logarithmic slope `m_sup r / 2`.

`packet_lift`: passed. The Hermitian packet formula matches the Schur
complement, support rank is preserved above the isotropic floor, hidden packet
energy screens monotonically, and the rank-one obstruction is exact.

`empirical_weighted_family`: passed. Iris, leaderboard, and Bell micro-real
bundles satisfy the weighted-family energy split. These tests do not make the
bundles exact closure cases; they only show the finite weighted identity works
on data-shaped families.

## Provisional Surfaces Already Present

These were added before the planning freeze. Do not widen them without
authorization.

| Surface | Status | Release Readiness | Remaining Risk |
|---|---|---:|---|
| `intrinsic_local_quadratic_ensemble` | Provisional module surface | High | Keep hidden-load/clock absent; only `Phi` summaries are intrinsic. |
| `ceiling_mediated_local_quadratic_ensemble` | Provisional module surface | High | Must reject ceilings that do not dominate `Phi`; docs must keep ceiling contract explicit. |
| `coordinate_local_quadratic_ensemble` | Provisional module surface | High | Safe as a coordinate convenience wrapper, not a general observer mode. |
| `rank_one_covariance_perturbation` | Provisional module surface | Medium-high | Correct theorem-local covariance diagnostic; one-channel interpretation must remain white/aligned only. |
| `residual_margin_ordering` | Provisional module surface | High | Only a strict certificate from a supplied residual bound; no residual estimation. |
| `simple_spectrum_closure_certificate` | Provisional module surface | Medium-high | Exact for simple-spectrum anchors; not a noncommuting optimizer and must continue rejecting degenerate anchors. |

Recommended action after authorization: review these surfaces for naming,
docstrings, result metadata, and theorem-map consistency, then either bless them
as 0.3.1 exact theorem-local additions or move them behind an internal/audit
surface. No new theory is needed for these, but release discipline is needed.

## Candidate Implementation Items

| Item | Theory Status | Stress Status | Implementation Status | Recommendation |
|---|---|---|---|---|
| Variable-precision affine-hidden reducer | Exact full-law special sector | Strong | Implemented as `variable_precision_affine_hidden_reduction` | Keep strict affine-hidden fibre scope. |
| Variable-precision tower law | Exact | Strong | Implemented as `staged_affine_hidden_elimination` | Helper only; no broader law claim. |
| Weighted-family energy/frontier evaluators | Exact quadratic frontier | Strong for identities | Implemented as `weighted_family_frontier_scores` | Evaluator only; no optimizer. |
| Exact-branch Hessian diagnostic | Exact at exact branches | Strong for sign/status | Implemented as `exact_branch_hessian` | Requires exact branch; no selector policy. |
| Residual certification ladder | Exact bounds from supplied residual data | Partial strong | Only score-margin wrapper exists | Extend cautiously with supplied Hessian/drift bounds only. |
| Support-event charge | Exact matrix-support stratum result | Synthetic strong | Already present through `semisimple_event_block.clock_log_coefficient` | Docs aligned with `q_sup=m_sup*r/2`; avoid probability-support language. |
| Positive packet lift | Exact Hermitian positive-cone bridge | Strong | Not implemented | Defer until complex/Hermitian dtype policy is explicit. |
| Fibre-cumulant branch diagnostics | Higher-order research layer | Pathology scripts only | Not implemented | Keep as research/audit, not public exact engine. |
| Probability-support events | Extra law data required | Known failure of finite jets | Not implemented | Do not implement from Hessian/kernel jets alone. |
| Laplace entropic quotient | Asymptotic | Numerically strong | Not implemented | Research helper only unless API states asymptotic order and small-noise assumptions. |

## Implemented Order

1. Stabilize the provisional 0.3.1 surfaces already present.
   Check names, docstrings, dataclass fields, theorem-map entries, and
   validation notes. Add no new mathematics.

2. Add a variable-precision affine-hidden reducer.
   Minimal API should accept already-evaluated `A`, `J`, and SPD `D` samples or
   callables over visible points. It should return variational action,
   fibre-volume term, exact visible action up to constants, and optional
   derivative terms only when derivatives are supplied. It must distinguish
   fixed-precision and variable-precision sectors.

3. Add exact staged elimination tests.
   These should prove Schur determinant associativity and completed-square
   associativity over random block orders. This is a release gate for the
   reducer, not a separate broad API.

4. Add weighted-family evaluator primitives.
   Start with `S_nu`, `L_nu`, `M_nu`, captured curvature, and stationarity
   residual for finite weighted families. Do not expose an optimizer initially.

5. Add exact-branch Hessian diagnostic.
   Input should require an exact branch projector or basis. Output should make
   the sign convention explicit: positive contract operator means local maximum
   of the penalized frontier because the second variation is negative. Status
   must be `strict_max`, `degenerate`, `saddle`, or `local_min`.

6. Extend residual certification only for supplied bounds.
   Add Hessian-margin and branch-drift certificates if they can be expressed
   without estimating unknown residuals. Do not infer residuals from data.

7. Align support-event charge documentation.
   The existing `semisimple_event_block` already reports the matrix-support
   charge through `clock_log_coefficient=order*dimension/2` on death-like
   blocks. Do not add probability-support claims.

8. Defer packet lift until dtype policy is settled.
   The exact bridge is real, but exposing complex Hermitian APIs would widen
   the package's numerical contract. It deserves its own small design pass.

9. Keep fibre-cumulant and probability-support layers out of the public exact
   engine.
   They remain research notes and stress-test pathologies until explicit law
   data structures exist.

## Required Gates Before Module Implementation

For each authorized item:

- theorem source and exact domain recorded in `docs/theorem_map.md`;
- API contract written before coding, including non-claims;
- unit tests for theorem identities and malformed inputs;
- edge tests for degeneracy, ill conditioning, and equality thresholds;
- stress numerics in `audit/` updated where applicable;
- README/API notes updated to keep Gaussian, quadratic, exact special-sector,
  and higher-order layers separate;
- full `pytest` green;
- no claims about arbitrary non-Gaussian law exactness unless law-level data
  are explicit inputs.

## Current Weak Points

The residual ladder is only as good as the supplied residual bound. The module
must not pretend to estimate the residual in generic non-Gaussian laws.

Weighted-family frontier implementation could be misread as a noncommuting
closure optimizer. The safe first public surface is evaluators and exact-branch
Hessian diagnostics, not a selector policy.

Branch-Hessian status is a fixed-rank, fixed-support local statement. Support
gap failure and probability-support failure are separate.

Variable-precision fibres are exact but not Gaussian-visible in general. The
right claim is exact hidden Gaussian fibre elimination with possibly
non-Gaussian visible marginal, not generic non-Gaussian closure.

Packet lift is exact but would introduce Hermitian/complex behaviour. That is
not just an extra formula; it changes dtype and positive-cone policy.

Connection hardening still shows conditioning stress: metric residuals grow at
`cond(H) ~ 1e8`, and adapted observer acceptance drops before `1e9`. This is a
numerical warning threshold issue, not a theorem issue.

## Readiness State

Implementation has now covered the narrow authorized batch:

- provisional theorem-local surfaces remain scoped and documented;
- the variable-precision affine-hidden reducer is implemented with strict scope;
- weighted-family evaluator primitives and exact-branch Hessian diagnostics are
  implemented without an optimizer.

The detailed API contracts and release gates for that batch are in
`audit/0_3_1_dev_plan.md`. The rest should remain research or design until a
second planning checkpoint and QA pass.
