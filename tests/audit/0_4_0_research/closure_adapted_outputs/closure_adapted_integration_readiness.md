# Closure-Adapted Observers: Pre-Integration Readiness

Date: 2026-04-09

Scope: hardening and release-gate planning before any promotion into `src/nomogeo`.

## Status

The exact closure-adapted core is numerically ready for module integration.

The revised [Closure_Adapted_Observers.tex](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/Closure_Adapted_Observers.tex) now aligns with the intended module boundary:
- noncommuting design is stated as a leakage-visibility frontier
- `sum_a D_a^2` is stated as a warm start, not a final observer rule
- tower refinement is stated through the quotient-family / invariant-flag criterion

What is ready:
- whitened observer normal form
- exact closure criterion by common invariant subspace
- commuting-family exact solver by aggregate common-mode energy
- bridge-adapted exact observer demo

What is not ready for promotion as public module behaviour:
- noncommuting approximate-design optimisation as a finished API
- `sum_a D_a^2` as a final observer rule
- tower/refinement language that omits the invariant-flag criterion

## Rerun Verification

The existing paper-side harnesses were rerun from the freeze.

Successful reruns:
- `python papers/closure_adapted_audit.py`
- `python -m papers.closure_adapted_use_cases`
- `python -m papers.closure_adapted_implications`

Current exact evidence:
- [closure_adapted_audit_results.json](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_outputs/closure_adapted_audit_results.json)
- [closure_adapted_use_cases.json](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_outputs/closure_adapted_use_cases.json)
- [closure_adapted_implications.json](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_outputs/closure_adapted_implications.json)

Numerically confirmed exact facts:
- normal-form residuals remain at machine scale:
  - max `Phi - I` error `5.14e-13`
  - max lift error `1.96e-13`
  - max projector error `1.06e-11`
  - max `V` formula error `7.39e-14`
  - max `Q` formula error `5.59e-14`
- exact closure criterion remains sharp:
  - adapted `Q` norm `1.02e-29`
  - adapted off-plane norm `1.91e-15`
  - random `eta` on same family `0.715`
- commuting solver stress remains exact:
  - exact-match rate against exhaustive search `1.0`
  - adapted `eta` max `1.46e-15`
  - random median `eta` across trials `0.718`
- bridge stress remains exact:
  - adapted `Q` max `1.34e-28`
  - adapted `eta` max `1.94e-15`
  - adapted scaled quartic residual median `1.14e-12`

## Supplemental Checks Added Here

These checks were run after the original paper-side sweeps to harden the newer conclusions that matter for module policy.

### 1. `Tr(Pi M) = S + L` identity check

For symmetric families with `M = sum_a D_a^2`, a 200-trial random check gave:
- max absolute residual in `Tr(Pi M) - (S + L)`: `2.84e-14`

Interpretation:
- the `sum_a D_a^2` rule is exact for total captured curvature
- this confirms its correct role as a warm start or frontier point, not as an approximate-closure theorem by itself

### 2. Reducible noncommuting witness

A controlled `5 x 5`, rank-2 example was built with:
- one exact low-energy invariant `2`-plane
- a high-energy noncommuting `3`-dimensional complement

Observed scores:
- exact invariant plane:
  - `S = 2.05`
  - `eta = 0`
- `M = sum_a D_a^2` eigenspace:
  - `S = 170.83`
  - `eta = 0.1929`
- best sampled single-matrix eigenspace:
  - `S = 154.74`
  - `eta = 0.2146`
- best sampled observer under `eta <= 0.08`:
  - `S = 134.07`
  - `eta = 0.0580`

Interpretation:
- pure exact closure can select an operationally weak observer
- `sum_a D_a^2` and single-matrix eigenspaces find energetic planes, but with much too much leakage
- the correct noncommuting object is the leakage-visibility frontier, not one scalar objective

### 3. Tower flag obstruction witness

Using
`D1 = diag(2, 1, -1)` and
`D2 = [[3,0,0],[0,0,1],[0,1,0]]`,
the rank-1 observer on `span(e1)` is exactly adapted:
- `S = 13.0`
- `eta = 0`

But the quotient family on `span(e2, e3)` has compressed commutator norm:
- `||[D1_hat, D2_hat]||_F = 2.8284`

A dense sweep over all rank-2 planes containing `e1` gave:
- minimum observed `eta = 0.125`

Interpretation:
- exact coarse closure does not automatically refine
- the correct exact refinement condition is the existence of a common invariant flag, equivalently a common invariant subspace in the quotient family

### 4. Stationarity and order-gain spot checks

The revised TeX makes two additional theorem-level claims especially visible:
- the Grassmannian stationarity equation for leakage
- the even-order gain on adapted bridge branches

Numerical spot checks now recorded:
- stationarity residual `|| [Pi, sum_a [D_a, [D_a, Pi]]] ||_F` on 40 theorem-built exact adapted commuting observers:
  - median `7.64e-16`
  - max `3.10e-15`
- bridge remainder order after subtracting the `eps^2 V` term:
  - seeded random observers show quartic scaling with median log-log slope `3.98`
  - adapted observers drive the same remainder to machine floor across the tested `eps` range, so the expected higher-order cancellation is numerically supported but not slope-identifiable beyond floating-point precision in the current harness

Interpretation:
- the stationarity equation is consistent with the exact adapted endpoint
- the adapted bridge order gain is numerically present, but the current bridge harness is already so exact that the upgraded slope is hidden by machine precision rather than by a visible finite-sample exponent

## Integration Boundary

The safe module promotion boundary is now clear.

Promote now:
- exact whitening helper for perturbations
- exact observer-from-subspace constructor
- exact closure scores `L`, `S`, `eta`
- exact commuting-family observer synthesis
- exact bridge-adapted example and tests

Do not promote yet:
- public noncommuting Grassmannian optimiser
- public constrained-`eta` selector policy
- `sum_a D_a^2` as a public final observer
- tower-refinement API beyond exact flag language

## Required Cleanup Before Integration

### 1. Script execution consistency

The paper-side scripts are not execution-consistent yet.

Observed:
- `python papers/closure_adapted_audit.py` works
- `python papers/closure_adapted_use_cases.py` fails as a script because it imports `from papers...`
- `python papers/closure_adapted_implications.py` fails for the same reason
- both work as modules with `python -m papers...`

Required fix:
- choose one execution mode and make all closure-adapted study files consistent with it

### 2. Approximate-regime scripts and reports still lag the revised TeX

The paper is now aligned, but the approximate study scripts and reports still overstate the design status of the aggregate-mode heuristic.

Files to update after the revised TeX lands:
- [closure_adapted_implications.py](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_implications.py)
- [closure_adapted_use_cases.py](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_use_cases.py)
- [closure_adapted_implications_report.md](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/closure_adapted_outputs/closure_adapted_implications_report.md)

Specifically, module-facing language must now reflect:
- noncommuting design is a frontier in `(L, S)` or `(eta, S)`
- pure `eta` minimisation is a closure endpoint and diagnostic, not the generic task-optimal target
- `sum_a D_a^2` is retained as a warm start only
- tower refinement requires the invariant-flag criterion

### 3. Documentation alignment

Now that the revised adapter paper is in place, the following must be brought into exact agreement:
- [Closure_Adapted_Observers.tex](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/papers/Closure_Adapted_Observers.tex)
- [theorem_map.md](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/docs/theorem_map.md)
- [api_notes.md](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/docs/api_notes.md)
- [validation_note.md](C:/Users/jrdun/family/observer_geometry_workspace_v0.2_freeze_2026-04-08/docs/validation_note.md)

## Required Release-Gate Tests

Before any `src/` promotion, add release-gate tests for:

1. exact whitening and observer parameterisation
- `whitened_delta(H, Delta)` symmetry and covariance checks
- `observer_from_subspace(H, B)` round-trips against the existing lift/projector identities

2. closure scores
- `L >= 0`, `S >= 0`, `eta = L / (L + S)` on the theorem domain
- `L = 0` exactly on constructed invariant-subspace families

3. commuting exact solver
- exhaustive agreement in small dimensions against all common eigenspace subsets
- exact `eta ~= 0` for theorem-built observers

4. bridge-adapted exact demo
- adapted observer suppresses `Q`
- adapted observer suppresses quartic residual after removing the `eps^2 V` term
- adapted observer beats seeded random observers on leakage and visible score

5. `Tr(Pi M) = S + L` identity
- explicit theorem-facing test on random symmetric families
- this is the clean test that justifies `sum_a D_a^2` as a warm start only

6. tower obstruction
- explicit `R^3` counterexample where rank-1 exact closure exists and rank-2 refinement fails
- quotient-family criterion checked directly

7. hostile misuse rejection
- reject malformed `B`
- reject non-symmetric family members
- make commuting-only solver fail loudly on noncommuting input rather than silently returning a heuristic answer

## Recommended Integration Order

1. align paper-side reports and study scripts with the new frontier and flag conclusions
2. implement only the exact closure-adapted kernel in `src/nomogeo`
3. add release-gate tests
4. update theorem map, API notes, and validation note
5. leave noncommuting optimisation public surfaces deferred until a frontier API is actually implemented

## Judgement

The integration decision is no longer blocked by the exact mathematics.

The remaining blockers are packaging and boundary discipline:
- keep the exact commuting/invariant-subspace center
- do not over-promote the noncommuting heuristic layer
- force the docs, tests, and public language to use the frontier and invariant-flag conclusions exactly
