# Closure-Adapted Observers Audit

Date: 2026-04-09

Source paper:
- `papers/Closure_Adapted_Observers.tex`

Audit harness:
- `papers/closure_adapted_audit.py`
- `papers/closure_adapted_outputs/closure_adapted_audit_results.json`

## What The Paper Claims

The paper turns observer design into an exact theorem-native problem.

Given:
- a latent precision `H`
- a visible dimension `m`
- a task family of perturbations `{Delta_a}`

it says:
- after `H`-whitening, an observer is just an `m`-plane `Ran(B)`
- the first visible response is the compressed whitened perturbation
- the quartic defect is exactly the off-plane action
- exact closure is equivalent to invariance of the visible plane under the whitened perturbation family
- in the commuting case, exact-closure observers are obtained by common eigenspace selection
- among exact-closure observers, the best one is the span of the `m` common modes of largest aggregate spectral energy

## Verification Outcome

The paper checked out on the exact layers I tested.

### 1. Whitened normal form

Across 50 random trials:
- max `Phi - I` error: `5.14e-13`
- max lift error: `1.96e-13`
- max projector error: `1.06e-11`
- max `V` formula error: `7.39e-14`
- max `Q` formula error: `5.59e-14`

So the normal form is not just suggestive. It matches the current kernel numerically.

### 2. Exact closure criterion

In a constructed invariant-subspace example:
- adapted observer `Q` norm: `1.02e-29`
- adapted off-plane action norm: `1.91e-15`

For a random observer on the same task:
- `Q` norm: `4.48`
- leakage ratio `eta`: `0.715`

So the paper's closure criterion is behaving exactly as advertised.

### 3. Commuting-family reduction

For an 8D, 3-task commuting family:
- the paper's top-`mu_i` rule selected indices `[0,1,2]`
- exhaustive search over all `8 choose 3 = 56` common-eigenspace observers gave the same subset as the best exact-closure observer
- adapted observer score:
  - `S = 29.14`
  - `H ~= 0`
  - `eta ~= 0`
- random observer comparison:
  - median `eta = 0.662`
  - minimum random `eta = 0.406`
  - median random `S = 7.22`

This is the strongest implication of the paper. In the commuting case, observer synthesis is effectively solved.

## Strong Demo: Bridge-Adaptive Observer

I used the repo's `dv_bridge(H0, Jhat).delta_dv` as the task perturbation and built the closure-adapted observer from the top eigenspace of the whitened bridge correction.

Results:
- adapted `Q` norm: `6.48e-31`
- adapted `eta`: `~ 0`
- adapted visible score norm `||V||`: `1.53`

Random 2D observers on the same bridge task:
- median `eta = 0.645`
- minimum random `eta = 0.032`
- median random `Q` norm: `0.272`
- maximum random `||V||`: `1.14`

Quartic-onset contrast at `eps = 0.2`:
- adapted observer residual after subtracting the `eps^2 V` term, scaled by `eps^4`: `1.42e-12`
- random observer median scaled quartic residual: `0.265`

Interpretation:
- the adapted observer nearly eliminates quartic hidden birth completely
- it also retains more visible signal than any random observer in the sample
- this is exactly the kind of capability the module does not currently expose, but now could

## Implications For The Python Module

This paper suggests a real new capability, not just a new note.

### What the module already has

- `visible_precision`
- `visible_geometry`
- `canonical_lift`
- `hidden_projector`
- `local_visible_calculus`
- `dv_bridge`

### What it does not have yet

- `H`-whitened observer parameterisation
- closure-adapted observer construction from a perturbation family
- commuting-family solver by common spectral mode selection
- leakage objective `eta_F(B)` as a first-class observer-design quantity

### Minimal useful API implied by the paper

Something like:
- `observer_from_subspace(H, B)`
- `whitened_perturbation(H, Delta)`
- `closure_scores(H, family, B)` returning `S`, `H`, `eta`
- `closure_adapted_observer(H, family, rank, mode=\"commuting\")`

### Best first release target

Implement the commuting-family case first.

Reason:
- it is exact
- it is now verified
- it gives a clean, impressive capability jump
- it has an immediate bridge demo already supported by the repo

## Bottom Line

`Closure_Adapted_Observers` looks correct on its exact claims and has direct module implications.

The most important consequence is:
- the package can plausibly grow from "given an observer, compute visible geometry" to
- "given a task family, synthesize an observer that exactly closes it when possible"

That is a meaningful capability increase.
