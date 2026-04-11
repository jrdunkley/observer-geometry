# Validation Note

## Implemented Now

- exact visible precision, variational character, tower law, and Schur form
- canonical lift and hidden projector
- local visible calculus with the positive Gram form of `Q`
- exact closure-adapted whitening, observer-from-subspace construction, closure
  scores `(L, S, eta)`, leakage-channel diagnostics, same-rank observer
  comparison, and commuting-family observer synthesis
- exact fixed-observer chart coordinates `(Phi, R, K)`, exact chart
  reconstruction, observer-transition law, and fixed-observer current / forcing
  identities from `Connection_Flatness`
- exact intrinsic local quadratic ensemble evaluator for arbitrary surjective
  observers, returning only samplewise visible precisions and functorial
  summaries
- exact ceiling-mediated local quadratic ensemble evaluator requiring declared
  ceilings before hidden loads and clocks are exposed
- coordinate-split local quadratic ensemble convenience wrapper with descriptive
  clock summaries over samplewise outputs
- support-aware hidden-load parametrisation beneath a ceiling
- visible reconstruction from hidden load with explicit representation rules in
  ambiguous full-rank cases
- first-class fixed-ceiling inverse theorem surface through
  `inverse_visible_class(T, Lambda, ...)`
- canonical and minimal hidden realisations
- hidden-load transport and additive determinant clock on the theorem domain
  `Lambda >= 0`
- hidden contraction factors and recovery of hidden loads from contraction
  factors, exposing the associative downstairs composition object
- thin Donsker-Varadhan bridge
- thin quotient-side Gaussian contraction checks and observer-collapse checks
- exact rank-one and rank-k covariance/Fisher perturbation diagnostics for coordinate visible splits, including the white/aligned one-channel sector, coloured two-rank boundary, and rank-at-most-`2k` update bound
- residual-margin ordering certificate for bounded branch/observer residuals
- simple-spectrum closure certificate for exact common-closure existence or obstruction among anchor-coordinate eigenspans
- exact variable-precision affine-hidden reducer with fibre-volume term
- exact staged affine-hidden Schur elimination helper
- exact finite affine-hidden branch-reversal diagnostic and guarded fibre-dominance tuple
- finite weighted-family leakage/visibility frontier evaluator
- exact-branch Hessian diagnostic for already-exact weighted-family branches
- exact general graph-frontier Hessian at declared observers, declared-ladder dimension-cost intervals, and sufficient declared local frontier certificates
- deterministic batch wrapper

## Deliberately Deferred

- full finite-observer optimisation layer from `Finite_Observation`
- noncommuting closure-adapted frontier optimisation and selector policy
- generic closure obstruction certificates without a simple-spectrum anchor
- full quotient-observation framework beyond the thin contraction adapter
- generic non-Gaussian branch or observer selector
- residual estimation for full-law corrections
- noncommuting closure irreducibility/frontier API
- noncommuting weighted-family optimiser or selector policy
- arbitrary fibre-cumulant full-law branch engine
- probability-support event engine
- Hermitian packet API
- hidden-load or clock ensemble summaries for general observers without an explicit ceiling contract
- downstream consciousness analysis and empirical pipelines
- plotting, notebooks, data ingestion, benchmarks, and application-specific
  wrappers

## Tolerances

Central defaults:

- `atol = 1e-10`
- `rtol = 1e-8`

Rank and PSD decisions use the spectrum-aware cutoff described in
[docs/numerical_conventions.md](numerical_conventions.md).

## Hidden-Load Domain Rules

- `hidden_load(T, X)` works on the theorem domain beneath a ceiling and returns:
  - `reduced_lambda`: active-support coordinates
  - `lambda_`: ambient embedding
- the exact inverse theorem is ceiling-conditioned: once `T` is fixed, PSD
  hidden loads on `Ran(T)` are in bijection with support-preserving visible
  objects beneath `T`
- `inverse_visible_class(T, Lambda, ...)` is the theorem-facing alias for the
  same inverse map implemented by `visible_from_hidden_load`
- `visible_from_hidden_load(T, Lambda, ...)` requires an explicit coordinate
  declaration whenever full-rank ceilings make reduced and ambient shapes
  coincide
- `clock(Lambda)` and `transport_hidden_load(Lambda, M)` enforce the theorem
  domain and reject indefinite loads
- zero-support ceilings are handled explicitly: the inverse map returns the zero
  matrix and the zero-size clock is `0`
- `transport_hidden_load` is retained as the theorem's exact two-step
  Gram-shadow transport law
- long associative composition is exposed through `hidden_contraction(Lambda)`
  and `load_from_hidden_contraction(K)`, matching the paper's statement that
  the true associative object lives downstairs rather than in the visible
  product alone

## Release-Gate Tests Added

Track A:

- full lift/projector identity suite including `P^T H = H P`, variational
  energy splitting, and `Ran(L)` being `H`-orthogonal to `ker(C)`
- pure-visible and hidden-annihilating local-calculus families
- determinant-curvature theorem check and `O(t^3)` local remainder control
- three-stage and four-stage tower/composition checks
- latent-basis and visible-basis covariance checks
- hidden-load representation consistency, support-gauge covariance,
  Loewner monotonicity, rank/clock identities, and near-ceiling recovery
- direct fixed-ceiling inverse-theorem round-trips, order reversal, zero-support
  behaviour, full-rank ambiguity rejection, and large-load conditioning checks
- two-step hidden-load transport checks plus downstairs associative-contraction
  checks over long chains
- strengthened local hidden-birth checks across step ladders
- DV bridge trivial-current, sign symmetry, basis covariance, and
  small-parameter scaling checks
- seeded quotient/divergence contraction checks with both strict and equality
  cases

Track B:

- rejection of malformed visible, hidden, and batch inputs with explicit
  exception types
- near-boundary tolerance tests
- input immutability checks
- batch determinism, failure propagation, and forced process-fallback checks
- seeded validation-sweep harness smoke test
- malformed closure-adapted basis and family rejection
- commuting-only rejection for noncommuting closure-adapted input
- malformed fixed-observer chart, transition, and current inputs with explicit
  exception types

Closure-adapted release-gate checks:

- exact whitening and observer normal-form round-trips
- exact closure on invariant-subspace families
- exact leakage-channel / Gram-operator recovery checks
- exact same-rank observer comparison / dominance checks
- exhaustive commuting-family agreement against all common-eigenspace subsets in
  small dimensions
- bridge-adapted exact observer checks for vanishing `Q`, vanishing `eta`, and
  quartic-residual suppression
- explicit `Tr(Pi M) = L + S` identity checks
- explicit `R^3` tower-obstruction witness checks
- empirical hardening checks against the local micro-real Iris, leaderboard,
  and Bell bundles, verifying that real noncommuting families are rejected
  loudly while data-shaped exact commuting families still solve exactly

Connection release-gate checks:

- exact fixed-observer chart round-trips and `Phi_C(H)` recovery
- exact observer-transition-law recovery against direct recomputation
- exact current-to-forcing identity checks `Q = Phi J R^(-1) J^T Phi`
- fixed-`K` zero-forcing edge checks
- explicit full-visible zero-hidden-sector checks
- empirical hardening checks against the local micro-real Iris, leaderboard,
  and Bell bundles
- hostile conditioning checks up through rotated families with
  `cond(H) ~ 1e8`, with explicit warning thresholds rather than silent failure

Perturbation and residual-margin checks:

- exact white-sector hidden-gap formula and one-channel rank
- coloured-background rank-two witness and Sherman-Morrison formula residual
- block-aligned coloured sector remaining one-channel
- rank-k covariance/Fisher Woodbury formula and rank-at-most-`2k` bound
- malformed covariance/signal/visible-dimension/epsilon rejection
- strict margin certificate, equality boundary, and adversarial reversal case
- simple-spectrum closure certificate with exact block-invariant witness
- simple-spectrum generic obstruction certificate with all coordinate cross-blocks nonzero
- degenerate-anchor rejection

Coordinate ensemble checks:

- samplewise agreement with `visible_precision` and `hidden_load`
- intrinsic `Phi` covariance under latent basis changes
- ceiling-mediated clock and hidden-load spectrum invariance under visible basis changes when ceilings transform by congruence
- exact agreement between the coordinate convenience wrapper and explicit ceiling-mediated mode
- descriptive clock summary agreement with direct NumPy summaries
- witness that the clock of the mean Hessian need not equal the mean sample clock
- malformed shape, empty sample, visible-dimension, and non-SPD rejection
- malformed or too-small ceilings rejected

Affine-hidden and weighted-frontier checks:

- numerical hidden-fibre integration against the variable-precision exact
  action formula
- explicit log-determinant branch-flip witness where variational action is flat
- exact finite affine-hidden branch-reversal witness
- guarded fibre-dominance denominator-floor witness
- staged single-coordinate Schur eliminations against one-step elimination
- malformed affine-hidden shape and non-SPD hidden-precision rejection
- weighted-family energy split against random finite families and agreement
  with unweighted `closure_scores`
- exact-branch Hessian sign/status convention against graph-chart finite
  differences
- declared-ladder dimension-cost interval phase diagram
- general graph-frontier Hessian exact-branch reduction and non-exact stationary
  breakpoint at `sqrt(3)`
- declared local frontier certificate max/min/vacuity checks
- noncommuting internal block exact-branch Hessian witness
- rejection of non-exact branches and invalid full-rank Hessian requests

## Fixes Triggered By This Campaign

- added `hidden_contraction` and `load_from_hidden_contraction` in
  [src/nomogeo/hidden.py](../src/nomogeo/hidden.py)
  - the release-gate campaign exposed that raw load transport is the exact
    two-step Gram-shadow law, but not the associative object for long chains
  - the paper already states this boundary; the fix was to expose the actual
    downstairs contraction object rather than to loosen `transport_hidden_load`
- added `inverse_visible_class` in
  [src/nomogeo/hidden.py](../src/nomogeo/hidden.py)
  - this does not widen the kernel; it promotes the already-valid fixed-ceiling
    inverse theorem into an explicit public surface
- added exact closure-adapted kernel surfaces in
  [src/nomogeo/adapted.py](../src/nomogeo/adapted.py)
  - this promotion is intentionally limited to the exact commuting/invariant-subspace center
  - noncommuting frontier optimisation remains deferred
  - the layer now also exposes exact leakage-channel diagnostics and same-rank
    observer inefficiency witnesses
- added empirical and conditioning hardening harness in
  [tools/adapted_hardening.py](../tools/adapted_hardening.py)
  - this pushes the exact closure-adapted slice against local real bundles and hostile conditioning rather than only synthetic theorem fixtures
- added exact fixed-observer chart surfaces in
  [src/nomogeo/connection.py](../src/nomogeo/connection.py)
  - this promotion is intentionally limited to the fixed-observer exact chart
    and current identities from `Connection_Flatness`
  - varying-observer geometry is still deferred
- added connection hardening harness in
  [tools/connection_flatness_hardening.py](../tools/connection_flatness_hardening.py)
  - this pushes the fixed-observer chart against synthetic gauge checks,
    hostile conditioning, and local micro-real bundles
- made `BatchTaskError` process-pool round-trip safe in
  [src/nomogeo/exceptions.py](../src/nomogeo/exceptions.py)
  - the full-suite rerun exposed that remote batch failures were breaking process-backend failure propagation
  - the fix preserves the existing exception surface while making it pickle-safe
- added install-surface smoke in
  [tools/install_surface_smoke.py](../tools/install_surface_smoke.py)
  - this verifies that a user-style import resolves to this workspace and that
    the public fixed-observer connection API survives direct import/use
- no theorem domains were widened
- no application-layer features were added

## Verification Run

Current smoke result:

```text
100 passed
```

Current collection check:

```text
100 tests collected
```

0.3.1 planning stress result:

```text
overall_passed: true
```

This planning stress result is not a module release gate. It checks the
no-canonical-ceiling obstruction, intrinsic quotient covariance,
variable-precision affine-hidden elimination, staged entropic tower
composition, weighted-family energy/stationarity identities, branch-Hessian
status semantics, typed residual margins, support-event charge scaling, the
positive packet lift, and finite weighted-family identities on the local Iris,
leaderboard, and Bell micro-real bundles.

Install surface verified:

- `python -m pip install -e .`
- `pytest`
- `python tools/install_surface_smoke.py`
- direct execution of all files in
  [examples](../examples) against the
  installed package, without `sys.path` edits
- `python -m pip wheel . -w dist --no-deps`
- clean wheel install into `.wheeltest` with `--system-site-packages`
- import and minimal example execution from the wheel-installed environment

Residual sweep verified:

- `python tools/validation_sweep.py`

Worst observed sweep values with seed `20260408`, dimensions `2..8`, and `4`
draws per dimension:

- max lift/projector residual: `2.674941282677552e-14`
- max `H`-projector symmetry residual: `5.886197714175869e-14`
- max energy-split residual: `2.8421709430404007e-13`
- max tower-law residual: `3.284872374109682e-14`
- max curvature-split residual: `3.180913144220955e-04`
- max normalized local `O(t^3)` remainder `||R_t|| / t^3`:
  `6.940752273512302e-06`
- max hidden rank mismatch: `0`
- max hidden clock residual: `3.552713678800501e-15`
- max hidden Loewner-order violation: `1.0195833073384965e-16`
- max inverse-theorem `Lambda -> X -> Lambda` relative residual:
  `2.587909939792647e-15`
- max inverse-theorem `X -> Lambda -> X` relative residual:
  `2.2201725879984535e-15`
- max inverse-theorem order-reversal violation: `1.0195833073384965e-16`
- max two-step transport residual against downstairs factor composition:
  `1.8991764023475337e-14`
- max downstairs associativity residual: `4.292205267163458e-15`
- max long-chain transport residual: `0.0`
- max long-chain clock residual: `1.2434497875801753e-14`
- max local hidden-birth first-order residual: `1.778784990776859e-02`
- max DV visible residual: `1.1797712401204017e-07`
- max DV hidden normalized residual: `1.1723687933257149e-01`
- max divergence contraction violation: `0.0`

Covered checks:

- theorem identities for visible precision, lift, projector, variational
  splitting, and quotient towers
- coordinate covariance in latent and visible bases
- exact structural local-calculus families and determinant-curvature checks
- hidden-load representation, support, rank, clock, near-ceiling, and
  monotonicity checks
- explicit fixed-ceiling inverse-theorem round-trips, order reversal, ambiguity
  rejection, and conditioning checks
- two-step hidden-load transport plus downstairs associative-composition checks
- local hidden-birth asymptotics
- DV bridge symmetry, covariance, and small-parameter checks
- seeded quotient-divergence contraction checks
- hostile misuse rejection
- near-boundary stability
- input immutability
- deterministic batch behaviour and failure propagation
- seeded property-sweep dashboard generation

