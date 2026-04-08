# Validation Note

## Implemented Now

- exact visible precision, variational character, tower law, and Schur form
- canonical lift and hidden projector
- local visible calculus with the positive Gram form of `Q`
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
- deterministic batch wrapper

## Deliberately Deferred

- full finite-observer optimisation layer from `Finite_Observation`
- full quotient-observation framework beyond the thin contraction adapter
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
- no theorem domains were widened
- no application-layer features were added

## Verification Run

Current smoke result:

```text
47 passed
```

Install surface verified:

- `python -m pip install -e .`
- `pytest`
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

