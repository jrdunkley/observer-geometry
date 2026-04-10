# API Notes

## Public Entry Points

- `visible_precision(H, C)` returns the visible precision matrix.
- `canonical_lift(H, C)` returns `L_{C,H}`.
- `hidden_projector(H, C)` returns `P_{C,H}`.
- `visible_geometry(H, C)` returns all three together as a `VisibleGeometryResult`.
- `local_visible_calculus(H, C, Delta)` returns `phi`, `V`, `Q`, the lift/projector, and the determinant-curvature split.
- `whitened_perturbation(H, Delta)` returns the exact whitened perturbation `H^(-1/2) Delta H^(-1/2)`.
- `observer_from_subspace(H, B)` returns the observer `C = B^T H^(1/2)` attached to an orthonormal visible plane.
- `closure_scores(H, family, B)` returns exact leakage `L`, retained visible score `S`, total captured curvature `L + S`, and hidden fraction `eta`.
- `leakage_channels(H, Delta, B)` returns the exact leakage Gram operator and its singular leakage channels for one perturbation relative to one visible plane.
- `compare_observers(H, family, B_left, B_right)` returns an exact same-rank comparison witness, including dominance flags when one observer has no more leakage and no less visible score than the other.
- `closure_adapted_observer(H, family, rank, mode="commuting_exact")` returns the exact commuting-family observer selected by aggregate common-mode energy.
- `fixed_observer_coordinates(H, C)` returns the exact fixed-observer chart `(Phi, R, K)` together with the observer-fixed adapted basis used to define it.
- `reconstruct_precision_from_fixed_observer_coordinates(phi, hidden_block, coupling, adapted_basis)` reconstructs `H` exactly from the public fixed-observer chart.
- `observer_transition(H, C_left, C_right)` returns the exact block transition law between two fixed-observer charts of the same visible rank.
- `connection_current(H, C, Delta)` returns the exact fixed-observer coupling velocity, current `J`, and forcing `Q` attached to one tangent direction.
- `forcing_from_current(phi, hidden_block, current)` evaluates the exact current-to-forcing identity `Q = Phi J R^(-1) J^T Phi`.
- `pi_from_hidden_load(Lambda)` and `hidden_load_from_pi(Pi)` convert between reduced hidden load and the precision-side contraction variable.
- `pi_rhs(Pi, A_cpl)`, `lambda_rhs(Lambda, A_cpl)`, `clock_rate(A_cpl)`, and `support_stratum_transport(Lambda, A_cpl)` evaluate the exact support-stable stratum transport diagnostics.
- `restart_hidden_load_birth(lambda_before, old_basis, new_basis)` and `restart_hidden_load_death(lambda_before, old_basis, survivor_basis)` apply the forced finite positive restart maps in explicit support coordinates.
- `kernel_schur_jet_from_coefficients(coefficients)` classifies finite-order support events from the kernel Schur-complement jet.
- `semisimple_event_block(jet, ...)` returns the universal semisimple pole and clock diagnostics for an active event block.
- `local_coupled_birth(H, H_dot, H_ddot, C, C_dot, Z=None)` extracts the hidden-basis-invariant coupled local birth tensor from explicit derivative data.
- `sampled_interval_leakage`, `sampled_interval_stationarity`, `sampled_interval_closure_check`, and `interval_hessian_at_exact_family` expose sampled interval-family leakage diagnostics.
- `dv_bridge(H0, Jhat)` returns `H_DV`, `Delta_DV`, and the Gram factor.
- `hidden_load(T, X)` returns the support-aware hidden load and metadata.
- `visible_from_hidden_load(T, Lambda, ..., lambda_representation=...)` reconstructs `X`.
- `inverse_visible_class(T, Lambda, ..., lambda_representation=...)` is the theorem-facing alias for the same fixed-ceiling inverse map.
- `hidden_contraction(Lambda)` returns the canonical contraction factor `K = (I + Lambda)^(-1/2)`.
- `load_from_hidden_contraction(K)` recovers the hidden load from `K^T K`.
- `canonical_hidden_realisation(T, X)` and `minimal_hidden_realisation(T, X)` work on active-support coordinates and return explicit block matrices.
- `transport_hidden_load(Lambda, M)` returns the theorem-native two-step transported load.
- `clock(Lambda)` returns `log det(I + Lambda)`.
- `observed_covariance`, `gaussian_forward_kl`, `gaussian_reverse_kl`, `gaussian_hellinger_squared`, and `gaussian_data_processing_contraction` provide the thin quotient adapter layer.
- `batch_map(func, tasks, ...)` is the thin deterministic batch wrapper.

## Support-Aware Convention

The hidden-load layer follows the papers: the real operator lives on `S = Ran(T)`.

- `hidden_load` returns:
  - `lambda_`: ambient embedding with zeros off support
  - `reduced_lambda`: the active-support operator
  - `support_basis`: the basis used for the active support
- `visible_from_hidden_load` accepts either:
- `inverse_visible_class` accepts the same coordinate choices and exists to make the fixed-ceiling inverse theorem explicit at the API surface
  - an ambient `n x n` load with `lambda_representation="ambient"`, or
  - a reduced `rank(T) x rank(T)` load with `lambda_representation="reduced"`
- when `rank(T) = n`, shape alone is ambiguous and the function requires the representation to be stated explicitly
- this inverse surface is exact only after fixing a ceiling/reference `T`; it does not invert the global forward map `(H, C) -> Phi_C(H)`
- `clock` and `transport_hidden_load` are theorem-domain functions and reject indefinite hidden loads
- `clock(Lambda)` is an exact hidden-burden scalar, not a universal observer-quality or panel-quality score
- the associative downstairs object is the contraction factor `K = (I + Lambda)^(-1/2)`, not the raw visible-effect or raw load coordinate
- long chains should therefore be composed in contraction-factor coordinates and converted back with `load_from_hidden_contraction`

For applied panel or observer ranking, use task-family-aware closure or perturbation-family analysis in addition to hidden-load summaries.

`canonical_hidden_realisation` and `minimal_hidden_realisation` return support-coordinate block matrices, because that is the canonical domain of the theorem.

## Result Objects

The package uses small frozen dataclasses instead of mutable object graphs:

- [`src/nomogeo/types.py`](../src/nomogeo/types.py)

Each result carries `LinearAlgebraMetadata` with tolerance, support rank, and method information.

Closure-adapted results added:

- [`ClosureScoresResult`](../src/nomogeo/types.py) carries exact leakage, visible score, hidden fraction, total captured curvature, and the visible projector.
- [`LeakageChannelsResult`](../src/nomogeo/types.py) carries the exact coupling map, leakage Gram operator, and the visible/hidden singular leakage channels.
- [`ClosureAdaptedObserverResult`](../src/nomogeo/types.py) carries the selected basis `B`, observer `C`, projector, exact scores, the computed common basis, and per-mode spectral energies.
- [`ObserverComparisonResult`](../src/nomogeo/types.py) carries two exact score objects plus same-rank dominance / inefficiency verdicts.
- [`FixedObserverCoordinatesResult`](../src/nomogeo/types.py) carries the exact fixed-observer chart `(phi, hidden_block, coupling)`, the adapted basis used to define it, and conditioning metadata.
- [`ObserverTransitionResult`](../src/nomogeo/types.py) carries the left/right charts together with the exact block transition data and reconstruction residual.
- [`ConnectionCurrentResult`](../src/nomogeo/types.py) carries the chart data, coupling velocity, current `J`, forcing, and the directly computed local `Q`.
- [`SupportStratumTransportResult`](../src/nomogeo/types.py) carries reduced `Lambda`, `Pi`, `A_cpl`, the two transport right-hand sides, clock rate, eigenvalue bounds, and PSD-domain status.
- [`SupportRestartResult`](../src/nomogeo/types.py) carries birth/death restart data, old/new support bases, and the coordinate map between them.
- [`KernelJetResult`](../src/nomogeo/types.py) carries the kernel/gap split, effective Schur-complement coefficients, first nonzero order, leading eigenvalues, and event classification.
- [`SemisimpleEventBlockResult`](../src/nomogeo/types.py) carries finite-order semisimple pole, clock, and desingularisation diagnostics.
- [`LocalCoupledBirthResult`](../src/nomogeo/types.py) carries the local extractor outputs `Phi`, `L`, `R`, `V`, `B`, `beta`, `Q`, observer tensor, `W`, support basis, and `A_cpl`.
- [`SampledIntervalLeakageResult`](../src/nomogeo/types.py) carries sampled leakage, visible score, stationarity residual, projector, weights, and sampled exact-closure status.
- [`IntervalHessianResult`](../src/nomogeo/types.py) carries sampled interval Hessian and spectral-gap rigidity diagnostics.

## Closure-Adapted Boundary

The current closure-adapted public surface is intentionally exact and narrow.

- `closure_adapted_observer(..., mode="commuting_exact")` is implemented only for pairwise commuting symmetric whitened families.
- noncommuting frontier optimisation is not yet a public API.
- the `sum_a D_a^2` rule is not exposed as a final observer rule; its paper-level status is warm start only.
- tower refinement language is documentation-only at present; no public refinement API is exposed.
- adapted-layer metadata now includes an estimated condition number for `H`; high-condition notes are diagnostic warnings, not theorem defects.

## Connection Boundary

The current connection-facing public surface is also intentionally narrow.

- only the fixed-observer chart from `Connection_Flatness` is public
- the chart is defined using an observer-fixed adapted basis, not an `H`-dependent canonical lift basis
- `phi` is the intrinsic visible object; `hidden_block` and `coupling` are exact chart coordinates relative to the chosen observer/basis
- the zero-hidden full-visible edge case is supported explicitly; in that regime `hidden_block` is `0 x 0`, `coupling` is empty, and the exact forcing term vanishes
- full varying-observer connection geometry and global non-fixed-observer transport are not public API; the only moving-observer public surface is the explicit-derivative local extractor `local_coupled_birth`
- conditioning metadata on `H` is diagnostic only; very ill-conditioned inputs can degrade chart numerics before any theorem fails

## Field Boundary

The `0.30.0` field surface is exact but narrow.

- support-stratum transport diagnostics require a reduced support coordinate system
- `A_cpl` may be indefinite; when it is not PSD, the generator is reported but positive hidden-load transport is not licensed
- restart maps require explicit nested orthonormal support bases and return coordinate maps
- kernel-jet classification uses Taylor coefficients, where `coefficients[n]` multiplies `epsilon^n`
- kernel jets govern leading small-eigenvalue behaviour and near-zero inertia, not exact equality with the full Schur-complement spectrum
- sampled interval diagnostics certify only the supplied samples; they are not continuum certificates without extra analytic assumptions
- no noncommuting optimiser, sampled branch-continuation solver, global stratified integrator, or Hermitian packet API is public in this release

