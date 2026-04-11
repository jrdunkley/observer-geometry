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
- `simple_spectrum_closure_certificate(family, rank, anchor_index=0)` checks whether a simple-spectrum family member leaves any rank-`rank` coordinate eigenspan exactly invariant for the rest of the family. It can certify absence of exact common closure in that simple-spectrum setting, but it is not a noncommuting optimiser.
- `fixed_observer_coordinates(H, C)` returns the exact fixed-observer chart `(Phi, R, K)` together with the observer-fixed adapted basis used to define it.
- `reconstruct_precision_from_fixed_observer_coordinates(phi, hidden_block, coupling, adapted_basis)` reconstructs `H` exactly from the public fixed-observer chart.
- `observer_transition(H, C_left, C_right)` returns the exact block transition law between two fixed-observer charts of the same visible rank.
- `connection_current(H, C, Delta)` returns the exact fixed-observer coupling velocity, current `J`, and forcing `Q` attached to one tangent direction.
- `forcing_from_current(phi, hidden_block, current)` evaluates the exact current-to-forcing identity `Q = Phi J R^(-1) J^T Phi`.
- `intrinsic_local_quadratic_ensemble(hessians, C)` pushes a stack of local SPD Hessian/Fisher samples through the intrinsic quotient map and returns samplewise `Phi_C(H_i)` plus log-determinant summaries. It deliberately does not return hidden load or clock.
- `ceiling_mediated_local_quadratic_ensemble(hessians, C, ceilings)` adds an explicit ceiling family and returns samplewise hidden loads, ranks, clocks, and descriptive clock summaries.
- `coordinate_local_quadratic_ensemble(hessians, visible_dim)` is the coordinate-split convenience wrapper where `C=[I 0]` and `T_i=(H_i)[:visible_dim, :visible_dim]`.
- `pi_from_hidden_load(Lambda)` and `hidden_load_from_pi(Pi)` convert between reduced hidden load and the precision-side contraction variable.
- `pi_rhs(Pi, A_cpl)`, `lambda_rhs(Lambda, A_cpl)`, `clock_rate(A_cpl)`, and `support_stratum_transport(Lambda, A_cpl)` evaluate the exact support-stable stratum transport diagnostics.
- `restart_hidden_load_birth(lambda_before, old_basis, new_basis)` and `restart_hidden_load_death(lambda_before, old_basis, survivor_basis)` apply the forced finite positive restart maps in explicit support coordinates.
- `kernel_schur_jet_from_coefficients(coefficients)` classifies finite-order support events from the kernel Schur-complement jet.
- `semisimple_event_block(jet, ...)` returns the universal semisimple pole and clock diagnostics for an active event block, including the matrix-support event charge `q_sup = order * dimension / 2` as `clock_log_coefficient` on death-like blocks.
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
- `rank_one_covariance_perturbation(Sigma0, f, visible_dim, epsilon)` diagnoses the exact rank-one covariance/Fisher perturbation under a coordinate visible/hidden split. It reports the hidden-gap increment, the Sherman-Morrison formula increment, rank, channel alignment, and whether the update is one-channel. This is not a generic non-Gaussian law selector.
- `residual_margin_ordering(quadratic_gap, residual_bound)` returns the strict margin certificate for whether a quadratic branch/observer score gap survives residuals bounded by `residual_bound`.
- `variable_precision_affine_hidden_reduction(A, J, D)` evaluates the exact visible action for the affine-hidden sector `exp(-A(v) - 1/2 h^T D(v)h - J(v)^T h)`, including the fibre-volume term `1/2 log det D(v)`.
- `staged_affine_hidden_elimination(A, J, D, eliminate)` performs one exact Schur elimination stage in the affine-hidden sector.
- `weighted_family_frontier_scores(family, B, weights=None, mu=0.0)` evaluates finite weighted-family leakage, visibility, captured curvature, and penalised frontier score for supplied symmetric operators and an already-declared observer plane.
- `exact_branch_hessian(family, B, weights=None, mu=0.0)` returns the fixed-rank exact-branch Hessian contract on `Hom(U,U_perp)` when `B` is already an exact branch.
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
- [`SimpleSpectrumClosureCertificateResult`](../src/nomogeo/types.py) carries the simple-spectrum anchor eigendata, best coordinate eigenspan, minimum cross-block norm, checked subset count, and exact-closure/obstruction flags.
- [`FixedObserverCoordinatesResult`](../src/nomogeo/types.py) carries the exact fixed-observer chart `(phi, hidden_block, coupling)`, the adapted basis used to define it, and conditioning metadata.
- [`ObserverTransitionResult`](../src/nomogeo/types.py) carries the left/right charts together with the exact block transition data and reconstruction residual.
- [`ConnectionCurrentResult`](../src/nomogeo/types.py) carries the chart data, coupling velocity, current `J`, forcing, and the directly computed local `Q`.
- [`IntrinsicLocalQuadraticEnsembleResult`](../src/nomogeo/types.py) carries samplewise intrinsic visible precisions and log-determinant summaries.
- [`CeilingMediatedLocalQuadraticEnsembleResult`](../src/nomogeo/types.py) carries samplewise visible precisions, declared ceilings, hidden loads, clocks, ranks, and descriptive clock summaries.
- [`CoordinateLocalQuadraticEnsembleResult`](../src/nomogeo/types.py) carries the same coordinate-split convenience outputs.
- [`AffineHiddenReductionResult`](../src/nomogeo/types.py) carries the supplied affine-hidden action, coupling, hidden precision, hidden mean, variational action, fibre-volume term, and exact visible action up to an additive constant.
- [`AffineHiddenStageResult`](../src/nomogeo/types.py) carries one Schur elimination stage for the affine-hidden sector.
- [`WeightedFamilyFrontierResult`](../src/nomogeo/types.py) carries weighted-family leakage, visible score, captured curvature, energy-split residual, penalised score, moment operator, projector, and weights.
- [`ExactBranchHessianResult`](../src/nomogeo/types.py) carries the exact-branch Hessian contract, second-variation operator, spectrum, status, off-block norm, branch basis, complement basis, and weights.
- [`SupportStratumTransportResult`](../src/nomogeo/types.py) carries reduced `Lambda`, `Pi`, `A_cpl`, the two transport right-hand sides, clock rate, eigenvalue bounds, and PSD-domain status.
- [`SupportRestartResult`](../src/nomogeo/types.py) carries birth/death restart data, old/new support bases, and the coordinate map between them.
- [`KernelJetResult`](../src/nomogeo/types.py) carries the kernel/gap split, effective Schur-complement coefficients, first nonzero order, leading eigenvalues, and event classification.
- [`SemisimpleEventBlockResult`](../src/nomogeo/types.py) carries finite-order semisimple pole, clock, and desingularisation diagnostics.
- [`LocalCoupledBirthResult`](../src/nomogeo/types.py) carries the local extractor outputs `Phi`, `L`, `R`, `V`, `B`, `beta`, `Q`, observer tensor, `W`, support basis, and `A_cpl`.
- [`SampledIntervalLeakageResult`](../src/nomogeo/types.py) carries sampled leakage, visible score, stationarity residual, projector, weights, and sampled exact-closure status.
- [`IntervalHessianResult`](../src/nomogeo/types.py) carries sampled interval Hessian and spectral-gap rigidity diagnostics.
- [`RankOneCovariancePerturbationResult`](../src/nomogeo/types.py) carries the base/perturbed covariance and precision matrices, direct/formula hidden-gap increments, formula residual, two visible directions, singular values, update rank, and one-channel flag.
- [`ResidualMarginResult`](../src/nomogeo/types.py) carries the quadratic gap, residual bound, required gap, strict margin, worst-case gap, and robustness flags.

## Perturbation And Margin Boundary

The rank-one covariance perturbation diagnostic is exact in its stated theorem domain:

```text
Sigma_epsilon = Sigma_0 + epsilon f f^T
```

with a coordinate visible block. In a white or aligned background, the hidden-gap increment is one-channel. In a generic coloured background, the increment is the difference of two rank-one terms and can have rank two. The API reports this distinction; it does not turn a covariance/Fisher perturbation into a full non-Gaussian law claim.

`residual_margin_ordering` implements only the elementary strict certificate:

```text
quadratic_gap > 2 * residual_bound
```

Equality is not a strict certificate. The function does not estimate the residual; the caller must supply a valid residual bound.

## Affine-Hidden Exact Sector

`variable_precision_affine_hidden_reduction` is an exact special full-law sector, not a generic non-Gaussian marginalisation engine. It assumes the hidden fibre is conditionally Gaussian with supplied precision `D(v)` and affine coupling `J(v)`:

```text
p(v,h) proportional to exp(-A(v) - 1/2 h^T D(v) h - J(v)^T h).
```

The returned visible action is defined up to an additive constant independent of `v`:

```text
S_vis(v) = A(v) + 1/2 log det D(v) - 1/2 J(v)^T D(v)^(-1) J(v).
```

The `variational_action` omits the fibre-volume term. Branch comparisons in the variable-precision sector must use `visible_action`, because the log-determinant term can change the verdict.

`staged_affine_hidden_elimination` exposes one Schur stage. It is exact for the supplied block and is tested by comparing all single-coordinate elimination orders against one-step elimination.

## Weighted Frontier Boundary

`weighted_family_frontier_scores` and `exact_branch_hessian` work on supplied finite families of symmetric local quadratic operators. They do not whiten raw precision perturbations and do not infer a non-Gaussian full law.

`weighted_family_frontier_scores` evaluates:

```text
L_nu(P) = sum_i w_i 1/2 ||[A_i,P]||_F^2
S_nu(P) = sum_i w_i ||P A_i P||_F^2
M_nu = sum_i w_i A_i^2
Tr(P M_nu) = S_nu(P) + L_nu(P)
```

`exact_branch_hessian` requires the supplied observer plane to be an exact branch: all off-blocks must vanish within tolerance. It returns the contract operator

```text
H_mu,U = (1 + mu) C_U^* C_U - G_U.
```

The sign convention is deliberate: positive eigenvalues of `hessian_contract` mean a local maximum of the penalised frontier, because the second variation is `-2 * hessian_contract`. The status is a fixed-rank, fixed-support local statement, not a branch selector and not a probability-support event classifier.

## Ensemble Boundary

`intrinsic_local_quadratic_ensemble` is the observer-intrinsic mode. It is functorial in `(H,C)` and returns only `Phi_C(H_i)` and summaries built from `Phi`.

`ceiling_mediated_local_quadratic_ensemble` is the hidden-load mode. It requires explicit ceilings:

```text
Phi_i = Phi_C(H_i),   Lambda_i = hidden_load(T_i, Phi_i).
```

`coordinate_local_quadratic_ensemble` is a convenience special case. It accepts only a coordinate visible split, because that gives a canonical ceiling block for each sample:

```text
T_i = H_i[:m, :m],   Phi_i = Phi_C(H_i),   Lambda_i = hidden_load(T_i, Phi_i).
```

The returned mean, standard deviation, minimum, and maximum clock are descriptive statistics over exact samplewise outputs. They are not higher cumulants of a full non-Gaussian law, and these functions do not infer mixture mass, support events, or branch probabilities.

## Closure-Adapted Boundary

The current closure-adapted public surface is intentionally exact and narrow.

- `closure_adapted_observer(..., mode="commuting_exact")` is implemented only for pairwise commuting symmetric whitened families.
- noncommuting frontier optimisation is not yet a public API.
- `simple_spectrum_closure_certificate` is only a certificate for the simple-spectrum anchor case. It proves exact common closure exists or is absent among the anchor coordinate eigenspans; it does not search arbitrary noncommuting frontiers when the anchor has spectral degeneracy.
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

