# Theorem Map

This map uses exact theorem / proposition / definition labels from the TeX sources where they exist. Where the paper statement is an unlabeled displayed formula or a code-relevant internal formula rather than a theorem-like environment, the fallback cites the exact equation label or section formula.

## Geometry Core

- `visible_precision`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `[Visible precision]`, label `prop:visible-precision`
  - formulas: `eq:visible-variational`, `eq:tower-law`, `eq:schur-form`
  - quotient confirmation: `Geometric_Consciousness_MASTER.tex`, Definition `def:induced-visible-object`, Proposition `prop:variational`, Proposition `prop:composition`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `canonical_lift`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:lift-projector`
  - theorem using it: `thm:local-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `hidden_projector`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:lift-projector`
  - theorem using it: `thm:local-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `local_visible_calculus`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Local visible calculus]`, label `thm:local-calculus`
  - formulas: `eq:local-expansion`, `eq:VQ-def`, `eq:Q-gram`, `eq:det-split`
  - quotient descent confirmation: `Quotient_Descent_MASTER.tex`, Proposition `prop:hilbert-quotient-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

## Hidden-Load Layer

- `hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Hidden-load geometry]`, label `thm:hidden-load`
  - formulas: `eq:hidden-load-def`, `eq:X-from-Lambda`, `eq:rank-and-clock`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load`, Corollary `cor:hidden-load-rank-volume`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Definition `def:latent-return`, Proposition `prop:hidden-load`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `visible_from_hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `thm:hidden-load`
  - formula: `eq:X-from-Lambda`
  - inverse-theorem role: fixed-ceiling bijection `Lambda <-> X` on active support `S = Ran(T)`
  - domain note: `Lambda` must be positive semidefinite on active support; full-rank ceilings require explicit coordinate interpretation via `lambda_representation='reduced'` or `'ambient'`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_inverse_theorem.py`](../tests/test_inverse_theorem.py)

- `inverse_visible_class`
  - source: same fixed-ceiling inverse theorem as `visible_from_hidden_load`, promoted from `eq:X-from-Lambda` in `thm:hidden-load`
  - role: theorem-facing alias stating the ceiling-conditioned inverse map explicitly
  - non-claim: this is not an inverse of the global fibre map `(H, C) -> Phi_C(H)`
  - tests: [`tests/test_inverse_theorem.py`](../tests/test_inverse_theorem.py)

- `canonical_hidden_realisation`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:canonical-hidden`
  - theorem context: `thm:hidden-load`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

- `minimal_hidden_realisation`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:min-hidden`
  - rank / gauge confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:min-hidden-dim`, Proposition `prop:bridge-hidden-factor-gauge`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

- `transport_hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `[Transport law and determinant clock]`, label `prop:transport`
  - formulas: `eq:transport-law`, `eq:clock-add`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport`
  - domain note: theorem-domain kernel only; both inputs must satisfy `Lambda >= 0` on active support
  - structural note: this is the two-step Gram-shadow transport law, not the associative downstairs object
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `hidden_contraction`
  - source fallback: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport` together with the preceding remark `[The visible law is not the associative object]`
  - formulas: `Pi = (I_S + Lambda)^(-1)` from `thm:hidden-load-transport`, exposed as the canonical factor `K = Pi^(1/2) = (I + Lambda)^(-1/2)`
  - role: exposes the actual associative downstairs object used for long-chain composition
  - tests: [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `load_from_hidden_contraction`
  - source fallback: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport` with `Pi = K^T K`
  - formulas: `Lambda = Pi^(-1) - I`
  - role: converts the associative downstairs object back into hidden-load coordinates
  - tests: [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `clock`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `prop:transport`
  - formula: `eq:clock-def`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Corollary `cor:det-clock`, Theorem `thm:clock-uniqueness`
  - domain note: theorem-domain kernel only; `Lambda` must satisfy `Lambda >= 0` on active support
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- local hidden birth smoke
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Local hidden return as hidden-load birth]`, label `thm:birth`
  - formula: `eq:birth-law`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:local-birth`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

## Finite Bridge

- `dv_bridge`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Donsker-Varadhan bridge]`, label `thm:dv-bridge`
  - formulas: `eq:dv-bridge`, `eq:dv-gram`
  - finite source confirmation: `Finite_Observation_MASTER.tex`, Theorem `[Quadraticity of the DV correction]`, label `thm:exact_quadraticity`
  - quotient descent confirmation: `Quotient_Descent_MASTER.tex`, Proposition `prop:dv-hilbert-lift`
  - tests: [`tests/test_bridge.py`](../tests/test_bridge.py), [`tests/test_release_gate_bridge.py`](../tests/test_release_gate_bridge.py)

- DV quadratic onset and quartic hidden defect smoke
  - source: `Geometric_Observation_MASTER.tex`, Corollary `[Quadratic onset and quartic defect]`, label `cor:quartic-defect`
  - source: `Geometric_Observation_MASTER.tex`, Corollary `[The quartic defect is the first hidden birth]`, label `cor:final`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Corollary `cor:quadratic-dv`
  - tests: [`tests/test_bridge.py`](../tests/test_bridge.py)

## Closure-Adapted Layer

- `whitened_perturbation`
  - source: `Closure_Adapted_Observers.tex`, Theorem `[H-whitened observer normal form]`, label `thm:whitened-normal-form`
  - formulas: `eq:Delta-tilde`, `eq:normal-form-VQ`
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py), [`tests/test_release_gate_adapted.py`](../tests/test_release_gate_adapted.py)

- `observer_from_subspace`
  - source: `Closure_Adapted_Observers.tex`, Theorem `thm:whitened-normal-form`
  - formula: `eq:C-from-B`
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py)

- `closure_scores`
  - source: `Closure_Adapted_Observers.tex`, Definition `[Closure-adapted observer]`, label `def:closure-adapted`
  - formulas: `eq:L-pi-def`, `eq:eta-def`, `eq:energy-split`
  - theorem role: exact leakage `L`, retained visible score `S`, and hidden fraction `eta` on the whitened family
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py), [`tests/test_release_gate_adapted.py`](../tests/test_release_gate_adapted.py)

- `leakage_channels`
  - source: `Closure_Adapted_Observers.tex`, section `Operational consequences`
  - formula: `eq:Qa-op`
  - theorem role: exposes the canonical leakage Gram operator and its singular leakage channels for one perturbation
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py), [`tests/test_release_gate_adapted.py`](../tests/test_release_gate_adapted.py)

- `closure_adapted_observer`
  - source: `Closure_Adapted_Observers.tex`, Theorem `[Commuting family]`, label `thm:commuting-family`
  - source: `Closure_Adapted_Observers.tex`, Corollary `[Optimal exact closure in the commuting case]`, label `cor:optimal-commuting`
  - domain note: current public mode is exact commuting-family synthesis only
  - non-claim: this does not expose the noncommuting frontier optimiser or the `sum_a D_a^2` rule as a final observer
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py), [`tests/test_release_gate_adapted.py`](../tests/test_release_gate_adapted.py)

- `simple_spectrum_closure_certificate`
  - source: `Claude_Notes_Audit_Addendum.tex`, Proposition `[Exact closure can fail generically]`
  - formula: for a simple-spectrum anchor, exact common closure at rank `m` is equivalent to zero cross-blocks on one anchor-coordinate `m`-eigenspan
  - role: exact certificate for common closure or obstruction in the simple-spectrum anchor case
  - non-claim: this is not a noncommuting closure frontier optimiser and rejects degenerate anchors
  - tests: [`tests/test_perturbation.py`](../tests/test_perturbation.py)

- `compare_observers`
  - source fallback: exact score comparison built from `def:closure-adapted`, `eq:L-pi-def`, and `eq:eta-def`
  - role: same-rank exact comparison / inefficiency witness on one task family
  - non-claim: this is not a frontier optimiser; it compares already-specified observers only
  - tests: [`tests/test_adapted.py`](../tests/test_adapted.py), [`tests/test_release_gate_adapted.py`](../tests/test_release_gate_adapted.py)

## Fixed-Observer Connection Layer

- `fixed_observer_coordinates`
  - source: `Connection_Flatness.tex`, fixed-observer factorisation / chart section
  - formulas: observer-fixed adapted basis with exact coordinates `(Phi, R, K)` and reconstruction `H = U^{-T} T_K diag(Phi, R) T_K^T U^{-1}`
  - domain note: the public chart uses an observer-fixed adapted basis; it does not rebuild the basis from `H`
  - tests: [`tests/test_connection.py`](../tests/test_connection.py), [`tests/test_connection_flatness_hardening.py`](../tests/test_connection_flatness_hardening.py)

- `reconstruct_precision_from_fixed_observer_coordinates`
  - source: `Connection_Flatness.tex`, same fixed-observer factorisation formulas
  - role: exact inverse of the public fixed-observer chart when the adapted basis is supplied
  - tests: [`tests/test_connection.py`](../tests/test_connection.py)

- `observer_transition`
  - source: `Connection_Flatness.tex`, observer-to-observer transition law for fixed observers
  - formulas: block transform law for `(Phi, R, K)` under the adapted-basis change `U_1^{-1} U_2`
  - tests: [`tests/test_connection.py`](../tests/test_connection.py), [`tests/test_connection_flatness_hardening.py`](../tests/test_connection_flatness_hardening.py)

- `connection_current`
  - source: `Connection_Flatness.tex`, fixed-observer current / forcing identities
  - formulas: `J = Phi^{-1} dot(K) R`, `Q = Phi J R^{-1} J^T Phi`
  - role: exact current and forcing diagnostics for one tangent direction at fixed observer
  - tests: [`tests/test_connection.py`](../tests/test_connection.py), [`tests/test_connection_flatness_hardening.py`](../tests/test_connection_flatness_hardening.py)

- `forcing_from_current`
  - source: `Connection_Flatness.tex`, current-to-forcing identity
  - formula: `Q = Phi J R^{-1} J^T Phi`
  - tests: [`tests/test_connection.py`](../tests/test_connection.py)

## Local Quadratic Ensemble Layer

- `intrinsic_local_quadratic_ensemble`
  - source: `Claude_Notes_Audit_Addendum.tex`, section `[Local-geometry ensembles]`
  - stronger source: `0.3.1_Technical_Note_1.tex`, Definition `[Intrinsic quotient ensemble]`, Proposition `[Covariance of the intrinsic layer]`
  - formulas: samplewise `Phi_i = Phi_C(H_i)`
  - role: exact observer-intrinsic local quadratic pushforward plus log-determinant summaries
  - non-claim: this does not define hidden load or clock without a ceiling
  - tests: [`tests/test_ensemble.py`](../tests/test_ensemble.py)

- `ceiling_mediated_local_quadratic_ensemble`
  - source: `0.3.1_Technical_Note_1.tex`, Definition `[Ceiling-mediated ensemble]`, Theorem `[No canonical nontrivial ceiling from (H,C) alone]`
  - formulas: samplewise `Phi_i = Phi_C(H_i)`, `Lambda_i = hidden_load(T_i, Phi_i)`, `tau_i = log det(T_i) - log det(Phi_i)`
  - role: exact hidden-load and clock pushforward once a ceiling family is explicitly declared
  - non-claim: hidden load and clock are not intrinsic to `(H,C)` alone
  - tests: [`tests/test_ensemble.py`](../tests/test_ensemble.py)

- `coordinate_local_quadratic_ensemble`
  - source: `0.3.1_Technical_Note_1.tex`, Definition `[Ceiling-mediated ensemble]`, as the coordinate split special case
  - formulas: samplewise coordinate split `T_i = H_i[:m, :m]`, `Phi_i = Phi_C(H_i)`, `Lambda_i = hidden_load(T_i, Phi_i)`
  - role: exact samplewise local quadratic pushforward plus descriptive clock summaries
  - non-claim: this does not compute full-law cumulants, mixture mass, support events, or branch probabilities
  - tests: [`tests/test_ensemble.py`](../tests/test_ensemble.py)

## Affine-Hidden Exact Sector

- `variable_precision_affine_hidden_reduction`
  - source: `0.3.1_Technical_Note_1.tex`, Theorem `[Variable-precision affine-hidden elimination]`, Corollary `[Branch-change criterion]`
  - formulas: `S_vis(v) = A(v) + 1/2 log det D(v) - 1/2 J(v)^T D(v)^(-1) J(v) + const`
  - role: exact hidden Gaussian-fibre elimination with possibly non-Gaussian visible action
  - non-claim: this is not arbitrary non-Gaussian marginalisation and does not determine additive constants independent of `v`
  - tests: [`tests/test_affine.py`](../tests/test_affine.py)

- `staged_affine_hidden_elimination`
  - source: `0.3.1_Technical_Note_1.tex`, Proposition `[Staged elimination tower law]`
  - formulas: Schur complement hidden precision, shifted coupling, and action shift `1/2 log det D_ee - 1/2 J_e^T D_ee^(-1) J_e`
  - role: exact one-stage Schur elimination in the affine-hidden sector
  - tests: [`tests/test_affine.py`](../tests/test_affine.py)

- `affine_hidden_branch_reversal`
  - source: `0.3.2_Technical_Note_1.tex`, Proposition `[Affine-hidden branch reversal]`
  - formula: a variational winner `a` is reversed by branch `b` exactly when `Fib_a - Fib_b > Var_b - Var_a`
  - role: exact finite branch comparison once variational and fibre-volume terms are supplied in the affine-hidden sector
  - non-claim: this is not a generic marginalisation routine or branch probability
  - tests: [`tests/test_affine.py`](../tests/test_affine.py)

- `guarded_fibre_dominance`
  - source: `0.3.2_Technical_Note_1.tex`, Definition `[Guarded fibre-dominance diagnostic]`
  - formulas: centered fibre norm, centered variational norm, and ratio only when the denominator exceeds the declared floor
  - role: diagnostic tuple for declared finite affine-hidden samples
  - non-claim: this is not an invariant naked scalar or law-level branch selector
  - tests: [`tests/test_affine.py`](../tests/test_affine.py)

## Weighted-Family Frontier Layer

- `weighted_family_frontier_scores`
  - source: `0.3.1_Technical_Note_1.tex`, Definition `[Weighted-family leakage and visibility]`, Theorem `[Weighted-family energy split]`
  - formulas: `L_nu(P) = sum_i w_i 1/2 ||[A_i,P]||_F^2`, `S_nu(P) = sum_i w_i ||P A_i P||_F^2`, `Tr(P M_nu) = S_nu(P) + L_nu(P)`
  - role: finite weighted-family local quadratic frontier evaluator
  - non-claim: this is not a noncommuting optimiser or full-law selector
  - tests: [`tests/test_frontier.py`](../tests/test_frontier.py)

- `exact_branch_hessian`
  - source: `0.3.1_Technical_Note_1.tex`, Theorem `[Exact-branch Hessian for weighted families]`, Definition `[Exact-branch Hessian API contract]`, Proposition `[Status semantics]`
  - formulas: `H_mu,U = (1+mu) C_U^* C_U - G_U`; second variation is `-2 H_mu,U`
  - role: exact fixed-rank Hessian diagnostic at an already-exact branch
  - non-claim: this is not a branch optimiser, support-gap classifier, or probability-support event engine
  - tests: [`tests/test_frontier.py`](../tests/test_frontier.py)

## Observation Field Layer

- `pi_from_hidden_load`, `hidden_load_from_pi`
  - source: `Observation_Fields_ADDENDUM.tex`, stratumwise comparison section
  - formulas: `Pi = (I + Lambda)^(-1)`, `Lambda = Pi^(-1) - I`
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `pi_rhs`, `lambda_rhs`, `clock_rate`, `support_stratum_transport`
  - source: `Observation_Fields_UPDATED.tex`, support-stable transport theorem; `Observation_Fields_ADDENDUM.tex`, Theorem `[Stratumwise comparison theorem]`
  - formulas: `dot Pi = -Pi^(1/2) A_cpl Pi^(1/2)`, `dot Lambda = (I + Lambda)^(1/2) A_cpl (I + Lambda)^(1/2)`, `dot tau = Tr(A_cpl)`
  - domain note: the generator may exist when indefinite, but positive hidden-load transport is licensed only when `A_cpl >= 0`
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `restart_hidden_load_birth`, `restart_hidden_load_death`
  - source: `Observation_Fields_UPDATED.tex`, birth/death restart laws; `Observation_Fields_ADDENDUM.tex`, Corollary `[Forced finite restart law]`
  - role: zero extension at birth and survivor compression at death in explicit nested support bases
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `kernel_schur_jet_from_coefficients`, `classify_support_event_from_jet`
  - source: `Observation_Fields_ADDENDUM.tex`, Definition `[Kernel Schur-complement jet]`, Theorem `[Finite-order kernel-event classifier]`
  - formulas: `F(epsilon) = A(epsilon) - B(epsilon)(P + C(epsilon))^(-1)B(epsilon)^T`, recursive `E_n`
  - non-claim: the kernel jet controls leading small-eigenvalue behaviour and near-zero inertia, not exact full-spectrum equality
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `semisimple_event_block`
  - source: `Observation_Fields_ADDENDUM.tex`, Theorem `[Universal finite-order semisimple pole law]`, Corollary `[Finite-order birth, death, and clock laws]`
  - formulas: `A_U = -sigma m/(2s) I + bounded`, death clock coefficient `mr/2`
  - 0.3.1 alignment: `clock_log_coefficient` is the matrix-support event charge `q_sup = m_sup r / 2` on death-like blocks; this is not a probability-support event classifier
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `local_coupled_birth`
  - source: `Observation_Fields_UPDATED.tex`, Proposition `[Exact second covariant visible derivative]`, Theorem `[Support-restricted coupled local birth theorem]`; `Observation_Fields_ADDENDUM.tex`, Theorem `[Module-grade local extractor]`
  - formulas: `beta = -dot C Z`, `W = L^T ddot H L - 2Q - Phi beta R^(-1)B^T - B R^(-1) beta^T Phi`, `A_cpl = -1/2 V_S^(-1/2) W_S V_S^(-1/2)`
  - domain note: requires explicit derivative data; no finite-difference estimation is performed by the API
  - tests: [`tests/test_field.py`](../tests/test_field.py)

- `sampled_interval_leakage`, `sampled_interval_stationarity`, `sampled_interval_closure_check`, `interval_hessian_at_exact_family`
  - source: `Closure_Adapted_Observers.tex`, stationarity equation; `Observation_Fields_ADDENDUM.tex`, Section `[Continuous interval-family leakage calculus]`
  - formulas: sampled versions of `L_I(P)`, `S_I(P)`, `[P, integral [A,[A,P]]]`, and the exact-interval Hessian
  - non-claim: sampled diagnostics are not a public noncommuting optimiser and not a continuum certificate without extra assumptions
  - tests: [`tests/test_field.py`](../tests/test_field.py)

## Quotient Adapter Layer

- `observed_covariance` and Gaussian divergences
  - source context: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:gaussian-dp`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_quotient.py`](../tests/test_release_gate_quotient.py)

- `rank_one_covariance_perturbation`
  - source: `Claude_Notes_Audit_Addendum.tex`, Proposition `[White rank-one covariance sector]`, Proposition `[Coloured background is generally two-rank]`
  - formulas: `Sigma_epsilon = Sigma_0 + epsilon f f^T`; `Delta R = gamma w_V w_V^T - beta u_V u_V^T`
  - role: exact covariance/Fisher perturbation diagnostic for a coordinate visible split; one-channel in the white/aligned sector and generally two-rank in the coloured sector
  - non-claim: this is not a generic non-Gaussian full-law branch or observer selector
  - tests: [`tests/test_perturbation.py`](../tests/test_perturbation.py)

- `rank_k_covariance_perturbation`
  - source: `0.3.2_Technical_Note_1.tex`, Proposition `[Rank-k covariance/Fisher perturbation]`
  - formulas: `Sigma_1 = Sigma_0 + F F^T`; hidden-gap increment is the difference of the full-space and visible-space Woodbury rank-k terms, with rank at most `2k`
  - role: exact covariance/Fisher perturbation diagnostic for a coordinate visible split
  - non-claim: this is not a generic non-Gaussian full-law branch or observer selector
  - tests: [`tests/test_perturbation.py`](../tests/test_perturbation.py)

- `residual_margin_ordering`
  - source: `Branch_Selection_NonGaussian_Findings.tex`, Proposition `[Quadratic verdict with residual margin]`; `Claude_Notes_Audit_Addendum.tex`, section `[Robust branch and observer claims]`
  - formula: strict certificate `quadratic_gap > 2 * residual_bound`
  - role: deterministic margin check once a residual bound is supplied
  - non-claim: the API does not estimate law-level residuals or certify generic non-Gaussian law exactness
  - tests: [`tests/test_perturbation.py`](../tests/test_perturbation.py)

- `declared_ladder_dimension_cost_intervals`
  - source: `0.3.2_Technical_Note_1.tex`, Proposition `[Declared-ladder dimension-cost intervals]`
  - formula: candidate `a` wins exactly on the intersection of pairwise half-lines `s_a - c d_a >= s_b - c d_b`, `c >= 0`
  - role: exact finite phase diagram for a supplied observer ladder under scalar dimension cost
  - non-claim: this is not global Grassmannian optimisation or observer discovery
  - tests: [`tests/test_frontier.py`](../tests/test_frontier.py)

- `general_graph_frontier_hessian`
  - source: `0.3.2_Technical_Note_1.tex`, Theorem `[General graph-chart first variation]`, Theorem `[General graph-chart second variation]`, Proposition `[Frame contract]`
  - formulas: graph-chart gradient `2 sum_i w_i ((2+mu)E_i U_i - mu D_i E_i)` and polarized Hessian from `S_{i,2}` and `T_{i,2}`
  - role: exact declared-observer local quadratic frontier Hessian; recovers `exact_branch_hessian` when all off-blocks vanish
  - non-claim: this is not a global branch optimiser, probability-support classifier, or full-law non-Gaussian selector
  - tests: [`tests/test_frontier.py`](../tests/test_frontier.py)

- `declared_frontier_local_certificate`
  - source: `0.3.2_Technical_Note_1.tex`, Lemma `[Abstract local certificate]`, Proposition `[Conservative graph-chart certificate constant]`
  - formula: sufficient certificate from stationarity residual `eps`, Hessian margin `lambda`, chart radius `rho`, and conservative `L_cert(rho)`
  - role: sufficient local graph-chart max/min certificate for a declared observer
  - non-claim: can be vacuous; when `eps > 0` it certifies a nearby local optimizer, not global optimality of the supplied observer
  - tests: [`tests/test_frontier.py`](../tests/test_frontier.py)

- `observer_collapse_descends`
  - source: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:collapse-asymmetry`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_quotient.py`](../tests/test_release_gate_quotient.py)

- coarse observer misses a real change smoke
  - source: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:toy-example`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py)

## Batch Layer

- `batch_map`
  - source fallback: implementation requirement from `PLAN.txt`, section `Multicore support`
  - no theorem label exists for this engineering layer
  - tests: [`tests/test_batch.py`](../tests/test_batch.py)

