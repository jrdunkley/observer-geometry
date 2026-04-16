# Gaussian Scope Classification Table

Tags:

- `G`: exact full Gaussian law
- `Q`: exact local quadratic geometry
- `GQ`: Gaussian-dependent interpretation of a quadratic object
- `A`: approximate or empirical survival outside Gaussianity
- `X`: genuine failure outside Gaussianity
- `O`: open extension target

## Module Primitives

| Object / primitive | Tag | Classification reason |
|---|---:|---|
| `visible_precision(H, C) = (C H^{-1} C^T)^{-1}` | Q / GQ | Exact as constrained quadratic energy and Hilbert quotient. Exact marginal precision only for centred Gaussian laws, or for a user-supplied local Hessian/Fisher metric read as a quotient metric. |
| `visible_geometry`, `canonical_lift`, `hidden_projector` | Q | Exact minimal-lift and hidden-fibre projector algebra for a supplied SPD form. No full-law claim. |
| `local_visible_calculus(H, C, Delta)` | Q / GQ | Exact Taylor calculus of the SPD quotient map. It is the exact local Hessian/Fisher response if `H(t)` is the supplied local metric path. It is not the full non-Gaussian marginal-law response. |
| `V` first visible response | Q | Exact first derivative of `Phi_C`; can be a local Fisher/Hessian response. Full visible law may have unchanged `V` but changed skew/tails. |
| `Q` quartic defect | Q / GQ | Exact second nonlinear quotient defect and positive Gram form. Its Fisher/leakage reading is local quadratic; Gaussian full-law "first hidden birth" is not a generic full-law statement. |
| determinant split `D^2(-log det Phi)=||V||^2+2Tr Q` | Q / GQ | Exact log-determinant curvature of the quotient SPD map; Gaussian information-volume readings are Gaussian/local-quadratic only. |
| `whitened_perturbation`, `observer_from_subspace` | Q | Exact whitening normal form for SPD geometry. |
| `closure_scores`, leakage `L`, visible score `S`, hidden fraction `eta` | Q / GQ | Exact commutator leakage and retained local curvature for a supplied whitened perturbation family. Exact Fisher-tightness only at the quadratic/Fisher level. |
| `leakage_channels` | Q | Exact singular channel decomposition of the quadratic coupling map. Full non-Gaussian leakage channels require higher cumulant/tensor data. |
| `closure_adapted_observer(..., commuting_exact)` | Q | Exact common-eigenbasis solution for commuting whitened Hessian/Fisher perturbations. Not a full-law observer optimiser. |
| `compare_observers` | Q | Exact same-rank comparison of specified observers under the local score pair `(L,S)`. Full-law dominance can reverse. |
| `fixed_observer_coordinates`, reconstruction, `observer_transition` | Q | Exact fixed-observer chart, inverse chart, and transition law on the SPD cone. |
| `connection_current`, `forcing_from_current` | Q | Exact fixed-observer current/forcing identity for one tangent direction in SPD/Fisher-Rao geometry. |
| `hidden_load(T, X)`, `visible_from_hidden_load`, `inverse_visible_class` | Q / GQ | Exact fixed-ceiling bijection between `X` beneath `T` and a positive load `Lambda`. Gaussian latent-factor, canonical-correlation, or mutual-information readings are Gaussian. |
| `canonical_hidden_realisation`, `minimal_hidden_realisation` | Q / GQ | Exact positive block-matrix realisations and Gram rank. Minimal hidden random factor interpretation is Gaussian/model-class dependent. |
| `hidden_contraction`, `load_from_hidden_contraction`, `transport_hidden_load` | Q | Exact contraction/transport algebra on a fixed support stratum. Not generic composition of non-Gaussian law effects. |
| `clock(Lambda)=log det(I+Lambda)` | Q / GQ | Exact additive determinant clock for hidden-load transport. Equals Gaussian volume/mutual-information quantities only in Gaussian sector. |
| `dv_bridge` | Q | Exact local Donsker-Varadhan Hessian bridge in the finite Markov quadratic expansion. It is not a full path-law Gaussian theorem. |
| `observed_covariance` | G / Q | Exact covariance pushforward if `H^{-1}` is a covariance matrix; exact Gaussian observed covariance in the law-level adapter. |
| `gaussian_forward_kl`, `gaussian_reverse_kl`, `gaussian_hellinger_squared`, `gaussian_data_processing_contraction` | G | Exact only for centred Gaussian covariance laws. For non-Gaussian laws with the same covariance they are surrogate quantities. |
| `observer_collapse_descends` | Q | Equality of quadratic visible objects descends under coarsening. The analogous full-law equality descends only if equality is equality of the actual richer pushed-forward laws. |
| `pi_from_hidden_load`, `hidden_load_from_pi` | Q | Exact reduced positive-cone coordinate conversion. |
| `pi_rhs`, `lambda_rhs`, `clock_rate`, `support_stratum_transport` | Q | Exact support-stable hidden-load transport diagnostics for supplied `A_cpl`. Positive transport is licensed only when `A_cpl >= 0`. |
| `comparison_envelope_bounds` | Q | Exact Loewner envelope for the support-stratum ODE assumptions. |
| `restart_hidden_load_birth`, `restart_hidden_load_death` | Q | Exact finite positive restart maps for rank/support of the quadratic visible block in explicit nested bases. Not probability-support restart laws. |
| `kernel_schur_jet_from_coefficients`, `classify_support_event_from_jet` | Q | Exact leading small-eigenvalue event classifier from supplied Taylor coefficients. It does not classify full law support or multimodal branch events. |
| `semisimple_event_block` | Q | Exact finite-order pole/clock law for a semisimple matrix event block. |
| `local_coupled_birth` | Q | Exact extractor for the coupled local birth generator from explicit derivatives. It is not finite-difference estimation and not a global field simulator. |
| `sampled_interval_leakage`, `sampled_interval_stationarity`, `sampled_interval_closure_check`, `interval_hessian_at_exact_family` | Q / A | Exact for the supplied finite sample family; continuum interval closure requires analytic/quadrature assumptions. |
| `weighted_family_frontier_scores` | Q | Implemented exact finite/weighted family identity `Tr(P M_nu)=S_nu(P)+L_nu(P)` and penalised frontier score. It is a local quadratic observer-frontier evaluator, not a generic full-law selector or optimiser. |
| `exact_branch_hessian` | Q | Implemented exact local graph-Hessian contract for already-exact branches, including noncommuting internal blocks. A positive returned contract operator means a local maximum of the penalised frontier because the second variation is negative. Not a full-law branch optimiser. |
| fixed-precision affine-hidden reducer (research target) | non-Gaussian exact special sector | Exact full-law marginalisation for `p(v,h) proportional to exp(-A(v)-1/2 h^T D h-J(v)^T h)` with fixed `D`. Staged elimination is exact by Schur associativity. |
| `variable_precision_affine_hidden_reduction` | non-Gaussian exact special sector | Implemented exact hidden Gaussian-fibre reduction including the entropic term `1/2 log det D(v)`. Stress tests show this term can reverse a flat variational branch verdict. It must not be collapsed into variational affine-hidden elimination. |
| `staged_affine_hidden_elimination` | non-Gaussian exact special sector | Implemented exact one-stage Schur elimination helper for the affine-hidden tower law. |
| fibre-cumulant diagnostics (research target) | O | Correct higher-order law object for hidden perturbations, but not a current exact generic module engine. |
| `rank_one_covariance_perturbation` | Q / non-Gaussian exact special sector | Implemented theorem-local covariance/Fisher diagnostic. In a white/aligned covariance background, a rank-one differential-correlation update gives an exact one-channel hidden-gap formula. In a coloured background the update is generally the difference of two rank-one terms and has rank at most two, so the one-channel interpretation is not generic. |
| `residual_margin_ordering` | Q / O | Implemented strict margin certificate for bounded score residuals: a quadratic gap is certified only when `quadratic_gap > 2 * residual_bound`. It does not estimate residuals or certify generic full-law exactness. |
| `simple_spectrum_closure_certificate` | Q | Implemented exact closure/obstruction certificate for the simple-spectrum anchor case. It checks anchor-coordinate eigenspans; it is not a noncommuting closure optimiser and rejects degenerate anchors. |
| `intrinsic_local_quadratic_ensemble` | Q / A | Implemented general-observer intrinsic ensemble evaluator. It pushes local Hessian/Fisher samples through exact samplewise `Phi_C(H_i)` and reports functorial visible summaries only. It does not expose hidden load or clock without ceiling data. |
| `ceiling_mediated_local_quadratic_ensemble` | Q / A | Implemented declared-ceiling ensemble evaluator. It computes hidden load, rank, and clock only after explicit ceilings `T_i` are supplied with `Phi_i <= T_i`. |
| `coordinate_local_quadratic_ensemble` | Q / A | Implemented coordinate-split convenience wrapper. It uses the coordinate ceiling block `H_i[:m,:m]`, then reports exact samplewise hidden-load and clock summaries. It is not a full non-Gaussian law engine. |
| `batch_map`, validation helpers, metadata | Q / engineering | Deterministic engineering helpers; no law-level claim. |

## Paper Theorem Clusters

| Paper theorem / identity cluster | Tag | Classification reason |
|---|---:|---|
| `Finite_Observation`: reduced conductance identity, DV quadraticity, positivity/rank, detailed-balance characterisation | Q | Exact local quadratic Hessian theorem for the DV expansion of finite irreducible Markov chains. |
| `Finite_Observation`: exact block formula, fast-block coercive bound | Q | Exact algebra inside the Hessian/weighted-signal geometry. |
| `Finite_Observation`: projector formula, weighted Ky Fan envelope, retention/hidden fractions, stable-rank obstruction, alignment theorem | Q | Exact spectral geometry of a positive quadratic correction under orthogonal observers. |
| `Finite_Observation`: shadow operator and backbone factorisation | Q / GQ | Exact quadratic shadow object; becomes exact statistical distinguishability only after Gaussian observation is imposed. |
| `Finite_Observation`: exact Gaussian comparison formulas, Gaussian equilibrium-shadow criterion, noisy-channel shadow formulas | G | Exact centred Gaussian laws determined by precision/covariance. Not valid for arbitrary laws sharing the same Hessian or covariance. |
| `Finite_Observation`: cumulative envelope spectroscopy and observability dominance | Q | Exact spectral readout of the quadratic correction. |
| `Finite_Observation`: homogeneous ring / large-system claims | A / O | Claimed as verified/asymptotic or model-specific where the paper labels them so; not generic theorem-grade non-Gaussian geometry. |
| `Geometric_Observation`: visible precision, variational law, tower law, Schur form | Q / GQ | Exact positive-cone quotient geometry; exact Gaussian marginal precision when `H` is Gaussian precision. |
| `Geometric_Observation`: local visible calculus, determinant split, quadratic onset/quartic defect | Q / GQ | Exact SPD Taylor calculus; Gaussian law interpretation only in Gaussian/local-quadratic sector. |
| `Geometric_Observation`: hidden-load geometry, minimal realisation, rank/clock | Q / GQ | Exact fixed-ceiling positive-cone parametrisation; Gaussian latent-factor and mutual-information readings are extra. |
| `Geometric_Observation`: transport law and determinant clock | Q | Exact hidden-load/contraction transport on fixed support. |
| `Geometric_Observation`: local hidden return and "quartic defect is first hidden birth" | Q / GQ | Exact for local quotient-Hessian branch after ceiling normalisation. Full non-Gaussian hidden return can begin at higher/lower non-quadratic order. |
| `Closure_Adapted`: whitened normal form, closure criterion, commutator identity | Q / GQ | Exact local quadratic/Fisher observer geometry. |
| `Closure_Adapted`: single-operator and commuting-family exact solutions | Q | Exact common invariant-subspace/eigenbasis solution for supplied symmetric perturbations. |
| `Closure_Adapted`: insufficiency of pure `eta`, stationarity equation, total captured curvature | Q / O | Exact local objective identities; noncommuting frontier optimiser remains open/publicly unimplemented. |
| `Closure_Adapted`: quotient-family refinement and flag obstruction | Q | Exact invariant-flag criterion for quadratic closure towers. |
| `Closure_Adapted`: order upgrade on adapted branches | Q / GQ | Exact matrix-order improvement for smooth spectral functionals of the local shadow; Gaussian divergence use is Gaussian/local. |
| `Connection_Flatness`: global trivialisation, fixed-observer chart | Q | Exact SPD-cone diffeomorphism for fixed observer. |
| `Connection_Flatness`: Riemannian submersion and flat fixed-observer connection | Q | Exact Fisher-Rao geometry on `SPD`; not a theorem about arbitrary statistical manifolds without supplying the metric. |
| `Connection_Flatness`: spectral/trace identities, current/forcing, visible EOM, fixed-K slices, metric decomposition | Q | Exact SPD/Fisher-Rao differential geometry. |
| `Connection_Flatness`: hidden-mediated commutator as curvature | Q / O | The defect tensor identity is exact; its status as a curvature tensor is explicitly open. |
| `Observation_Fields_UPDATED`: exact pullback action, conserved current, split-frame theorem, gauge laws | Q | Exact field calculus on finite real SPD/symmetric coordinates. |
| `Observation_Fields_UPDATED`: coupled local birth, support-stable transport, support bifurcations, geodesic energy, canonical realisation | Q / O | Exact local/stratumwise. Global smooth PDE across support changes and universal `Q` replacement are explicitly not claimed. |
| `Observation_Fields_ADDENDUM`: module-grade extractor and invariance | Q | Exact extractor from supplied derivatives. |
| `Observation_Fields_ADDENDUM`: stratumwise comparison, bounded-generator obstruction, forced restarts | Q | Exact support-stratum theorem and sharp singular boundary. |
| `Observation_Fields_ADDENDUM`: kernel Schur-complement jet, finite-order classifier, semisimple pole law | Q | Exact leading small-eigenvalue matrix event theory. |
| `Observation_Fields_ADDENDUM`: intervalwise closure/refinement/leakage/stationarity/Hessian/sampled consistency | Q / A | Exact continuum when assumptions hold; sampled APIs certify samples only. |
| `Observation_Fields_ADDENDUM`: Hermitian-positive packet lift | Q / O | Exact finite Hermitian positive-cone bridge; not raw single-amplitude field ontology. |
| `Quotient_Observation`: variational characterisation, positive-cone calculus, Frechet derivative, tangent ceiling, composition, shadow transport, Schur gap | Q / GQ | Exact quotient-positive geometry; Gaussian marginal precision if read as Gaussian law. |
| `Quotient_Observation`: visible algebraic object, signed quartic bridge, sextic divergence robustness | Q / GQ | Exact quotient-Hessian/Schur geometry; divergence robustness becomes Gaussian/local spectral when divergences are Gaussian. |
| `Quotient_Observation`: latent contraction, support theorem, interval theorem, hidden load, monotonicity, minimal hidden dimension, transport, clock uniqueness | Q / GQ | Exact fixed-ceiling positive-cone geometry; Gaussian hidden-variable interpretation is not generic. |
| `Quotient_Observation`: conservation split, paired-spectrum criterion, hidden-component mediation, graph-distance decay, corridor locality | Q / A | Exact where stated as Schur/locality algebra; verified/asymptotic ladder claims remain `A/O` as labelled. |
| `Quotient_Observation`: Bell-square collapse, graph Gaussian gluing, chordal clique gluing, temporal triangle, sign-shadow tetrahedron | G | Full law statements for centred Gaussian completion/sign-shadow formulas. Outside Gaussian, covariance/correlation compatibility is not law compatibility. |
| `Quotient_Observation`: strict refinement over Bell-locality / fibre crossing | G / X | Exact Gaussian counterexamples showing lower-dimensional shadows do not determine full pair-law gluing. They also illustrate non-Gaussian danger: shadows are not full laws. |
| `Quotient_Observation`: Bisker parity collapse, two-kernel timing collapse, finite CTMC completion | Q / non-Gaussian exact finite-law | Exact finite Markov/phase-type construction, but not part of the generic `nomogeo` Gaussian kernel. |
| `Quotient_Observation`: structural identifiability bridge and Fisher refinement | Q | Exact regular-quotient Fisher descent for statistical models. This is broader than Gaussian but local Fisher. |
| `Quotient_Observation`: bridge hidden-load split and Gram gauge | Q / GQ | Exact fixed-ceiling split; hidden factor gauge has Gaussian factor-analysis resonance but is matrix geometry. |
| `Quotient_Descent`: QDP definitions and synthesis | Q / O | Broad organising doctrine; exact only instance-by-instance. |
| `Quotient_Descent` Instance A structural Fisher descent | Q | Exact for regular statistical quotients, not Gaussian-specific. |
| `Quotient_Descent` Instance B CME/CLE quadratic visible shadow | Q / X | Exact quadratic shadow; cubic remainder is the sharp non-quadratic boundary and sign-indefinite. |
| `Quotient_Descent` Instance C dissipative hidden elimination | Q / O | Exact in stated linear/resolvent cases; nonlinear memory embedding remains frontier. |
| `Quotient_Descent` Instance D precision descent under marginalisation | G / Q | Exact Gaussian graphical marginal precision; Schur complement algebra is Q but law-level marginal precision is Gaussian. |
| `Quotient_Descent` Instance E observability Gramian, F Dirichlet-to-Neumann | Q | Exact quadratic/operator energy descent. |
| `Quotient_Descent` Instance G KL descent under marginalisation | non-Gaussian exact / Q | Exact chain rule/data processing for all laws with KL defined; its second-order shadow is Fisher/quadratic. |
| `Quotient_Descent` Instance H total covariance | non-Gaussian exact / Q | Exact law of total covariance for finite-second-moment laws; Gaussian duality to precision descent is special. |
| `Quotient_Descent` Hilbert projection root and Hilbert quotient root | Q | Exact abstract quadratic geometry; strongest evidence that the backbone is not merely Gaussian. |
| `Geometric_Consciousness`: latent/access definitions, variational object, quotient geometry, tower composition | Q / GQ | Exact positive visible-object geometry; consciousness language is interpretation. |
| `Geometric_Consciousness`: hidden-load parametrisation, local response, onset, no interior reversal | Q / GQ | Exact local quadratic/hidden-load statements; empirical reversals need extra ingredients. |
| `Geometric_Consciousness`: collapse asymmetry | Q | Exact for visible quadratic objects; full-law equality under richer observation descends only if actual pushed-forward laws are equal. |
| `Geometric_Consciousness`: Gaussian data processing and observer gain | G | Exact for observed centred Gaussian distributions. |
| `Geometric_Consciousness`: empirical worked examples and predictions | A / O | Application-layer, not theorem-grade module claims. |
| `Branch_Selection_NonGaussian_Findings`: leakage convention, branch parity, selector Hessian, exchange mode, dominance crossing | Q | Exact local branch calculus for supplied symmetric Hessian/Fisher perturbation families. Gaussian law promotion is separate. |
| `Branch_Selection_NonGaussian_Findings`: fixed-shape quadratic tilt theorem | non-Gaussian exact special sector | Exact symmetric-KL law for a fixed non-Gaussian shape under quadratic precision tilts, assuming integrability/moment hypotheses. |
| `Branch_Selection_NonGaussian_Findings`: fixed-precision affine-hidden sector and tower law | non-Gaussian exact special sector | Exact full-law hidden marginalisation and staged elimination for arbitrary visible action with conditionally Gaussian fixed-precision hidden fibre. |
| `Branch_Selection_NonGaussian_Findings`: variational fourth-order quotient, entropic derivative calculus, Laplace bridge | O / A | Correct next-layer theorem targets and asymptotic bridge; not generic closed full-law engine. |
| `Branch_Selection_NonGaussian_Weakpoints`: variable hidden precision boundary | X / O | Shows the exact fibre-volume term `1/2 log det D(v)` can change verdicts and is absent from variational elimination. |
| `Claude_Notes_Audit_Addendum`: rank-one/two-rank differential-correlation calculus | Q / non-Gaussian exact special sector | Exact covariance/Fisher perturbation calculus. The one-channel result is white/aligned; the coloured-background update is generally rank two. |
| `Claude_Notes_Audit_Addendum`: noncommuting closure obstruction and residual margin doctrine | Q / O | Exact local quadratic invariant-subspace obstruction and branch-margin theorem. Generic full-law branch selection still requires residual/cumulant bounds. |
| `Claude_Notes_Audit_Addendum`: local-geometry ensemble proposal | O / A | Good next-layer research target for heterogeneous non-Gaussian systems; not a derivation of full-law tails or application exponents. |
| `nomodescent` core | G / Q / A | Exact finite linear/Gaussian/quadratic observer relation tests; deterministic PSD search is audited approximate where enabled. |
| `evidence` layer | A | Encoding/audit layer with exact/parsed/inferred/ambiguous status, not a source of new theorem-grade geometry. |

## Cross-Cutting Breakpoints

| Breakpoint | Tag | Classification reason |
|---|---:|---|
| Same local Hessian, different tails/skew/mixture branches | X | Quadratic objects agree exactly while full visible laws differ. |
| Quadratic observer ranking vs full-law ranking | X | A higher-order visible change can dominate a small Hessian perturbation. |
| Quadratic closure vs cubic hidden-visible coupling | X / O | `Q=0` only closes order-two geometry; full-law marginal can change through cubic/higher terms. |
| Gaussian common gluing vs non-Gaussian gluing | G / O | Covariance/correlation completion solves Gaussian gluing, not arbitrary law compatibility. |
| Probability support boundaries vs matrix support strata | X / O | Matrix rank events do not classify hard support, atoms, truncation, or threshold events in arbitrary laws. |
| Visible-dependent hidden fibre precision | X / O | Exact marginal law contains a branch-active log-determinant term that is invisible to fixed-precision affine-hidden variational reduction. |
| Quadratic branch tie vs fibre-cumulant branch separation | X / O | The quadratic branch gap can be zero while the visible law separates at first cumulant order. |
| White one-channel differential correlation vs coloured two-rank update | Q / X | A rank-one covariance perturbation gives one visible hidden-gap channel only when the full and visible precision-weighted signal directions align. Generic coloured backgrounds split it into two directions. |
| Commuting exact closure vs generic noncommuting closure | Q / X | Exact closure at intermediate rank is easy in a commuting family and generically absent for noncommuting symmetric families. |
| Single averaged Hessian vs local-geometry ensemble tails | O / X | A single representative quadratic object can miss tail-sensitive variation in the distribution of local clocks or hidden loads. |
| Elliptical and weak higher-order perturbations | A | Quadratic rankings often survive empirically when higher cumulants are controlled or shape family is fixed, but this is not theorem-grade full-law exactness. |
