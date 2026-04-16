# Research Frontier Assessment

Date: 2026-04-11

This file records the current unknowns, opportunities, weak points, and low-hanging next steps after the Gaussian/non-Gaussian audit, branch-selection consolidation, weakpoint pass, `papers/claude_notes` audit, and the `0.3.1_Technical_Note_1.tex` planning checkpoint.

## Current Position

The project is past Gaussian toy territory in the following precise sense:

- the core positive-cone quotient geometry is exact without probabilistic Gaussian assumptions once a local Hessian/Fisher/SPD object is supplied;
- Gaussian law mode is the first broad full-law closure sector, not the whole backbone;
- exact special sectors now exist beyond pure Gaussian law mode, especially fixed-shape quadratic tilts, fixed-precision affine-hidden fibres, KL chain-rule descent, total covariance, and white/aligned rank-one covariance perturbations;
- the new weakpoints are structural and sharp, not vague "non-Gaussian might break things" warnings.

The project is not yet a generic full non-Gaussian observation law theory. The highest-risk overreach remains branch/observer verdict promotion from order two to full law without a residual, cumulant, support, or mixture margin.

## Biggest Wins

1. **Affine-hidden exact sector.** Fixed hidden precision makes hidden integration exact for arbitrary visible action `A(v)` and coupling `J(v)`, with staged elimination governed by Schur associativity. This is the strongest full-law sector beyond Gaussian visible laws.

2. **Variable-precision boundary.** The exact `1/2 log det D(v)` fibre-volume term is a clean breakpoint. It explains exactly when affine-hidden minimisation and law marginalisation split.

3. **Rank-one differential-correlation calculus.** In a white/aligned covariance background, the hidden gap is exactly one-channel and signal-aligned. In a coloured background, the update is explicitly two-rank. This preserves the useful Claude neural insight while preventing overstatement.

4. **Generic noncommuting closure obstruction.** Exact closure is common invariant-subspace geometry. Commuting-family wins do not generalise; generic noncommuting pairs can have no exact intermediate-rank closure observer.

5. **Residual margin doctrine.** The right bridge from quadratic verdicts to fuller laws is a certified margin, not a semantic upgrade. If the branch gap is `gamma` and residuals are bounded by `R`, the verdict is robust only when `gamma > 2R`.

6. **Local-geometry ensemble frontier.** For heterogeneous non-Gaussian systems, the next natural layer may be a distribution over local Hessian/Fisher objects. This is exact samplewise and statistically meaningful over the ensemble, but it is not a substitute for cumulants.

## Weak Points

- **Quadratic ties are fragile.** A branch tie at the Hessian layer can be separated by first fibre cumulants or remote mass.
- **Closure is order-two unless explicitly upgraded.** `Q=0` and `eta=0` do not imply full-law marginal closure.
- **Support strata are matrix-rank strata.** They do not classify probability supports, atoms, hard thresholds, or truncation boundaries.
- **Gaussian gluing is covariance gluing.** It does not solve arbitrary law compatibility.
- **Averaging local geometries can erase tails.** `tau(E H_local)` can differ strongly from `E tau(H_local)`.
- **Noncommuting observer synthesis is still open.** The exact commuting mode is solved; the general Pareto frontier is not.

## Low-Hanging Theory Steps

1. **Write the rank-one/two-rank covariance perturbation lemma into the branch/non-Gaussian findings stack.** It is short, exact, and useful for Fisher applications.

2. **Name the entropic fibre-volume sector.** Fixed-precision affine-hidden and variable-precision Gaussian-fibre reductions are both exact, but they are different sectors.

3. **Formalise residual-margin selectors.** State a theorem for score gaps and Hessian gaps with explicit residual norms. This is easy and prevents overclaiming.

4. **Add a closure irreducibility proposition.** The simple-spectrum plus generic second operator proof is short and gives a clean counterweight to commuting examples.

5. **Define local-geometry ensembles.** Start with a definition and samplewise pushforward:

```text
H ~ mu, then Phi_C(H), Lambda_T(H), tau_T(H) are random quadratic observables.
```

Add only what is theorem-grade: Jensen/delta-method warnings, concentration assumptions, and robust observer-selection conditions under ensemble spread.

## Module Scoping

After user approval to proceed, the safe module-facing work is additive and theorem-local:

- rank-one/two-rank covariance perturbation diagnostic; implemented as `rank_one_covariance_perturbation`;
- residual-margin wrapper for branch or observer comparison; implemented as `residual_margin_ordering`;
- closure irreducibility certificate for simple-spectrum families; implemented as `simple_spectrum_closure_certificate`;
- intrinsic ensemble evaluator for samples of local Hessians and arbitrary surjective observers; implemented as `intrinsic_local_quadratic_ensemble`.
- ceiling-mediated ensemble evaluator with explicit supplied ceilings; implemented as `ceiling_mediated_local_quadratic_ensemble`.
- coordinate-split ensemble convenience wrapper; implemented as `coordinate_local_quadratic_ensemble`.

Each should be labelled as quadratic/Fisher or exact special-sector tooling. None should be described as a generic non-Gaussian full-law engine.

Implementation note: the implemented surfaces above are theorem-local and opt-in. They do not widen the Gaussian or quadratic scope of the module. The ensemble surface now follows the 0.3.1 technical split: intrinsic `Phi` summaries for general observers, and hidden-load/clock summaries only in a mediated mode with explicit ceilings. This avoids implying that hidden load or clock are intrinsic to `(H,C)` alone.

## Next Research Frontier

The clean next frontier is not "extend Gaussian formulas to non-Gaussian laws." It is a layered theory:

1. exact quadratic quotient geometry;
2. exact special law sectors where integration/minimisation remains closed;
3. residual-margin transfer from quadratic geometry to fuller laws;
4. higher-order fibre cumulants, fibre-volume terms, support data, and local-geometry ensembles.

The programme's backbone looks right. Gaussian is not the true boundary. But the verdict layer needs discipline: observer choice, collapse, closure, branch selection, and event structure are stable outside Gaussianity only under named hypotheses or margins.

## 0.3.1 Planning Checkpoint

The detailed readiness matrix is in `audit/0_3_1_readiness_matrix.md`, and the
implementation contract plan is in `audit/0_3_1_dev_plan.md`. The new planning
stress script `audit/0_3_1_planning_stress_tests.py` passes and writes
`audit/outputs/0_3_1_planning_stress_tests.json`.

The checkpoint upgrades the roadmap in three ways:

- the intrinsic-vs-ceiling split is now a hard API doctrine, not a wording
  preference;
- the variable-precision affine-hidden reducer was the first post-freeze
  implementation target because it is exact, small, and branch-relevant, and
  it is now implemented;
- weighted-family frontier and branch-Hessian tooling should start as
  evaluator/diagnostic surfaces, not as a noncommuting selector policy.

The following remain research-only until further design: generic fibre-cumulant
branch selectors, probability-support event engines, and Hermitian packet APIs
without an explicit complex dtype policy.
