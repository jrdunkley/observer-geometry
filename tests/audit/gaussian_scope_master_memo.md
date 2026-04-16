# Gaussian Dependence vs Broader Observation Geometry: Master Audit Memo

## Executive Judgement

The project is not merely a Gaussian exact theory. Its mathematical backbone is an exact local quadratic observation geometry on finite-dimensional positive cones: quotient precision, canonical lift, local Taylor calculus, hidden-load coordinates, determinant clock, commutator leakage, closure-adapted observer design, fixed-observer Fisher-Rao geometry, and support-stratum transport are all exact once an SPD Hessian/Fisher/precision object and the required derivatives are supplied.

The project is also not yet an exact full non-Gaussian law theory. Gaussian laws are the first fully closed probabilistic sector because a centred Gaussian law is determined by its first two moments, marginalisation is Schur/quotient precision, Gaussian divergences are spectral functions of covariance/precision, and Gaussian gluing is covariance completion. Outside that sector, higher cumulants, mixtures, hard supports, atoms, and threshold events can change the observational verdict while leaving the entire local quadratic layer unchanged.

The post-audit branch notes sharpen the middle ground. There are now two additional exact non-Gaussian law sectors worth foregrounding:

- **Fixed-shape quadratic tilts** `p_H(x) proportional to exp(-1/2 x^T H x - psi(x))`, where symmetric KL along precision rays is exactly controlled by the variance of the quadratic statistic.
- **Fixed-precision affine-hidden fibres** `p(v,h) proportional to exp(-A(v)-1/2 h^T D h-J(v)^T h)`, where hidden marginalisation is exact and gives `A(v)-1/2 J(v)^T D^{-1}J(v)+const`.

These sectors strengthen the claim that Gaussianity is not the conceptual boundary. They do not make the module a generic non-Gaussian law engine. The sharp weakpoint is visible-dependent hidden precision: if `D=D(v)`, the exact visible action gains `1/2 log det D(v)`, a fibre-volume term that can change branch selection and is not present in the variational quotient.

The correct classification is therefore:

```text
exact local quadratic observation geometry,
with exact Gaussian full-law theory as the first closed law sector,
plus several broader law-level quotient identities and sectors
such as fixed-shape tilts, fixed-precision affine-hidden fibres,
KL chain rule, and total covariance,
but no generic full non-Gaussian law engine.
```

The README clarification is directionally correct and necessary. It does not undersell the module. If anything, it should be read strictly: outside exact Gaussian law mode, `nomogeo` is exact for supplied SPD/Hessian/Fisher geometry, not for arbitrary non-Gaussian pushed-forward laws. The phrase "or Fisher geometry" is correct only when the user really supplies the local Fisher object whose quotient is being studied. It should not be read as "the module computes the observed Fisher information of an arbitrary non-Gaussian observation channel."

## Theorem Layer Classification

The full dense table is in `audit/classification_table.md`. The theorem stack separates into these layers.

### Q: Exact Local Quadratic Geometry

The core exact `Q` layer contains:

- `Phi_C(H)` as constrained quadratic energy / Hilbert quotient.
- Tower law and Schur complement formulas.
- `L_{C,H}` and `P_{C,H}` as minimal lift and hidden projector.
- `V`, `Q`, determinant curvature, and quartic hidden defect as Taylor calculus of the SPD quotient map.
- Hidden-load parametrisation beneath a fixed ceiling, rank/clock identities, minimal block realisations, contraction factors, transport law, and additive determinant clock.
- Closure-adapted whitening, commutator leakage, leakage channels, common invariant-subspace closure, commuting-family observer synthesis, and invariant-flag refinement.
- Fixed-observer global coordinates, Fisher-Rao submersion, metric decomposition, conserved current, and exact visible forcing.
- Observation-field support-stratum transport, local coupled birth extraction, kernel Schur-jet event classification, finite restarts, and intervalwise closure calculus.
- Hilbert quotient root and Hilbert projection root in `Quotient_Descent`.

These results do not require a full Gaussian law. They require a valid SPD/local quadratic object and the stated regularity/support hypotheses.

### G: Exact Full Gaussian Law

The full-law Gaussian layer contains:

- observed Gaussian covariance and precision under linear maps;
- Gaussian KL, reverse KL, Hellinger, Bhattacharyya, and data-processing contraction formulas;
- Gaussian equilibrium-shadow comparison formulas;
- Gaussian noisy-channel formulas;
- common Gaussian gluing, Bell-square Gaussian completion, graph Gaussian gluing, chordal clique gluing, temporal Gaussian triangle, and Gaussian sign-shadow formulas;
- Gaussian canonical-correlation, factor-analysis, and mutual-information interpretations of hidden load/clock where invoked.

These are exact because Gaussian laws are closed under linear maps and determined by covariance. The same formulas are not exact for arbitrary non-Gaussian laws with the same Hessian or covariance.

### GQ: Gaussian Interpretation of a Deeper Quadratic Object

Many objects are algebraically exact as `Q` but have a stronger probabilistic interpretation only in Gaussian mode:

- `Phi_C(H)` is a quotient quadratic form always, but a marginal precision law only for Gaussian precision `H`.
- `tau=log det(I+Lambda)` is an exact determinant clock always, but Gaussian mutual-information/volume language only under Gaussian assumptions.
- hidden rank is exact Gram rank beneath a ceiling; "minimal latent Gaussian factor dimension" is a model-class interpretation.
- `Q` is exact local leakage; "hidden Fisher curvature" is local Fisher language, and "full hidden law leakage" is not licensed.
- Gaussian gluing/canonical correlation interpretations are not covariance-only facts outside Gaussianity.

### A: Approximate / Empirical / Verified-But-Not-Theorem

The papers correctly label several items as verified asymptotic or empirical:

- fast-memory/corridor ladder claims in `Quotient_Observation`;
- large-`M` odd-density class and related scaling statements;
- application-layer examples in `Geometric_Consciousness`;
- `nomodescent` deterministic PSD search when exact constraints do not close the problem;
- sampled interval diagnostics when used as continuum evidence without a quadrature theorem.

These may survive outside Gaussianity in specific models, but they are not general theorem-level survival.

### X/O: Failures and Open Extensions

The sharp non-Gaussian failures are:

- same local Hessian but different full visible law;
- quadratic observer ranking reversed by higher-order visible shape;
- `Q=0` but cubic hidden-visible coupling changes the marginal;
- Gaussian covariance/gluing compatible but non-Gaussian law compatibility undecided or false;
- matrix-rank support strata confused with probability support;
- branch/event structure controlled by remote mixture mass invisible at a mode.
- visible-dependent hidden fibre precision, where the exact law has a `1/2 log det D(v)` branch-active term missing from variational elimination.

The open frontier is a higher-order observation geometry carrying at least cubic tensors, quartic terms, fibre-volume/log-determinant data, support-boundary data, and law-level mixture/branch structure.

## Object Layer Classification

### Visible Precision `Phi_C(H)`

Exact local meaning: yes. It is exactly the quotient metric

```text
y^T Phi_C(H) y = min { x^T H x : Cx=y }.
```

Probabilistic meaning outside Gaussian: conditional. If `H` is a local Hessian/Fisher metric, `Phi_C(H)` is the quotient local metric. It is not generally the precision of the pushed-forward non-Gaussian law. Two laws can share `H` and have different pushed-forward tails/branches.

### Hidden Load `Lambda`

Exact local meaning: yes, after fixing a ceiling `T` and a visible SPD/PSD object `X` with `0<X<=T` on the active support.

Probabilistic meaning outside Gaussian: limited. It is an exact attenuation coordinate and determinant-rank invariant of the supplied matrices. Gaussian latent-factor, canonical-correlation, and mutual-information readings do not transfer generically.

### Determinant Clock `tau`

Exact local meaning: yes, `tau=log det(I+Lambda)` is additive under the hidden-load transport algebra.

Probabilistic meaning outside Gaussian: not as a universal entropy or mutual-information clock. It can still be a useful volume/barrier scalar of the local quadratic geometry.

### First Visible Response `V`

Exact local meaning: yes, `V=L^T Delta L` is the first derivative of the quotient SPD map.

Probabilistic meaning outside Gaussian: only as a local Hessian/Fisher response. Full non-Gaussian visible laws can change with zero `Delta`, as the cubic coupling stress test shows.

### Quartic Leakage `Q`

Exact local meaning: yes, `Q=L^T Delta P H^{-1} Delta L >=0` and its Gram/commutator forms are exact.

Probabilistic meaning outside Gaussian: only local. `Q=0` means second-order closure of the supplied quadratic perturbation, not full-law closure.

### Closure Leakage and Commutator Functionals

Exact local meaning: yes. They are exact norms of off-plane action of symmetric whitened perturbation operators.

Probabilistic meaning outside Gaussian: local Fisher/second-order only. A non-Gaussian task family needs higher-order closure tensors before "closed observer" can mean full-law closed observer.

### Closure-Adapted Observer Design

Exact local meaning: yes for commuting families and specified finite families. The noncommuting stationarity and total-curvature identities are exact local objective identities, but the public module does not expose a noncommuting optimiser.

Probabilistic meaning outside Gaussian: a local task-Hessian observer. Full-law optimal observers can differ.

### Quotient Descent Relations

Exact local meaning: yes for the Hilbert quotient root and linear observer relations.

Broader law meaning: mixed. KL chain rule and total covariance are genuinely non-Gaussian exact in their own settings. Gaussian precision marginalisation is Gaussian. CME/CLE is exact through quadratic test functions with a sharp cubic boundary.

### Branch / Observer Selection Objects

Exact local meaning: yes when branch means local spectral/eigen/support branch of a matrix path.

Probabilistic meaning outside Gaussian: not generally. Remote mixture branches, threshold events, atoms, and support boundaries are law-level data absent from the SPD Hessian alone.

### Support Strata and Restarts

Exact local meaning: yes for rank/support of supplied PSD visible blocks and their Taylor jets.

Probabilistic meaning outside Gaussian: not probability support. The truncated-Gaussian stress test has the same interior Hessian and different hard support.

### Canonical Correlation / Gaussian Gluing

Exact local meaning: the matrix cone and Schur formulas are exact.

Full-law meaning outside Gaussian: no. Gaussian gluing is a covariance-completion problem. Non-Gaussian gluing needs more than covariance/correlation data.

## Stress Test Results

The script `audit/non_gaussian_stress_tests.py` writes `audit/outputs/non_gaussian_stress_tests.json` and `.csv`.
The detailed proof sketches and counterexample derivations are in `audit/technical_derivations.md`.
The raw theorem-like TeX inventory is in `audit/tex_theorem_inventory.md`.

| Case | Quadratic verdict | Full-law verdict | Numeric result |
|---|---|---|---|
| Same local Hessian, heavy-tail law | collapse | visible laws differ | symmetric KL `0.255837`, Hellinger squared `0.0195813` |
| Quadratic ranks `x`, full law ranks `y` | `x` wins because `Delta_H=(0.08,0)` | `y` wins by quartic tail change | KL ratio `y/x ≈ 305.17` |
| Cubic hidden-visible coupling | no order-two response/leakage | visible marginal changes | symmetric KL `0.0139792`, mean shift `-0.087432` |
| Remote symmetric branches | central Hessian collapse | threshold event changes | `P(|x|>3)` rises from `0.0026998` to `0.0525648` |
| Same interior Hessian, different support | same `H=1` | support event differs; one KL infinite | Gaussian mass below `-2` is `0.0227501`, truncated mass `0` |
| Variable hidden precision | variational affine-hidden quotient is flat | exact marginal selects by `1/2 log det D(v)` | missing term range `0.935`; hostile exact argmin at boundary |
| Quadratic branch tie, fibre cumulant perturbation | branch gap `0` | visible law separates immediately | symmetric KL `0.0027934`, mean shift `-0.0350396` |
| Student-t scale perturbations | strong perturbation outranks weak | same ranking | KL ratio `≈11.62` |
| Small shared quartic family | strong Hessian perturbation outranks weak | same ranking | KL ratio `≈12.47` |

The pathology details are in `audit/pathology_library.md`; survival details are in `audit/survival_library.md`.

## Failure Modes

1. **False collapse by matched Hessian.** A local Hessian can agree while tails, branches, and thresholds differ.
2. **Observer ranking reversal.** A small Hessian difference can be dominated by a zero-Hessian higher-order visible-law change.
3. **Quadratic closure failure at full law.** `Q=0` eliminates the quadratic off-plane term, but cubic hidden-visible terms can alter the visible marginal after hidden integration.
4. **Gaussian gluing overreach.** Covariance/correlation completion is not law completion outside Gaussianity.
5. **Support confusion.** Matrix support/rank events are not hard probability supports, atoms, or truncation boundaries.
6. **Branch selection overreach.** Spectral branches of local forms do not see remote mixture branches.
7. **One-channel overreach.** Rank-one differential-correlation geometry is one-channel in a white/aligned covariance background, but in a coloured background the hidden-gap update is generally a two-rank signed update.
8. **Commuting-closure overreach.** Exact intermediate-rank closure can be absent for generic noncommuting symmetric perturbation families.
9. **Averaged-Hessian overreach.** Tail-sensitive distributions of local clocks and hidden loads can be lost by replacing a local-geometry ensemble with a single averaged Hessian.

## What Survives Cleanly

The structural geometric content survives almost effortlessly:

- quotient precision as a constrained quadratic form;
- lift/projector and tower composition;
- all local derivative identities of `Phi_C`;
- determinant curvature split;
- hidden-load fixed-ceiling cone, rank, and determinant clock;
- hidden-load transport on fixed support;
- closure-adapted commutator calculus;
- fixed-observer Fisher-Rao geometry;
- support-stratum ODE diagnostics and event-jet classification;
- Hilbert quotient root;
- regular Fisher quotient descent;
- KL chain rule and total covariance in their separate law-level settings.
- fixed-shape quadratic-tilt KL laws;
- fixed-precision affine-hidden marginalisation and tower elimination.
- white/aligned rank-one covariance perturbation calculus, including the exact one-channel hidden-gap formula and Fisher ceiling.
- residual-margin selection theorems, provided the residual bound is an actual assumption or estimate.

A clean broader restatement is:

> Given a finite-dimensional positive quadratic object `H` and a surjective linear observation `C`, observation induces the quotient quadratic form `Phi_C(H)`. The induced form, its local Taylor calculus, and fixed-ceiling hidden-load attenuation geometry are exact positive-cone geometry. Gaussian probability is a closed law sector in which this geometry also determines full observed laws and divergences.

## What Remains Gaussian Dependent

The following should keep Gaussian labels:

- exact marginal precision law from `H`;
- closed-form Gaussian divergences and data-processing values computed from covariances;
- Gaussian sign-shadow arcsine formulas;
- Bell/temporal/graph Gaussian gluing;
- canonical-correlation and mutual-information readings of hidden load;
- factor-analysis hidden realisation as a law-level latent model;
- Gaussian noisy-channel shadow formulas;
- any claim that covariance/Hessian fully determines visible law.

## Current Labels

Recommended doctrine:

1. **Exact Gaussian Law Mode.** Inputs are Gaussian precisions/covariances or explicitly Gaussian pair/context laws. Output claims can be full law claims.
2. **Exact Local Quadratic Mode.** Inputs are supplied Hessian/Fisher/SPD matrices and derivatives. Output claims are exact about quotient quadratic geometry, not full law.
3. **Higher-Order Extension Mode.** Inputs include third and higher local tensors, cumulants, mixture/support data, or law families. This is mostly open, with isolated exact non-Gaussian identities from `Quotient_Descent` as models.
4. **Exact Special Law Sectors.** Inputs certify a non-Gaussian family with an additional closed structure, such as fixed-shape quadratic tilt or fixed-precision affine-hidden fibre. Output claims are exact only inside that named sector.
5. **Local-Geometry Ensemble Mode.** Inputs are distributions or samples of local Hessian/Fisher objects. Output claims are exact samplewise at the quadratic layer and statistical only over the ensemble. This can capture heterogeneity of `Phi`, `Lambda`, and `tau`, but it is still not a full-law cumulant theory.

This is not merely relabelling. The stress tests show that higher-order terms can reverse observer choice, invalidate closure as a full-law statement, create branch events, and produce support singularities invisible to the Hessian. A higher-order theory would need new objects, not just softer wording.

## Next Theory Frontier

The next frontier is a higher-order observation geometry. Minimal ingredients:

- a projected cubic tensor calculus for visible skew and hidden-visible cubic couplings;
- quartic/tail tensors capable of distinguishing matched-Hessian laws;
- law-level branch and mixture stratification;
- probability-support boundary data separate from matrix-rank support;
- non-Gaussian gluing criteria that use more than covariance;
- full-law observer objective frontiers combining quadratic Fisher curvature with higher cumulant divergences;
- continuity/restart laws for support and atoms in actual distributions.
- ensemble laws for local quadratic geometries, including concentration, tails, and stability of observer selection under local heterogeneity.
- closure-frontier algorithms for noncommuting families, with explicit irreducibility certificates and residual margins.

Gaussian is therefore not the true boundary of the programme. It is the first fully solved full-law sector and a useful closed laboratory. It becomes misleading only if Gaussian exactness is allowed to stand in for full non-Gaussian law exactness.

## Final Answer to the Primary Question

The project's core is closest to option 2:

> an exact local quadratic observation geometry with Gaussian laws as its first fully closed sector.

There is "something in between" only in the sense that the paper stack also contains exact non-Gaussian law-level quotient identities outside the module, especially KL chain rule and total covariance. But the `nomogeo` kernel itself should not be read as a full non-Gaussian law engine.
