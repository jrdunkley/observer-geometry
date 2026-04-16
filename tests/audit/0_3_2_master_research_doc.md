# 0.3.2 Master Research Doc

Date: 2026-04-11

Status: research planning artifact. This document digests `papers/claude_notes`
as an external/quarantined source. It is not a paper update and not a module
implementation plan. The module source should remain closed until the theory
targets below are separately approved.

Companion probe:

- `audit/0_3_2_research_probe.py`
- `audit/outputs/0_3_2_research_probe.json`
- `audit/0_3_2_second_pass_probe.py`
- `audit/outputs/0_3_2_second_pass_probe.json`
- `audit/0_3_2_near_branch_variation_check.py`
- `audit/outputs/0_3_2_near_branch_variation_check.json`
- theorem research note: `audit/0_3_2_theory_research_note.md`
- frontier push note: `audit/0_3_2_frontier_push_note.md`
- weighted-frontier certificate note:
  `audit/0_3_2_weighted_frontier_certificate_note.md`
- projector-bounds refinement:
  `audit/0_3_2_projector_bounds_note.md`
- non-exact stationary frontier note:
  `audit/0_3_2_non_exact_stationary_frontier_note.md`
- non-exact stationary frontier check:
  `audit/0_3_2_non_exact_stationary_frontier_check.py`
- `audit/outputs/0_3_2_non_exact_stationary_frontier_check.json`
- general graph Hessian note:
  `audit/0_3_2_general_graph_hessian_note.md`
- general graph Hessian check:
  `audit/0_3_2_general_graph_hessian_check.py`
- `audit/outputs/0_3_2_general_graph_hessian_check.json`
- graph Hessian invariance check:
  `audit/0_3_2_graph_hessian_invariance_check.py`
- `audit/outputs/0_3_2_graph_hessian_invariance_check.json`
- declared frontier certificate note:
  `audit/0_3_2_declared_frontier_certificate_note.md`
- declared frontier certificate check:
  `audit/0_3_2_declared_frontier_certificate_check.py`
- `audit/outputs/0_3_2_declared_frontier_certificate_check.json`
- shelved feature recovery note:
  `audit/0_3_2_shelved_feature_recovery_note.md`
- shelved feature recovery check:
  `audit/0_3_2_shelved_feature_recovery_check.py`
- `audit/outputs/0_3_2_shelved_feature_recovery_check.json`

## Executive Judgement

The updated Claude notes are useful. They show that the 0.3.1 tools are not
merely decorative: the affine-hidden reducer, staged elimination, weighted
frontier evaluator, exact-branch Hessian, rank-one covariance perturbation
diagnostic, residual margin wrapper, and local-geometry ensemble ideas all
find real structural use.

The notes do not justify a generic non-Gaussian branch or observer selector.
Their strongest exact claims remain covariance/Fisher/quadratic or special
affine-hidden full-law claims. The right 0.3.2 move is therefore not to widen
the module's public law claims. The right move is to build a disciplined
research layer around:

1. general graph-chart frontier Hessians at declared observers;
2. near-exact branch diagnostics as a special residual-controlled case;
3. normalized fibre-volume diagnostics;
4. declared-ladder observer selection versus true Grassmannian optimization;
5. low-rank covariance perturbation classification;
6. local-geometry ensemble statistics and robust verdict transfer.

The biggest hidden win is that 0.3.1 made the project operationally broader
than "Gaussian toy theory" without losing scope discipline. The exact
affine-hidden fibre sector is a genuine full-law sector beyond Gaussian
visible laws. The weighted-frontier and branch-Hessian tools are exact local
quadratic tools. The local-geometry ensemble layer is the right way to start
talking about heterogeneous non-Gaussian systems without pretending that a
single Hessian is a full law.

A later stress pass adds an important correction: exact-branch geometry is not
the parent object for all weighted-frontier critical points. Exact branches
are sufficient for frontier stationarity, but weighted-family first variations
can cancel even when individual family members do not preserve the observer.
The broader order-two parent object is the general graph-chart first and
second variation at a declared observer; the exact-branch Hessian is the
invariant closed sector inside it.

## Quarantine Assessment Of Claude Notes

What to keep:

- Gaussian/covariance bridge identities: visible precision as marginal
  precision, Schur complement identities, tower law, log-det clock as Gaussian
  log-det correlation, and collapse descent.
- The affine-hidden exact sector, especially the state-dependent
  fibre-volume term `0.5 log det D(v)`.
- Staged affine-hidden action shifts as an exact Schur-associative tower.
- Weighted-family frontier as a declared-observer evaluator.
- Exact-branch Hessian as a local fixed-rank diagnostic at exact branches.
- Residual margins as the correct transfer mechanism from quadratic verdicts
  to fuller law claims.
- Local-geometry ensembles as a serious research layer for heterogeneous
  systems.

What to quarantine or rewrite before use:

- "Rank one iff signal has zero visible projection" is false as stated.
- "Observer selection" must mean comparison over a declared ladder unless a
  global search domain and optimizer are supplied.
- Fibre-volume/variational ratios need a declared centered norm and a
  small-denominator guard.
- Closure-adapted neural readouts should not be described as "selecting
  neurons" if the observer is an arbitrary linear subspace.
- KL/Hellinger contraction ratios can be compared locally, but there is no
  invariant divergence-independent fragility score yet.
- Black-hole Page/island language must remain an external-penalty lower
  envelope analogy unless a real gravitational area/topology sector is added.
- Turbulence/intermittency language must remain local-geometry ensemble
  structure, not a derivation of anomalous exponents.

The note verification harnesses were run as demonstrations:

- `python papers/claude_notes/turbulence/verify_all_claims.py`: 11/11 passed.
- `python -X utf8 papers/claude_notes/info-limiting correlations/verify_neural_claims.py`: 13/13 passed.
- The first neural run failed only at console encoding while printing a Unicode
  arrow under Windows `cp1252`; it was not a mathematical failure.

These harnesses mostly verify established kernel identities and toy-model
calculations. They are not enough to promote the speculative readings.

## Probe Findings

### 1. Rank-One Boundary

The exact formula for a covariance perturbation

```text
Sigma_eps = Sigma_0 + eps f f^T
```

under a coordinate visible split is

```text
Delta R =
  gamma w_V w_V^T - beta u_V u_V^T,

u_V = (H_0 f)_V,
w_V = Phi_0 f_V,
beta = eps / (1 + eps f^T H_0 f),
gamma = eps / (1 + eps f_V^T Phi_0 f_V).
```

Therefore `rank(Delta R) <= 2`. One-channel behavior is controlled by
collinearity or vanishing of the two visible directions `u_V` and `w_V`, not by
hidden-only support of `f`.

The probe found:

| Case | Signal visible mass | Background | Rank |
|---|---:|---|---:|
| hidden-only, block diagonal | zero | no visible-hidden covariance coupling | 0 |
| hidden-only, correlated | zero | visible-hidden covariance coupling present | 1 |
| visible-containing, block-aligned | nonzero | aligned block background | 1 |
| generic coloured | nonzero | generic correlated background | 2 |

This is the clean correction to the Claude synthesis. The 0.3.2 theorem target
should be the full two-direction rank classification, not the hidden-only
folk theorem.

### 2. Near-Exact Branch Failure

`exact_branch_hessian` correctly rejects branches once the off-block norm
exceeds its exact-branch tolerance. The probe showed that small off-block
perturbations pass at `1e-10` and `1e-7`, but are rejected at `1e-5` and above
for the test family.

This suggests a high-value 0.3.2 research target: a non-throwing near-branch
diagnostic that reports:

- `off_block_norm`;
- family scale;
- exact-branch tolerance;
- maybe the zeroth-order exact-branch Hessian when a block-diagonal projection
  is explicitly requested;
- a warning that no exact Hessian theorem applies unless an error bound is
  supplied.

The exact Hessian API should remain strict. The new object should be a separate
diagnostic unless the result type is deliberately extended.

### 3. Fibre-Dominance Normalisation

The affine-hidden exact sector is real, but a scalar "fibre dominance ratio"
is not canonical until the norm is declared. In the probe, using a centered
L2 norm over a visible sample grid:

| Coupling amplitude | Variational norm | Fibre norm | Ratio |
|---:|---:|---:|---:|
| 0 | 0 | 0.2734 | undefined |
| 1e-4 | 7.97e-9 | 0.2734 | 3.43e7 |
| 1e-2 | 7.97e-5 | 0.2734 | 3.43e3 |
| 1 | 0.7968 | 0.2734 | 0.343 |

So the "fibre-volume/variational" ratio can be useful, but only with:

- a specified sample measure or domain;
- centered variation, range, or another declared seminorm;
- a denominator floor;
- possibly separate reporting of the numerator and denominator.

The safe research definition is a tuple, not a naked ratio:

```text
FibreDominance_N =
  (||fib - mean fib||_N,
   ||var - mean var||_N,
   ratio if denominator >= delta else undefined)
```

### 4. Gaussian Divergence Contraction Ratios

For small covariance perturbations, KL and Hellinger contraction ratios were
nearly identical in the probe. As perturbations grew, they separated. The
ratio gap moved from about `4e-5` at scale `1e-3` to about `7.6e-2` at scale
`1`.

Implication: a "fragility spectrum" is promising, but it must be indexed by
the chosen divergence. Local agreement follows from shared quadratic/Fisher
onset; global agreement is not an invariant.

### 5. Declared-Ladder Frontier

The weighted frontier evaluator can compare observers in a declared ladder.
The probe used a coordinate ladder and an unrelated rotated ladder. The score
curves and leakage levels differed substantially; both ladders happened to
prefer the largest rank in that toy run, but the absolute frontier values and
intermediate behavior were ladder-dependent.

Implication: the current tool supports "among these observers, with these
weights and this penalty, this one scores best." It does not support "this is
the globally optimal observer" without a Grassmannian optimization theory and
implementation.

## Second-Pass Expansion

After writing the first master plan, I ran a second pass aimed at the plan's
own weak points. Three additions matter.

### A. Near-Branch Hessians Need A Gradient Layer

The second-pass probe shows that near-exact branches generally have first-order
stationarity drift. In the tested family, the graph-score gradient norm scaled
linearly with the off-block norm; the ratio stabilized around `3.95`.

This was the first correction to the near-branch target. A non-throwing
off-block diagnostic is still useful, but a perturbative branch result cannot
only approximate the Hessian. It must report or bound the first variation as
well:

```text
first_variation(U) = gradient of frontier score in Hom(U,U_perp),
second_variation(U) = Hessian-like contract only after stationarity is handled.
```

Research consequence: 0.3.2 should split near-branch analysis into:

- stationarity residual;
- off-block inexactness;
- projected exact-branch Hessian;
- score-margin certification only if the gradient cannot move the optimizer
  across the claimed branch.

Without this split, a "near-branch Hessian" could look stable while the branch
is not even locally critical.

Follow-up research in `audit/0_3_2_theory_research_note.md` now derives both
the graph-chart first variation and the directional second variation for
non-exact branches. The remaining hard problem is no longer the raw Hessian
formula; it is the stability certificate controlling how far the optimizer
moves when the first variation is small but nonzero.

Further follow-up in `audit/0_3_2_frontier_push_note.md` proves an abstract
near-branch stability lemma. It reduces the weighted-frontier specialization
to Hessian-Lipschitz control in the graph chart.

`audit/0_3_2_weighted_frontier_certificate_note.md` now gives a conservative
explicit certificate candidate. The remaining issues are formalization and
sharpness: the projector derivative constants need a line-by-line formal lemma
and are probably too pessimistic for a public diagnostic without refinement.

`audit/0_3_2_projector_bounds_note.md` supplies a first line-by-line
rho-dependent refinement. Brutality checks found no projector-bound violations
and showed roughly a factor-two improvement in the certificate constant on
synthetic near-branch families, but the certificate remains stringent.

### A2. Exact Branch Is Not The Parent Of All Stationary Frontier Geometry

The non-exact stationary frontier pass in
`audit/0_3_2_non_exact_stationary_frontier_note.md` isolates a clean
counterexample to a tempting extension:

```text
exact-branch Hessian + off-block tolerance
```

is not the right general frontier Hessian.

For a one-dimensional observer in `R^2` and a family member

```text
A_i = [[a_i,e_i],
       [e_i,d_i]],
```

the weighted-frontier first variation at the coordinate observer is

```text
D F_mu(0) = 2 sum_i w_i e_i ((2+mu)a_i - mu d_i).
```

Thus exact branch, `e_i=0` for all `i`, implies stationarity but is not
necessary. Off-blocks can cancel across the family or through the leakage
penalty.

The sharp breakpoint example is

```text
A_+ = [[3,e],[ e,1]],
A_- = [[3,-e],[-e,1]],
mu = 0.
```

The first variation vanishes for every `e`, but

```text
D^2 F_0(0) = -24 + 8e^2.
```

The invalid exact-branch proxy stays at `-24`, while the true second
variation changes sign at `|e|=sqrt(3)`. Numerics in
`audit/outputs/0_3_2_non_exact_stationary_frontier_check.json` reproduce this
threshold and find random stationary non-exact local maxima, local minima, and
near-degenerate cases.

Research consequence: the next frontier object should be a general
declared-observer graph Hessian plus stationarity residual. Near-exact branch
certificates remain useful, but only as one residual-controlled route, not as
the full theory of frontier criticality.

Follow-up in `audit/0_3_2_general_graph_hessian_note.md` makes that parent
object explicit by polarizing the exact directional second-variation formula.
It reduces to the module's `exact_branch_hessian` in the invariant sector with
max error `3.55e-15` in a random check, while still matching finite-difference
Hessians off the exact-branch sector.

An additional invariance pass in
`audit/0_3_2_graph_hessian_invariance_check.py` removes the remaining
coordinate ambiguity. In arbitrary embedded observer frames, the first
variation matched finite differences at `9.52e-09`, the bilinear Hessian
matched finite differences with relative error `4.62e-06`, frame-change
invariance held to `5.68e-13`, and exact-branch reduction matched the module
to `1.42e-14` once expressed in the same complement frame.

Follow-up in `audit/0_3_2_declared_frontier_certificate_note.md` then widens
the certificate interpretation: the same abstract Lipschitz lemma applies to
any declared observer with a small stationarity residual and definite general
graph Hessian. This is broader than near-exact branches. The stress check
certifies stationary non-exact maxima for `e=1` and `e=1.7` in the symmetric
pair, fails at the degenerate threshold `e=sqrt(3)`, and certifies local
minima after the sign flip. Imbalanced cases show the expected conservatism:
small residuals pass, larger residuals become vacuous under `L_cert`.

### A3. Shelved Feature Recovery

`audit/0_3_2_shelved_feature_recovery_note.md` reviews older withheld ideas
under the new scope discipline. The important recovery rule is:

```text
declared finite object + exact local/special-sector theorem + explicit scope
```

The pass found several low-hanging or near-low-hanging candidates:

- rank-k covariance/Fisher perturbation classification is theorem-grade;
- declared-ladder dimension-cost intervals are theorem-grade for supplied
  ladders;
- affine-hidden stage signs and branch reversal are exact inside the supplied
  affine-hidden Gaussian-fibre sector;
- residual-margin and ensemble finite-candidate wrappers are safe when honest
  residual/sampling bounds are supplied;
- guarded fibre dominance is useful as a diagnostic tuple, not a naked ratio;
- declared local frontier certificates are now the right future surface, but
  remain gated by conservative `L_cert` semantics.

Still blocked: generic non-Gaussian branch probabilities, global
noncommuting Grassmannian optimization, and probability-support event engines
from Hessians alone.

### B. Variable-Rank Frontier Selection Needs A Dimension Budget

The weighted frontier score is exact for a supplied observer. It is not a
complete variable-rank selector. If full rank is allowed and no dimension cost
is charged, full rank is structurally privileged:

```text
P = I gives leakage = 0 and visible score = Tr(M_nu).
```

For any lower-rank `P`,

```text
penalized_score(P) = S(P) - mu L(P) <= S(P) + L(P) = Tr(P M_nu) <= Tr(M_nu).
```

So "choose the rank m maximizing penalized score" is ill-posed unless:

- rank is fixed;
- a maximum rank budget is fixed;
- a dimension cost is included;
- or the comparison is explicitly over a declared finite ladder with a
  predeclared stopping rule.

The second-pass probe found full rank best with zero cost; only a strong enough
dimension cost moved the best rank below full. This should be treated as a
theorem-grade correction to any future observer-selection language.

### C. Rank-k Covariance Perturbations Are A Natural Extension

The rank-one/two-rank formula appears to extend cleanly. For

```text
Sigma_1 = Sigma_0 + F F^T,  F in R^{n x k},
```

Woodbury gives:

```text
Delta R =
  W (I + F_V^T Phi_0 F_V)^(-1) W^T
  - U (I + F^T H_0 F)^(-1) U^T,

U = (H_0 F)_V,
W = Phi_0 F_V.
```

Therefore `rank(Delta R) <= 2k`. The second-pass probe with `k=2` found
generic rank `4`, formula residual `1.11e-16`, hidden-only correlated rank
`2`, and hidden-only block-diagonal rank `0`.

Research consequence: the 0.3.2 perturbation theorem should probably be
rank-k first, with rank-one as the corollary. That gives a stronger and cleaner
result for low-rank Fisher/covariance updates.

## Main 0.3.2 Research Targets

### Target 1: Low-Rank Covariance Perturbation Theorem

Prove and write the exact classification:

- for a rank-k covariance update, the hidden-gap increment has rank at most
  `2k`;
- rank zero occurs by vanishing/cancellation;
- in the rank-one case, one-channel behavior occurs exactly when the two
  weighted rank-one visible directions span a one-dimensional subspace or one
  term vanishes;
- generic coloured rank-one backgrounds give rank two;
- white/aligned rank-one backgrounds give the one-channel sector.

Deliverables:

- theorem statement;
- degenerate cases;
- numerical examples from `audit/0_3_2_research_probe.py`;
- API doc refinement if needed.

Do not state hidden-only support as necessary or sufficient for nonzero
one-channel behavior.

Second-pass upgrade: formulate this as the `rank <= 2k` low-rank covariance
perturbation theorem, then specialize to `k=1`. The rank-one theorem should
not be isolated if the rank-k proof is just Woodbury plus visible-block
Woodbury.

### Target 2: Near-Exact Branch Diagnostic

Research a strict diagnostic layer adjacent to `exact_branch_hessian`:

```text
branch_inexactness(family, B) =
  max_i ||B_perp^T A_i B||_F / family_scale.
```

Questions:

- What is the natural scale: max Frobenius norm, weighted RMS norm, or moment
  operator norm?
- Does a projected block-diagonal family give a meaningful zeroth-order
  Hessian?
- Can the difference between true graph objective and projected exact-branch
  Hessian be bounded by `O(off_block_norm)` or `O(off_block_norm^2)` under
  spectral-gap assumptions?
- Which result should be returned when the branch is nearly exact but not
  exact?

This is likely the lowest-risk 0.3.2 module candidate after the theorem is
written, because it improves failure diagnosis without widening claims.

Second-pass correction: include a stationarity residual. Off-block norm alone
is not a verdict. The graph-score gradient can be first order in the off-block
term, so the diagnostic should probably be named "branch inexactness and
stationarity", not just "near-branch Hessian".

### Target 3: Perturbative Branch Hessian

Go beyond diagnostics only if the error bound is proved. Candidate theorem:

Let `A_i` be decomposed into block-diagonal part `A_i^0` plus off-block part
`E_i` relative to `U`. The graph objective Hessian at `U` equals the exact
branch Hessian for `A_i^0` plus a controlled correction whose norm is bounded
by constants depending on family scale, weights, `mu`, and `max_i ||E_i||`.

Research questions:

- Is the first correction linear or quadratic in off-block norm?
- Does stationarity fail at first order even when the Hessian-like contract is
  close?
- Should the correction be reported as a bound, not as an approximate Hessian?
- How does degeneracy of the exact-branch Hessian affect certification?

Do not implement an approximate stability classifier without this theorem.

Second-pass correction: any perturbative Hessian theorem must be conditional
on a first-variation bound or must explicitly report that the supplied plane is
not stationary. A nonzero gradient can dominate Hessian status for branch
selection.

### Target 4: Fibre-Dominance Functional

Define a robust diagnostic for variable-precision affine-hidden fibres:

```text
var(v) = A(v) - 0.5 J(v)^T D(v)^(-1) J(v)
fib(v) = 0.5 log det D(v)
vis(v) = var(v) + fib(v)
```

Candidate outputs:

- centered norm of `fib`;
- centered norm of `var`;
- centered norm of `vis`;
- ratio with denominator guard;
- sign and magnitude of action-shift stages;
- branch-verdict reversal flag when `argmin(var)` differs from `argmin(vis)`.

Research hazards:

- sample-grid dependence;
- additive constants;
- denominators near zero;
- discrete versus continuous visible domains;
- whether ratio or dominance should be measured by range, L2, Linf, or
  user-supplied weights.

### Target 5: Staged Elimination Screening And Anti-Screening

Claude's "anti-screening" reading of negative action shifts is plausible in
the affine-hidden sector. The exact stage shift is

```text
Delta A_E = 0.5 log det D_EE - 0.5 J_E^T D_EE^(-1) J_E.
```

This directly gives a sign criterion:

```text
screening if log det D_EE > J_E^T D_EE^(-1) J_E,
anti-screening if log det D_EE < J_E^T D_EE^(-1) J_E.
```

Research tasks:

- check sign convention against the visible action convention;
- handle non-singleton stage blocks;
- study order-dependence of per-stage signs even when total tower law is exact;
- decide whether classification belongs in docs, audit, or API metadata.

This is an easy theory win if stated strictly inside the affine-hidden sector.

### Target 6: Declared-Ladder Observer Frontiers

Build a research doctrine for observer selection:

- `weighted_family_frontier_scores`: exact evaluator for supplied observers.
- declared ladder selection: exact comparison within an explicitly supplied
  finite list or nested flag.
- exact branch Hessian: local fixed-rank status at exact branches.
- noncommuting global optimization: open Grassmannian problem.

Research questions:

- How should penalties `mu` be calibrated?
- Can critical `mu` values be solved exactly for diagonal/exact-branch
  families?
- Can robust selection be certified by score gaps plus residual/inexactness
  bounds?
- What metadata should record the declared observer class?

Do not imply a global selector unless the search domain is part of the input
and the optimizer has convergence/failure semantics.

Second-pass correction: variable-rank selection without dimension cost is
mathematically degenerate toward full rank. A declared ladder is not enough if
the ladder includes full rank and the question is "which rank is best?" The
problem must include a rank budget, dimension penalty, target compression
level, or explicit stopping rule.

### Target 7: Local-Geometry Ensembles

The ensemble layer is a major bridge between non-Gaussian heterogeneity and
exact quadratic geometry:

```text
H_i sampled local SPD Hessians/Fisher objects
Phi_i = Phi_C(H_i)
optional ceiling T_i gives Lambda_i, tau_i
```

0.3.2 should research:

- tail-sensitive summaries of `tau_i`, ranks, and log-det visible summaries;
- Jensen/delta-method warnings for replacing `E f(H)` by `f(E H)`;
- robust observer comparison under ensemble spread;
- bootstrap or concentration certificates for finite samples;
- whether branch events are better detected by distributions of local scores
  than by a single average Hessian.

This is not a cumulant engine. It is exact samplewise local geometry plus
statistics over samples.

### Target 8: Residual-Margin Doctrine For All Verdicts

Generalize the existing residual-margin wrapper from scalar branch gaps to:

- observer scores;
- declared-ladder frontier comparisons;
- local Hessian stability;
- approximate branch-Hessian status;
- ensemble mean score comparisons.

The core theorem remains elementary:

```text
gap > residual_bound_for_candidate_A + residual_bound_for_candidate_B
```

The hard part is not the inequality. The hard part is obtaining honest
residual bounds from cumulants, sampling error, inexact branches, or model
misspecification.

### Target 9: Domain-Specific Research Lines

Black-hole/information analogy:

- treat Page/island-like selection as an external-penalty lower envelope
  `min_alpha(P_alpha + tau_alpha)`;
- keep area terms and topology outside the module;
- study branch-switch kinks as observer-prescription events only in a lifted
  prescription space.

Turbulence:

- local-geometry ensembles over windows or scales;
- clock and hidden-load tail distributions;
- stage shifts as cascade-like screening/anti-screening diagnostics;
- declared filter ladders rather than global LES-filter claims;
- no claim of anomalous exponent derivation without law-level cumulants.

Neural coding:

- rank-one/two-rank differential-correlation classification;
- critical recorded dimension as projection of signal direction into a
  declared subpopulation/observer ladder;
- strict distinction between arbitrary linear observer dimensions and actual
  neuron subset selection;
- collapse depth/fragility spectrum by divergence and observer ladder.

## Low-Hanging Wins

1. Write the corrected rank-k covariance perturbation lemma into the research
   stack, with rank one as a corollary.
2. Add a short theory note on screening versus anti-screening in
   affine-hidden staged elimination.
3. Define a guarded fibre-dominance diagnostic in research language, not as an
   API field yet.
4. Specify a near-branch diagnostic theorem that includes stationarity
   residual, not just off-block norm.
5. Add declared-ladder plus dimension-budget language to future
   observer-selection discussions.
6. Build a small empirical protocol for local-geometry ensembles:
   input local Hessian/Fisher samples, output exact samplewise geometry plus
   uncertainty summaries.
7. Write residual-margin rules for ensemble and frontier comparisons.

## Weak Points To Attack Before Implementation

- Near-exact branch theory may fail by stationarity drift even if Hessian
  shape is close. Need a theorem, not just an off-block number.
- Fibre-dominance can become meaningless when the variational term is flat.
  Guard denominators and report raw norms.
- Frontier maxima over rank can be trivial if the score increases with
  dimension. Penalty or dimension budget must be part of the problem.
- Divergence contraction curves are divergence-dependent away from the local
  Fisher regime.
- Local-geometry ensembles do not recover mixture mass, skew, tail, or support
  events unless the sampling scheme actually captures them.
- Closure-adapted linear subspaces are not necessarily sparse or biological
  coordinate subsets.
- Black-hole, turbulence, and neural interpretations can be valuable examples
  only if every statement is labelled covariance-level, local-quadratic,
  exact affine-hidden, or speculative analogy.

## Suggested Theory-To-Validation Pipeline For 0.3.2

1. Prove the rank-k low-rank covariance perturbation classification including
   degeneracies.
2. Prove the affine-hidden stage-shift sign criterion.
3. Define fibre-dominance as a guarded seminorm-based diagnostic.
4. Define branch inexactness plus stationarity residual and prove what they do
   and do not certify.
5. Attempt a perturbative branch Hessian bound; stop if the bound is not clean.
6. Define declared-ladder observer selection with rank budgets or dimension
   costs, then add residual-margin comparison.
7. Extend residual-margin statements to ensemble score estimates.
8. Stress test all statements on block-diagonal, correlated, nearly singular,
   degenerate-eigenvalue, and high-condition examples.
9. Only then decide which, if any, 0.3.2 items become module APIs.

## Possible 0.3.2 Module Candidates After Theory Approval

These are not approved by this document. They are candidates only.

- `branch_inexactness(...)`: non-throwing off-block diagnostic.
- `branch_stationarity(...)`: graph-score first-variation diagnostic, if a
  stable chart convention is chosen.
- `projected_branch_hessian(...)`: only if its approximation and warning
  semantics are theorem-backed.
- `fibre_dominance(...)`: only as a guarded diagnostic with declared weights
  or sample measure.
- `affine_stage_shift_classification(...)`: likely safe if it simply reports
  the exact stage-shift sign in the existing affine-hidden sector.
- `frontier_ladder_scores(...)`: a convenience wrapper over declared observer
  lists, with explicit rank budget or dimension-cost semantics, not a
  Grassmannian optimizer.
- `ensemble_score_margin(...)`: only with explicit supplied sampling/residual
  bounds.

Do not implement:

- generic non-Gaussian full-law branch selector;
- global noncommuting observer optimizer;
- probability-support event engine from Hessians alone;
- black-hole Page curve engine;
- turbulence intermittency exponent engine;
- neural biological subpopulation selector from arbitrary linear observers.

## Records And Validation State

Created:

- `audit/0_3_2_research_probe.py`
- `audit/outputs/0_3_2_research_probe.json`
- `audit/0_3_2_second_pass_probe.py`
- `audit/outputs/0_3_2_second_pass_probe.json`
- `audit/0_3_2_near_branch_variation_check.py`
- `audit/outputs/0_3_2_near_branch_variation_check.json`
- `audit/0_3_2_master_research_doc.md`
- `audit/0_3_2_theory_research_note.md`
- `audit/0_3_2_frontier_push_note.md`
- `audit/0_3_2_weighted_frontier_certificate_note.md`
- `audit/0_3_2_projector_bounds_note.md`
- `audit/0_3_2_non_exact_stationary_frontier_note.md`
- `audit/0_3_2_non_exact_stationary_frontier_check.py`
- `audit/outputs/0_3_2_non_exact_stationary_frontier_check.json`
- `audit/0_3_2_general_graph_hessian_note.md`
- `audit/0_3_2_general_graph_hessian_check.py`
- `audit/outputs/0_3_2_general_graph_hessian_check.json`
- `papers/0.3.2_Technical_Note_1.tex`
- `audit/0_3_2_implementation_plan.md`

Read and digested:

- `papers/claude_notes/v031_reapproach_synthesis.md`
- `papers/claude_notes/black hole analysis/black_hole_information_analysis.md`
- `papers/claude_notes/turbulence/turbulence_intermittency_analysis.md`
- `papers/claude_notes/info-limiting correlations/neural_coding_information_analysis.md`
- selected note verification scripts.

Validation commands run:

- `python audit/0_3_2_research_probe.py`
- `python audit/0_3_2_second_pass_probe.py`
- `python audit/0_3_2_near_branch_variation_check.py`
- `python audit/0_3_2_certificate_brutality_check.py`
- `python audit/0_3_2_non_exact_stationary_frontier_check.py`
- `python audit/0_3_2_general_graph_hessian_check.py`
- `python audit/0_3_2_graph_hessian_invariance_check.py`
- `python audit/0_3_2_declared_frontier_certificate_check.py`
- `python audit/0_3_2_shelved_feature_recovery_check.py`
- `python papers/claude_notes/turbulence/verify_all_claims.py`
- `python -X utf8 "papers/claude_notes/info-limiting correlations/verify_neural_claims.py"`
- `python -m pytest tests/test_perturbation.py tests/test_frontier.py tests/test_affine.py -q`
- `python -m pytest tests/test_frontier.py tests/test_perturbation.py -q`
- `python -m pytest -q`
- from `nomodescent/`: `python -m pytest -q`
- from `evidence/`: `python -m pytest -q`
- TeX compile check:
  `pdflatex -interaction=nonstopmode -halt-on-error -output-directory papers papers/0.3.2_Technical_Note_1.tex`

The TeX compile produced no document errors after cleanup. The generated
auxiliary/PDF build artifacts were removed, leaving the paper source as the
tracked research record.

No module source files were intentionally modified for this research pass.

## Final Research Position

The project is past Gaussian toy territory, but not because Gaussian claims can
be casually exported. It is past that point because it now has a layered stack:

- exact SPD/Fisher/local-quadratic observation geometry;
- exact Gaussian law mode;
- exact special full-law sectors such as affine-hidden Gaussian fibres;
- exact residual-margin transfer when residuals are supplied;
- an emerging ensemble and higher-jet layer for genuinely non-Gaussian
  structure.

For 0.3.2, the highest-value work is to make the transition layer precise:
general graph-chart frontier Hessians, declared local frontier certificates,
near-exact branches as a residual case, guarded fibre-volume diagnostics,
declared-ladder selection, and robust transfer from quadratic verdicts to
richer laws. That is where the new research value is. Gaussian is not the true
boundary of the programme, but full non-Gaussian verdicts still require
explicit higher-order, residual, support, or ensemble data.
