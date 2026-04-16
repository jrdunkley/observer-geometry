# 0.3.2 Theory Research Note

Date: 2026-04-11

Status: pure research note. This file proves or delimits the 0.3.2
opportunities identified in `audit/0_3_2_master_research_doc.md`. It does not
authorize implementation. Module source was not modified.

Numerical companions:

- `audit/0_3_2_research_probe.py`
- `audit/0_3_2_second_pass_probe.py`
- `audit/0_3_2_near_branch_variation_check.py`
- `audit/outputs/0_3_2_research_probe.json`
- `audit/outputs/0_3_2_second_pass_probe.json`
- `audit/outputs/0_3_2_near_branch_variation_check.json`

Further frontier push:

- `audit/0_3_2_frontier_push_note.md`
- `audit/0_3_2_weighted_frontier_certificate_note.md`
- `audit/0_3_2_projector_bounds_note.md`
- `audit/0_3_2_non_exact_stationary_frontier_note.md`
- `audit/0_3_2_general_graph_hessian_note.md`
- `audit/0_3_2_graph_hessian_invariance_check.py`
- `audit/0_3_2_declared_frontier_certificate_note.md`
- `audit/0_3_2_shelved_feature_recovery_note.md`

## Executive Research State

There are several theorem-grade wins ready to write formally:

1. a low-rank covariance perturbation theorem: rank-k covariance updates move
   the visible hidden-gap by a difference of two rank-k positive terms, hence
   rank at most `2k`;
2. a full-rank degeneracy theorem for variable-rank weighted-frontier
   selection: without a rank budget or dimension cost, the full observer
   dominates the penalized frontier for `mu >= 0`;
3. an affine-hidden staged-elimination sign identity: in a fixed measure and
   coordinate convention, the stage action shift is positive or negative
   exactly according to a log-det versus saddle-term inequality.
4. residual-margin and ensemble finite-candidate comparison theorems, once
   residual or sampling bounds are supplied.
5. rank-k covariance/Fisher perturbation classification, not only rank-one.
6. declared-ladder dimension-cost interval calculus for supplied finite
   ladders.

There are now three useful diagnostics or theorem targets that are not yet
module-complete enough for implementation:

1. near-branch analysis must include stationarity drift, not just off-block
   norm;
2. general graph-chart frontier Hessians must be treated as the parent object
   for stationary, possibly non-exact, observers;
3. declared local frontier certificates are possible, but currently rely on
   conservative Hessian-Lipschitz bounds and can be vacuous;
4. fibre-dominance needs an explicit seminorm, sample measure, and denominator
   guard.

The cleanest current 0.3.2 research package is therefore:

- theorem note for low-rank perturbations;
- theorem note for declared-ladder and dimension-budget observer selection;
- theorem note for affine-hidden stage signs with coordinate caveats;
- theorem note for residual and ensemble finite-candidate margins when bounds
  are supplied;
- theorem note for rank-k covariance/Fisher perturbations;
- theorem note for declared-ladder dimension-cost phase diagrams;
- theorem note for general graph-chart frontier first/second variation, with
  exact-branch Hessian recovered as the invariant special case;
- theorem note for declared local frontier certificates using stationarity
  residual, Hessian margin, and chart Lipschitz control;
- research-only diagnostic definitions for fibre dominance and near-branch
  first/second variation.

## 1. Low-Rank Covariance Perturbation Theorem

### Theorem

Let `Sigma_0` be SPD on `R^n` and let the first `m` coordinates be visible.
Write

```text
Sigma_1 = Sigma_0 + F F^T,
F in R^{n x k}.
```

This includes any PSD rank-k update `B R B^T` with `R >= 0` by taking
`F = B R^{1/2}` on the active support of `R`.

Let

```text
H_0 = Sigma_0^{-1},
H_1 = Sigma_1^{-1},
Phi_0 = (Sigma_0)_{VV}^{-1},
Phi_1 = (Sigma_1)_{VV}^{-1}.
```

Define the visible hidden-gap object

```text
R(Sigma) = H_{VV} - Phi,
```

where `Phi = Sigma_{VV}^{-1}`. Also define

```text
U = (H_0 F)_V,
W = Phi_0 F_V,
K = (I_k + F^T H_0 F)^{-1},
G = (I_k + F_V^T Phi_0 F_V)^{-1}.
```

Then

```text
R(Sigma_1) - R(Sigma_0) = W G W^T - U K U^T.
```

Consequently,

```text
rank(R(Sigma_1) - R(Sigma_0)) <= rank(W) + rank(U) <= 2k.
```

### Proof

Woodbury gives

```text
H_1 = H_0 - H_0 F (I_k + F^T H_0 F)^{-1} F^T H_0.
```

Taking the visible block,

```text
(H_1)_{VV} - (H_0)_{VV} = - U K U^T.
```

The visible covariance block is

```text
(Sigma_1)_{VV} = (Sigma_0)_{VV} + F_V F_V^T.
```

Applying Woodbury again to the visible covariance block,

```text
Phi_1 = Phi_0 - Phi_0 F_V (I_k + F_V^T Phi_0 F_V)^{-1} F_V^T Phi_0
      = Phi_0 - W G W^T.
```

Therefore

```text
R(Sigma_1) - R(Sigma_0)
= [(H_1)_{VV} - Phi_1] - [(H_0)_{VV} - Phi_0]
= [-U K U^T] - [-W G W^T]
= W G W^T - U K U^T.
```

The rank bound follows from `rank(A-B) <= rank(A)+rank(B)` and
`rank(WGW^T) <= rank(W) <= k`, `rank(UKU^T) <= rank(U) <= k`.

### Rank-One Corollary

For `k=1`, write `F=f`. The update is

```text
Delta R =
  gamma w_V w_V^T - beta u_V u_V^T,

u_V = (H_0 f)_V,
w_V = Phi_0 f_V,
beta = (1 + f^T H_0 f)^{-1},
gamma = (1 + f_V^T Phi_0 f_V)^{-1}.
```

With an explicit strength `epsilon`, replace `F` by `sqrt(epsilon) f`, giving
the same formula with the usual `epsilon` factors.

The rank-one update has rank:

- `0` if both terms vanish or cancel exactly;
- `1` if the two visible directions span a one-dimensional subspace or one
  term vanishes;
- generically `2` in a coloured background when `u_V` and `w_V` are not
  collinear.

The hidden-only support condition `f_V=0` is neither a full theorem nor a
nonzero one-channel criterion:

- if the background is block diagonal, `f_V=0` gives `W=0` and `U=0`, hence
  rank `0`;
- if the background has visible-hidden covariance coupling, `f_V=0` gives
  `W=0` but can give `U != 0`, hence rank `1`;
- if the signal has visible support but `u_V` and `w_V` are aligned, the rank
  can still be `1`.

### Numerical Status

`audit/0_3_2_research_probe.py` found all three rank-one cases:

- hidden-only, block diagonal: rank `0`;
- hidden-only, correlated: rank `1`;
- visible-containing, block aligned: rank `1`;
- generic coloured: rank `2`.

`audit/0_3_2_second_pass_probe.py` tested `k=2` and found:

- generic rank `4`;
- formula residual `1.11e-16`;
- hidden-only correlated rank `2`;
- hidden-only block-diagonal rank `0`.

### Research Consequence

The formal result should be rank-k first. The earlier rank-one story is a
special case. This is a clean theorem-grade result and likely one of the best
0.3.2 paper wins.

## 2. Declared-Ladder And Dimension-Budget Observer Selection

### Setup

Let `A_i` be symmetric operators, weights `w_i >= 0`, and

```text
M = sum_i w_i A_i^2.
```

For an orthogonal projector `P`, define

```text
S(P) = sum_i w_i ||P A_i P||_F^2,
L(P) = sum_i w_i 1/2 ||[A_i, P]||_F^2.
```

The finite weighted-family identity is

```text
Tr(P M) = S(P) + L(P).
```

The implemented frontier score is

```text
F_mu(P) = S(P) - mu L(P),   mu >= 0.
```

### Theorem: Full-Rank Degeneracy

If the admissible observer class includes the full-rank projector `I`, then
for every projector `P` and every `mu >= 0`,

```text
F_mu(P) <= F_mu(I).
```

### Proof

Since `L(P) >= 0` and `mu >= 0`,

```text
F_mu(P) = S(P) - mu L(P) <= S(P) + L(P) = Tr(P M).
```

Since `M` is PSD and `0 <= P <= I`,

```text
Tr(P M) <= Tr(M).
```

For `P=I`,

```text
L(I)=0,
S(I)=sum_i w_i ||A_i||_F^2 = Tr(M),
F_mu(I)=Tr(M).
```

Therefore `F_mu(P) <= F_mu(I)`.

### Consequence

Variable-rank observer selection is ill-posed unless the problem includes one
of the following:

- a fixed rank;
- a maximum rank budget excluding full rank;
- a dimension cost;
- a compression target;
- an explicit stopping rule on a declared ladder;
- an external application constraint.

"Choose the rank maximizing penalized score" is not a meaningful generic rule
when full rank is admissible and `mu >= 0`.

### Dimension-Cost Variant

A safe variable-rank objective is of the form

```text
F_{mu,c}(P) = S(P) - mu L(P) - c rank(P),
```

or a constrained form

```text
maximize F_mu(P) subject to rank(P) <= r_max.
```

The dimension cost `c` is not determined by the quadratic geometry alone. It
is a modelling or design parameter. This should be explicit in any future
observer-selection doctrine.

### Numerical Status

`audit/0_3_2_second_pass_probe.py` used a strong declared ladder built from
the top eigenspaces of `M`. With zero dimension cost, full rank won. Only after
adding a large enough dimension cost did a lower rank win.

### Research Consequence

0.3.2 should use the language:

```text
declared-ladder comparison under a declared rank budget or dimension cost.
```

It should not say:

```text
the frontier score selects the observer.
```

unless the admissible class and rank/cost policy are part of the theorem.

## 3. Affine-Hidden Staged Elimination: Screening Sign

### Setup

In one affine-hidden stage, split hidden coordinates into kept `K` and
eliminated `E`. The hidden fibre has action

```text
A + J_K^T h_K + J_E^T h_E
  + 1/2 h_K^T D_KK h_K
  + h_K^T D_KE h_E
  + 1/2 h_E^T D_EE h_E.
```

Assume `D_EE` is SPD. Completing the square in `h_E` gives the new action
shift

```text
Delta A_E = 1/2 log det D_EE - 1/2 J_E^T D_EE^{-1} J_E,
```

with shifted kept coupling and precision

```text
J_K' = J_K - D_KE D_EE^{-1} J_E,
D_KK' = D_KK - D_KE D_EE^{-1} D_EK.
```

### Sign Criterion

Within a fixed hidden measure and coordinate convention:

```text
Delta A_E > 0  iff  log det D_EE > J_E^T D_EE^{-1} J_E,
Delta A_E = 0  iff  log det D_EE = J_E^T D_EE^{-1} J_E,
Delta A_E < 0  iff  log det D_EE < J_E^T D_EE^{-1} J_E.
```

It is reasonable, as research language, to call the positive case
"screening" and the negative case "anti-screening" if the convention is stated.

### Caveat: Absolute Sign Needs A Measure Convention

The absolute log-determinant term depends on the hidden coordinate/measure
normalization. Since the visible action is only defined up to an additive
constant independent of visible variables, an isolated constant stage shift can
change under a constant hidden-coordinate rescaling.

Therefore the most invariant claims are:

- sign of variation of `Delta A_E(v)` over visible states after fixing a
  measure convention;
- branch differences between candidates evaluated in the same convention;
- total visible-action differences, not arbitrary absolute constants.

For a single scalar stage reported as a module diagnostic, the sign is useful
only if the caller understands the coordinate and measure normalization.

### Research Consequence

The screening/anti-screening result is mathematically simple, but the public
language must avoid implying a coordinate-free physical sign unless the
measure convention is part of the setup.

## 4. Guarded Fibre-Dominance Diagnostic

### Exact Sector

For the variable-precision affine-hidden fibre

```text
p(v,h) proportional to
exp(-A(v) - 1/2 h^T D(v) h - J(v)^T h),
```

the exact visible action, up to a global additive constant, is

```text
S_vis(v) =
A(v) - 1/2 J(v)^T D(v)^{-1} J(v) + 1/2 log det D(v).
```

Define

```text
Var(v) = A(v) - 1/2 J(v)^T D(v)^{-1} J(v),
Fib(v) = 1/2 log det D(v),
Vis(v) = Var(v) + Fib(v).
```

### Problem With A Naked Ratio

The ratio

```text
||Fib|| / ||Var||
```

is not canonical:

- actions are defined up to additive constants;
- the norm depends on a sample grid or measure;
- the denominator can be zero or tiny;
- range, L2, Linf, and weighted norms answer different questions.

`audit/0_3_2_research_probe.py` showed exactly this instability: when the
variational part was nearly flat, the fibre/variational ratio jumped from
undefined to `3.4e7`, while the actual fibre norm stayed fixed.

### Research Definition

Let `(Omega, nu)` be a declared finite sample measure or visible-domain
measure, and let `N_nu` be a seminorm insensitive to constants, for example:

```text
N_nu(f) = ||f - E_nu f||_{L2(nu)}
```

or a centered range:

```text
N_range(f) = max f - min f.
```

For a denominator floor `delta > 0`, define

```text
FibreDominance = {
  fibre_norm = N_nu(Fib),
  variational_norm = N_nu(Var),
  visible_norm = N_nu(Vis),
  ratio = fibre_norm / variational_norm
          if variational_norm >= delta
          else undefined,
  denominator_status = "regular" or "flat-variational"
}.
```

### Branch-Reversal Diagnostic

A more verdict-relevant diagnostic is:

```text
argmin Var != argmin Vis
```

on a declared discrete candidate set or under a declared continuous-domain
optimization theorem. This directly tests whether fibre volume changes the
branch verdict. It is often more meaningful than the ratio itself.

### Research Consequence

Fibre dominance is a good 0.3.2 concept, but it should enter as a guarded
research diagnostic:

- always report raw centered norms;
- make the ratio optional/undefined under a denominator floor;
- define the sample measure;
- separately report branch-verdict reversal.

## 5. Near-Branch Stationarity And Hessian Split

### Setup

Use a basis adapted to a candidate observer subspace `U` and complement `W`.
For each symmetric family member,

```text
A_i = [ a_i   e_i^T
        e_i   d_i   ],
```

where `e_i = W^T A_i U` is the off-block. The exact-branch condition is
`e_i=0` for all `i`.

Let `X in Hom(U,W)` be the graph-chart coordinate. At `X=0`, the tangent
projector variation is

```text
dot P_X = [0  X^T
           X  0  ].
```

The weighted frontier objective is

```text
F_mu(P) = S(P) - mu L(P)
        = (1+mu) S(P) - mu Tr(P M),
M = sum_i w_i A_i^2.
```

### First Variation Formula

For one operator,

```text
d S_i[X] = 4 <e_i a_i, X>_F.
```

Also,

```text
(A_i^2)_{WU} = e_i a_i + d_i e_i,
d Tr(P A_i^2)[X] = 2 <e_i a_i + d_i e_i, X>_F.
```

Therefore

```text
d F_mu[X]
= 2 sum_i w_i < (2+mu) e_i a_i - mu d_i e_i, X >_F.
```

The graph-chart gradient is therefore

```text
Grad F_mu(U)
= 2 sum_i w_i [ (2+mu) e_i a_i - mu d_i e_i ].
```

This vanishes when `e_i=0` for all `i`, recovering exact-branch stationarity.
But for a near branch it is generically first order in the off-blocks.

### Consequence

A near-branch diagnostic cannot be only:

```text
max_i ||e_i||_F.
```

It should include:

```text
stationarity_norm = ||Grad F_mu(U)||_F,
off_block_norm = max_i ||e_i||_F,
family_scale,
projected_exact_branch_hessian_if_requested,
```

plus explicit wording that the exact-branch Hessian theorem does not classify
local stability until the first variation is zero or controlled.

### Numerical Status

`audit/0_3_2_second_pass_probe.py` showed stationarity drift proportional to
off-block norm, with gradient/off-block ratio stabilizing around `3.95` in the
test family. That confirms the first-variation warning.

A separate finite-difference check in
`audit/0_3_2_near_branch_variation_check.py` gave:

```text
max_abs_error 8.24e-09
relative_error 3.05e-10
grad_norm 47.0596, fd_norm 47.0596
```

So the first-variation formula is numerically consistent with the graph-chart
objective to finite-difference precision.

### Directional Second Variation Formula

The second variation can also be written explicitly. For a single operator
`A_i`, define

```text
B_i(X) = e_i^T X + X^T e_i,
M_i = A_i^2,
M_i,UU = (M_i)_{UU},
M_i,WW = (M_i)_{WW}.
```

For the graph path `t -> graph(tX)`, the coefficient of `t^2` in `S_i` is

```text
S_{i,2}(X)
= Tr(B_i(X)^2)
  + 2 Tr(a_i X^T d_i X)
  - 2 Tr(X^T X a_i^2).
```

The coefficient of `t^2` in `Tr(P(t) M_i)` is

```text
T_{i,2}(X)
= Tr(X^T M_i,WW X) - Tr(X^T X M_i,UU).
```

Therefore the directional second derivative of the weighted frontier is

```text
D^2 F_mu(U)[X,X]
= 2 sum_i w_i [
     (1+mu) S_{i,2}(X)
     - mu T_{i,2}(X)
   ].
```

This formula is valid whether or not `U` is an exact branch. It is not by
itself a local optimality certificate when the first variation is nonzero.

Finite-difference check for this formula:

```text
h=3e-4:
finite difference = -22.5263351676
formula           = -22.5263367179
absolute error    = 1.55e-06
relative error    = 6.88e-08
```

The larger errors at smaller `h` were the expected floating-point subtraction
limit.

### Abstract Stability Certificate

The useful theorem is not just the formula above. The needed theorem is a
certification statement:

```text
if stationarity_norm <= eps_1
and Hessian margin >= eps_2
and higher chart terms are bounded,
then a local verdict survives within a radius r.
```

A standard finite-dimensional route is possible. Suppose in a graph-chart ball
`||X|| <= rho`:

```text
||Grad F_mu(0)|| <= eps,
D^2 F_mu(0) <= -lambda I        with lambda > 0,
||D^2 F_mu(X) - D^2 F_mu(0)|| <= L ||X||.
```

Then inside the radius where `L ||X|| <= lambda/2`, the objective is strongly
concave. Any local maximizer in that ball is displaced from the candidate
plane by at most approximately

```text
||X_*|| <= 2 eps / lambda.
```

This is the right shape of a near-branch stability theorem: it certifies a
nearby optimizer, not that the original nonstationary plane is itself the
optimizer.

What remains open is making this certificate sharp and module-grade for the
weighted frontier. A conservative family-norm bound now exists in
`audit/0_3_2_projector_bounds_note.md`; the remaining work is to fix the
public theorem conventions and reduce pessimism:

- choose the graph-chart norm and radius;
- decide which family-norm bound should be public rather than only
  research-grade;
- handle saddle or degenerate exact-branch Hessians;
- decide whether the certificate should be local around a declared plane or
  around a projected exact branch.

This remains the main heavy research item, but the first and second variation
formulas are now in hand.

## 6. Residual-Margin Extensions

### Finite Candidate Theorem

Let candidates `alpha` have exact quadratic scores `q_alpha` and unknown
residuals `r_alpha` satisfying

```text
|r_alpha| <= R_alpha.
```

Candidate `a` is certified to beat candidate `b` at the fuller-law level if

```text
q_a - q_b > R_a + R_b.
```

For a finite candidate set, `a` is certified globally best if

```text
q_a - q_b > R_a + R_b
```

for every `b != a`.

This is the direct extension of the existing scalar residual-margin doctrine.

### Ladder Version

For a declared observer ladder `P_1,...,P_N`, with scores `F(P_j)` and
residual bounds `R_j`, observer `P_a` is certified best in the declared ladder
if

```text
F(P_a) - F(P_b) > R_a + R_b
```

for all `b != a`.

This says nothing about observers outside the declared ladder.

### Ensemble Version

If the score is an ensemble mean and `hat F_j` is an estimator, a safe bound
has the form

```text
|F_j - hat F_j| <= sampling_bound_j + model_residual_bound_j.
```

Then the same margin theorem applies with

```text
R_j = sampling_bound_j + model_residual_bound_j.
```

The hard research problem is deriving honest sampling bounds under the
sampling scheme. The inequality itself is immediate.

## 7. Local-Geometry Ensemble Margin

### Setup

Let `H` be a random local Hessian/Fisher object sampled from an ensemble
`nu`. For a declared observer `P` or observation map `C`, let `s_P(H)` be a
samplewise exact quadratic score, for example a visible score, leakage score,
clock, log-det visible summary, or declared frontier score.

The exact samplewise layer is:

```text
H -> Phi_C(H)
```

and, only when a ceiling is supplied,

```text
H,T -> Lambda_T(Phi_C(H)), tau_T(H).
```

The ensemble score is

```text
S_P = E_nu s_P(H).
```

### Finite-Sample Margin Theorem

Suppose for each observer `P_j` we estimate `S_j` by the empirical mean
`hat S_j` from `N` independent local samples. If we have a valid confidence
bound

```text
|S_j - hat S_j| <= B_j(delta)
```

simultaneously for all declared observers with probability at least
`1-delta`, then observer `a` is certified best in the declared ensemble
comparison whenever

```text
hat S_a - hat S_b > B_a(delta) + B_b(delta)
```

for every `b != a`.

This is the same residual-margin theorem with sampling error as the residual.

### Example Bound

If each score is known to lie in `[l_j, u_j]`, Hoeffding plus a union bound gives

```text
B_j(delta) = (u_j-l_j) sqrt(log(2M/delta)/(2N)),
```

for `M` declared observers. This is not always tight, but it is theorem-grade
when the boundedness assumption is true.

### What This Does Not Do

This does not infer mixture mass, cumulants, skew, hidden wells, or probability
support events. It only certifies comparisons of declared exact samplewise
quadratic observables under an honest sampling model.

### Research Consequence

The ensemble layer is not a generic non-Gaussian law engine. It is an exact
samplewise local-quadratic pushforward plus ordinary statistical certification
over the sample distribution. This is enough to make it useful, but the
sampling assumptions must be explicit.

## 8. Research Priority Ranking

### Ready For Formal Note

1. Low-rank covariance perturbation theorem.
2. Full-rank degeneracy of variable-rank frontier selection.
3. Residual-margin finite-candidate and declared-ladder theorem.
4. Affine-hidden stage-shift sign identity, with coordinate/measure caveat.
5. Ensemble finite-sample margin theorem, conditional on supplied sampling
   bounds.

### Research-Ready But Needs Definitions

1. Guarded fibre-dominance diagnostic.
2. Branch-verdict reversal diagnostic for affine-hidden fibres.
3. Declared-ladder selection with dimension cost.

### Still Heavy

1. General graph-chart frontier Hessian at declared observers. The
   non-exact stationary pass in
   `audit/0_3_2_non_exact_stationary_frontier_note.md` shows why this must be
   the parent object: exact branches are sufficient for stationarity but not
   necessary, and the exact-branch proxy can get the second-variation sign
   wrong when off-block terms cancel in first variation.
2. Declared local frontier certificate sharpness. A research-grade version is
   recorded in `audit/0_3_2_declared_frontier_certificate_note.md`, but its
   usefulness is limited by conservative `L_cert` radii.
3. Sharper near-branch Hessian-Lipschitz constants for a practical certificate.
   A first rho-dependent improvement is recorded in
   `audit/0_3_2_projector_bounds_note.md`, but the certificate is still
   conservative.
4. Sharp ensemble concentration and robust branch selection beyond bounded
   score assumptions.
5. Any Grassmannian optimizer or noncommuting global frontier policy.

## 9. Recommended Next Research Moves

1. Write a TeX theorem note for the low-rank covariance perturbation theorem.
2. Write a TeX theorem note for declared-ladder observer selection and
   full-rank degeneracy.
3. Write a short affine-hidden addendum section on staged screening signs,
   but include the measure-convention caveat.
4. Keep fibre dominance in Markdown/audit until the sample-measure convention
   is settled.
5. Keep near-branch perturbative Hessian work in research only. A conservative
   stability-certificate candidate is now available, but its projector
   derivative constants need formalization and are likely too pessimistic for a
   public diagnostic.
6. Add a general graph-chart frontier Hessian theorem note before designing
   any approximate branch API. It should recover the exact-branch Hessian as
   the invariant special case and include stationary non-exact examples.
7. Treat declared local frontier certificates as the natural future surface,
   not approximate exact-branch Hessians. Do not promote until the conservative
   Lipschitz semantics are accepted.
8. Add ensemble certification only when a sampling model and score bounds are
   supplied.

## Final Position

The clean wins are real and not too heavy. The rank-k theorem and declared
observer-ladder theorem are especially crisp. The affine-hidden sign rule is
also exact, but its interpretation needs a fixed measure convention.

The main weakness is broader than near-branch analysis. We now understand why:
a near-exact branch is not merely an exact branch with a small Hessian error.
It usually has first-order stationarity drift. More importantly, the
weighted-frontier score can have stationary non-exact observers by
first-variation cancellation. That turns the problem from a simple branch
diagnostic into a general graph-chart critical-point theorem plus perturbation
certificate.

For 0.3.2, the right research posture is:

- bank the theorem-grade wins;
- keep the near-branch and ensemble layers honest as frontier work;
- do not implement any new surface until the theorem text and residual
  semantics are settled.
