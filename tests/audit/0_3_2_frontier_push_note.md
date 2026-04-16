# 0.3.2 Frontier Push Note

Date: 2026-04-11

Status: pure research note. This extends
`audit/0_3_2_theory_research_note.md` without implementing module changes.

Further specialization: `audit/0_3_2_weighted_frontier_certificate_note.md`
derives a conservative explicit Hessian-Lipschitz certificate for the weighted
frontier graph chart.

Sharper projector constants: `audit/0_3_2_projector_bounds_note.md` replaces
the crude graph-chart constants with rho-dependent bounds and records
non-vacuity stress tests.

General graph-Hessian refinement:
`audit/0_3_2_general_graph_hessian_note.md` now records the broader
declared-observer first/Hessian calculus, and
`audit/0_3_2_graph_hessian_invariance_check.py` checks arbitrary-frame
invariance. This makes exact-branch Hessian a strict invariant sector rather
than the parent object for all frontier criticality.

## Executive Additions

This pass pushes three items further, with later notes adding the general
graph-Hessian refinement:

1. an abstract near-branch stability lemma with a proof;
2. exact dimension-cost intervals for declared-ladder observer selection;
3. an affine-hidden fibre-volume branch-reversal criterion.

The main outcome: the frontier transition layer is now sharply localized. The
raw first and second variation formulas are known, exact-branch Hessian is a
closed sector of the general graph-Hessian object, and the remaining work is
to specialize a Hessian-Lipschitz bound for the weighted frontier graph chart.

## 1. Abstract Near-Branch Stability Lemma

### Lemma

Let `f` be a `C^2` scalar objective on a Euclidean ball `B_rho(0)` in
`R^d`. Assume:

```text
||grad f(0)|| <= eps,
D^2 f(0) <= -lambda I,        lambda > 0,
||D^2 f(x) - D^2 f(0)|| <= L ||x||
```

for all `x in B_rho(0)`. Let

```text
r0 = min(rho, lambda/(2L))
```

with the convention `lambda/(2L)=infinity` if `L=0`. If

```text
4 eps / lambda < r0,
```

then `f` has a unique local maximizer `x_*` in `B_{r0}(0)`, and

```text
||x_*|| <= 2 eps / lambda.
```

Moreover, `f` is strongly concave on `B_{r0}(0)` with curvature at least
`lambda/2`.

### Proof

Set `g=-f`. Then

```text
||grad g(0)|| <= eps,
D^2 g(0) >= lambda I,
||D^2 g(x) - D^2 g(0)|| <= L ||x||.
```

For `||x|| <= r0`, the Lipschitz condition gives

```text
D^2 g(x) >= (lambda - L ||x||) I >= (lambda/2) I.
```

So `g` is strongly convex on `B_{r0}(0)`.

For any `x` with `||x||=r`, `r <= r0`, Taylor's theorem with the strong
convexity lower bound gives

```text
g(x) >= g(0) + <grad g(0),x> + (lambda/4)||x||^2
     >= g(0) - eps r + (lambda/4) r^2.
```

Choose any `r` with

```text
4 eps / lambda < r < r0.
```

Then `g(x)>g(0)` on the boundary `||x||=r`, so the minimizer of `g` over the
closed ball `||x||<=r` lies in the interior. Hence it is a critical point:

```text
grad g(x_*)=0.
```

Strong convexity makes this critical point unique in `B_{r0}`. Finally, strong
monotonicity of `grad g` gives

```text
(lambda/2)||x_*|| <= ||grad g(x_*) - grad g(0)|| = ||grad g(0)|| <= eps,
```

hence

```text
||x_*|| <= 2 eps / lambda.
```

Since `g=-f`, this unique local minimizer of `g` is the unique local maximizer
of `f` in the ball.

### Interpretation For Near-Branch Geometry

Use `f = F_mu` in a graph chart around a candidate observer plane. Then:

- `eps` is the stationarity residual;
- `lambda` is the positive local maximum margin, i.e. `-D^2F_mu(0) >= lambda I`;
- `L` is a Hessian-Lipschitz constant in the graph chart;
- the result certifies a nearby optimizer within `2 eps/lambda`.

It does not certify that the original plane is itself optimal unless
`eps=0`.

### What Remains For The Weighted Frontier

The weighted frontier already has:

- exact first variation;
- exact directional second variation;
- exact-branch Hessian formula when the off-blocks vanish.

Further stress testing in
`audit/0_3_2_non_exact_stationary_frontier_note.md` shows that exact branches
are not the only stationary points of the weighted-frontier score. Off-block
terms can cancel in the first variation, and the off-block quadratic terms can
then change the second-variation sign. Therefore the exact-branch Hessian
should be treated as a closed invariant sector, not as the parent object for
all frontier criticality.

The remaining module-grade ingredient is not the existence of an explicit
bound, since `audit/0_3_2_weighted_frontier_certificate_note.md` now gives a
conservative certificate candidate and
`audit/0_3_2_projector_bounds_note.md` sharpens its graph-chart constants.
The remaining issue is a formal, norm-fixed, and less pessimistic bound on `L`
in terms of:

- `mu`;
- weights;
- operator norms or Frobenius norms of the family;
- graph-chart radius;
- visible and hidden dimensions.

This is now a concrete bounded-derivative problem, not an unclear conceptual
gap.

## 2. Declared-Ladder Dimension-Cost Intervals

### Setup

Let a declared finite observer ladder be `P_1,...,P_N`. Let

```text
s_j = F_mu(P_j),
r_j = rank(P_j).
```

For a dimension cost `c >= 0`, define

```text
score_j(c) = s_j - c r_j.
```

### Proposition

Candidate `a` is a maximizer of the costed declared-ladder objective iff `c`
lies in the interval

```text
I_a = intersection over b != a of I_{a,b},
```

where

```text
I_{a,b} = { c >= 0 : s_a - c r_a >= s_b - c r_b }.
```

Equivalently:

- if `r_a = r_b`, then `I_{a,b}` is all `c>=0` when `s_a >= s_b`, and empty
  otherwise;
- if `r_a > r_b`, then

```text
c <= (s_a - s_b)/(r_a - r_b);
```

- if `r_a < r_b`, then

```text
c >= (s_a - s_b)/(r_a - r_b),
```

with the inequality interpreted together with `c>=0`.

### Proof

The comparison between `a` and `b` is exactly

```text
s_a - c r_a >= s_b - c r_b.
```

Rearranging gives

```text
c(r_b-r_a) >= s_b-s_a.
```

The three cases above follow by the sign of `r_b-r_a`. Intersecting all
pairwise comparison intervals gives exactly the set of costs for which `a`
beats every other declared candidate.

### Consequence

The declared-ladder selection problem is not just "pick a cost." It has an
exact phase diagram in the scalar cost `c`. This is useful research output:

- it makes the rank-budget/cost policy explicit;
- it identifies cost ranges where each observer is optimal;
- it exposes empty intervals where a candidate is never optimal;
- it avoids global Grassmannian claims.

This theorem is finite and exact for a declared ladder. It says nothing about
observers outside the ladder.

Numerical sanity check: 100 random six-observer ladders checked over a cost
grid produced `0` interval/maximizer mismatches.

## 3. Fibre-Volume Branch-Reversal Criterion

### Setup

In the affine-hidden exact sector, define for branch/candidate `j`:

```text
Var_j = A_j - 1/2 J_j^T D_j^{-1} J_j,
Fib_j = 1/2 log det D_j,
Vis_j = Var_j + Fib_j.
```

Lower action is better.

### Two-Branch Criterion

Suppose branch `a` wins variationally against branch `b`:

```text
Var_a < Var_b.
```

Define the variational gap

```text
gamma_ab = Var_b - Var_a > 0.
```

Then branch `b` beats branch `a` in the exact visible action iff

```text
Fib_a - Fib_b > gamma_ab.
```

Equality gives a visible-action tie.

### Proof

Branch `b` wins exactly when

```text
Vis_b < Vis_a
```

i.e.

```text
Var_b + Fib_b < Var_a + Fib_a.
```

Rearranging gives

```text
Fib_a - Fib_b > Var_b - Var_a = gamma_ab.
```

### Multi-Branch Criterion

If `a` is the variational winner over a finite branch set `B`, then `a`
remains the exact visible-action winner iff for every `b != a`,

```text
Fib_a - Fib_b <= Var_b - Var_a.
```

It is strictly certified if the inequalities are strict in the robust
direction:

```text
Fib_a - Fib_b < Var_b - Var_a
```

for all `b != a`.

### Residual-Robust Version

If `Var_j` and `Fib_j` are known only up to branchwise residuals

```text
|delta Var_j| <= R^V_j,
|delta Fib_j| <= R^F_j,
```

then branch `a` is robustly preserved against `b` if

```text
Var_b - Var_a
>
(Fib_a - Fib_b) + R^V_a + R^V_b + R^F_a + R^F_b.
```

Branch `b` is robustly reversed against `a` if

```text
Fib_a - Fib_b
>
Var_b - Var_a + R^V_a + R^V_b + R^F_a + R^F_b.
```

### Consequence

The branch-reversal diagnostic is more directly meaningful than a fibre
dominance ratio. A large fibre norm is not the same as a changed verdict. The
verdict changes exactly when the fibre-volume branch difference exceeds the
variational branch gap.

This is a clean theorem-grade addition for the affine-hidden exact sector.

Numerical sanity check: 1000 random five-branch examples produced `0`
preservation-criterion mismatches against direct exact visible-action
minimization.

## 4. Updated Frontier Status

### Theorem-Ready After This Pass

1. Low-rank covariance perturbation theorem.
2. Full-rank degeneracy of variable-rank frontier selection.
3. Declared-ladder dimension-cost interval theorem.
4. Fibre-volume branch-reversal criterion.
5. Residual and ensemble finite-candidate margin theorem.
6. Abstract near-branch stability lemma.

### Still Frontier

1. Sharp Hessian-Lipschitz constants for the weighted frontier graph chart.
   A conservative explicit certificate now exists in
   `audit/0_3_2_weighted_frontier_certificate_note.md`, and a sharper
   rho-dependent projector-bound pass is in
   `audit/0_3_2_projector_bounds_note.md`.
2. Sharp ensemble concentration for unbounded scores such as log-det clocks
   near support boundaries.
3. Global noncommuting Grassmannian optimization.
4. Coordinate-free interpretation of affine-hidden stage-shift signs beyond a
   fixed hidden measure convention.

## Final Position

The research frontier has narrowed. Several items that looked like "planning"
are now theorem-grade:

- dimension-cost observer-ladder selection has an exact interval calculus;
- fibre-volume branch reversal has an exact gap criterion;
- near-branch stability has an abstract theorem.

The main remaining hard technical item is no longer existence of an explicit
Hessian-Lipschitz certificate candidate; a conservative one is now available in
`audit/0_3_2_weighted_frontier_certificate_note.md`. The remaining issues are
formalization and sharpness: the constants are deliberately pessimistic, so
further work should try to prove them line-by-line, reduce them, or replace
them with verified local derivative bounds.
