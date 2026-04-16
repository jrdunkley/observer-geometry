# 0.3.2 Weighted-Frontier Near-Branch Certificate Note

Date: 2026-04-11

Status: pure research note. This note pushes the near-branch stability problem
one step beyond `audit/0_3_2_frontier_push_note.md`. It gives a conservative
computable sufficient condition for a nearby local optimizer in the weighted
frontier graph chart. It is not an implementation plan.

Sharper follow-up: `audit/0_3_2_projector_bounds_note.md` derives
rho-dependent graph-chart projector bounds and should be preferred over the
older coarse constants below when doing further research.

## Executive Result

The weighted-frontier near-branch problem now has a conservative certificate
shape:

```text
stationarity residual eps
local Hessian maximum margin lambda
weighted-family scale A2
chart radius rho <= 1
```

produce a coarse Hessian-Lipschitz bound `L_cert`. If

```text
4 eps / lambda < min(rho, lambda/(2 L_cert)),
```

then the graph-chart objective has a unique nearby local maximizer, and its
distance from the declared candidate plane is at most

```text
2 eps / lambda.
```

The constants below are deliberately conservative. They are research-grade
sufficient constants under the derivative-bound sketches below, not sharp
constants and not a recommended public API. The projector derivative bounds
should be formalized line-by-line before this becomes paper or module text.

## 1. Frontier Objective As A Projector Functional

For symmetric family members `A_i`, nonnegative weights `w_i`, and `mu >= 0`,
define

```text
M = sum_i w_i A_i^2,
S(P) = sum_i w_i ||P A_i P||_F^2,
F_mu(P) = S(P) - mu L(P).
```

Using the weighted energy identity

```text
Tr(P M) = S(P) + L(P),
```

the same objective is

```text
F_mu(P) = (1+mu) S(P) - mu Tr(P M).
```

As a functional of an orthogonal projector `P`,

```text
S_i(P) = ||P A_i P||_F^2 = Tr(P A_i P A_i P).
```

So `F_mu(P)` is a cubic polynomial in `P` plus a linear term.

## 2. Derivative Bounds For The Projector Functional

Let

```text
A2 = sum_i w_i ||A_i||_op^2.
```

On rank-m projectors in `R^n`, with Frobenius norm on variations:

```text
||D Phi(P)|| <= G1,
||D^2 Phi(P)|| <= G2,
||D^3 Phi(P)|| <= G3,
```

for `Phi(P)=F_mu(P)`, where the following conservative bounds hold:

```text
G1 = [3(1+mu) sqrt(m) + mu sqrt(n)] A2,
G2 = 6(1+mu) A2,
G3 = 6(1+mu) A2.
```

### Proof Sketch

For one operator,

```text
S_i(P)=Tr(P A_i P A_i P).
```

The first derivative contains three terms, each bounded by
`||A_i||_op^2 sqrt(m) ||E||_F` on rank-m projectors. The linear term
`Tr(PM)` contributes at most `||M||_F ||E||_F`, and
`||M||_F <= sqrt(n) A2`.

The second derivative of the cubic term contains six terms with two variation
slots and one bounded projector slot. Each is bounded by
`||A_i||_op^2 ||E||_F ||F||_F`. The linear term has no second derivative.

The third derivative of the cubic term contains six terms with three variation
slots and no projector slot. Again each is bounded by
`||A_i||_op^2` times the three Frobenius norms. The linear term has no third
derivative.

Summing with weights gives the stated constants.

These constants are not sharp, but they are explicit and dimensionally stable.
They should be treated as a conservative research certificate until the
projector derivative bounds below are written as a formal lemma.

## 3. Coarse Graph-Chart Projector Constants

Use the graph chart

```text
Y(X) = [I; X],
P(X) = Y(X) (I + X^T X)^{-1} Y(X)^T.
```

On the Frobenius ball

```text
||X||_F <= rho <= 1,
```

one may use the conservative derivative bounds

```text
||DP(X)|| <= C1 = 8,
||D^2P(X)|| <= C2 = 64,
||D^3P(X)|| <= C3 = 512.
```

### Derivation Sketch

Write

```text
R(X) = (I + X^T X)^{-1}.
```

For `rho <= 1`, `||Y||_op <= sqrt(2)` and `||R||_op <= 1`. The derivatives of
`S(X)=I+X^T X` satisfy, for unit Frobenius directions,

```text
||DS||_F <= 2,
||D^2S||_F <= 2,
D^3S = 0.
```

Differentiating `R=S^{-1}` gives conservative bounds

```text
||DR|| <= 2,
||D^2R|| <= 10,
||D^3R|| <= 120.
```

Differentiating the product `P=YRY^T`, noting that `D^kY=0` for `k>=2`, gives
the coarse bounds `C1=8`, `C2=64`, `C3=512`.

The constants are intentionally rounded upward. Sharp constants are not needed
for a sufficient certificate; they would matter only for a less conservative
public diagnostic.

Numerical stress check: 300 random graph-chart points with `||X||_F<=1` and
random unit directions gave approximate directional derivative maxima

```text
max ||dP||  = 1.4143
max ||d2P|| = 2.8282
max ||d3P|| = 8.4843
```

against the conservative bounds `8,64,512`, with `0` bound violations. This is
not a proof of the constants, but it is a useful sanity check on the norm
conventions and derivative formulas.

## 4. Hessian-Lipschitz Bound For The Chart Objective

Let

```text
f(X) = F_mu(P(X)).
```

By the third-order chain rule,

```text
D^3 f
= D^3 Phi[DP,DP,DP]
 + 3 D^2 Phi[D^2P,DP]
 + D Phi[D^3P].
```

Therefore a Hessian-Lipschitz bound on `||X||_F <= rho <= 1` is

```text
L_cert = G3 C1^3 + 3 G2 C2 C1 + G1 C3.
```

Substituting the coarse constants gives

```text
L_cert
= 6(1+mu) A2 * 8^3
 + 18(1+mu) A2 * 64 * 8
 + [3(1+mu)sqrt(m) + mu sqrt(n)] A2 * 512.
```

This simplifies to

```text
L_cert
= A2 [
   3072(1+mu)
 + 9216(1+mu)
 + 512(3(1+mu)sqrt(m) + mu sqrt(n))
].
```

or

```text
L_cert
= A2 [
   12288(1+mu)
 + 1536(1+mu)sqrt(m)
 + 512 mu sqrt(n)
].
```

This is coarse but computable from declared inputs.

## 5. Conservative Near-Branch Certificate

Let the graph-chart gradient at the candidate plane be `g`, and let the
chart-Hessian at the candidate plane be `Hess`.

Suppose:

```text
eps = ||g||_F,
Hess <= -lambda I,   lambda > 0,
rho <= 1,
L_cert as above.
```

If

```text
4 eps / lambda < min(rho, lambda/(2 L_cert)),
```

then the weighted-frontier objective has a unique local maximizer in the graph
ball of radius `min(rho, lambda/(2L_cert))`, and the optimizer satisfies

```text
||X_*||_F <= 2 eps / lambda.
```

### Interpretation

This certifies that a nonstationary near-branch plane has a nearby true local
maximizer. It does not certify that the original plane itself is optimal.

If `eps=0`, the candidate plane is stationary and the displacement bound gives
`X_*=0`, recovering the exact-branch local maximum situation.

If `lambda <= 0`, the certificate is unavailable. This covers saddle,
degenerate, and local-minimum cases.

If the condition fails, that is not a proof of instability. It only means this
coarse sufficient certificate is too weak or inapplicable.

## 6. Practical Meaning

This result makes near-branch research concrete:

- compute or estimate the first-variation norm;
- compute the local Hessian margin;
- compute a conservative `A2` family scale;
- test the sufficient inequality.

The certificate will often be pessimistic because `C1,C2,C3` are deliberately
large. A useful next theory step is to sharpen the projector derivative bounds
or compute local derivative norms directly for a given candidate plane.

## 7. Remaining Research

This note closes the main conceptual gap but not the practical sharpness gap.
The next work is:

1. continue sharpening `C1,C2,C3` for the graph chart; a first rho-dependent
   improvement is now recorded in `audit/0_3_2_projector_bounds_note.md`;
2. compare graph-chart constants with exponential/geodesic chart constants;
3. test whether local numerical derivative bounds can safely replace global
   coarse constants in an audit-only setting;
4. extend the same certificate to local minima by applying the lemma to `-F`;
5. handle degenerate Hessian cases using center-manifold or higher-order
   terms, which is a separate problem.

## Final Position

Near-branch analysis is now much less vague. We have:

- exact first variation;
- exact directional second variation;
- an abstract stability lemma;
- a conservative explicit Lipschitz-certificate candidate.

The remaining question is no longer "is a theorem possible?" It is "how sharp
and how formal can the constants be made while staying simple enough to trust?"
