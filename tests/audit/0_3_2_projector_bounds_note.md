# 0.3.2 Graph-Chart Projector Bounds Note

Date: 2026-04-11

Status: pure research note. This sharpens the projector derivative part of
`audit/0_3_2_weighted_frontier_certificate_note.md`. It does not modify module
source and is not an implementation plan.

## Executive Result

For the graph-chart projector

```text
Y(X) = [I; X],
R(X) = (I + X^T X)^{-1},
P(X) = Y(X) R(X) Y(X)^T,
```

on the Frobenius ball `||X||_F <= rho <= 1`, with unit Frobenius tangent
directions, a conservative but much sharper derivative-bound family is:

```text
y(rho)  = sqrt(1 + rho^2),
C1(rho) = 2 y + 2 rho y^2,
C2(rho) = 2 + 8 rho y + (8 rho^2 + 2) y^2,
C3(rho) = (24 rho + 48 rho^3) y^2
          + 6(8 rho^2 + 2)y
          + 12 rho.
```

At `rho=1` this gives approximately:

```text
C1 = 6.8285,
C2 = 33.3137,
C3 = 240.8528,
```

instead of the earlier coarse constants:

```text
8, 64, 512.
```

The constants are still conservative, but they are now line-by-line traceable.

## 1. Basic Bounds

Let `||X||_F <= rho <= 1` and define

```text
S(X)=I+X^T X,
R(X)=S(X)^{-1},
Y(X)=[I;X].
```

Then

```text
||Y||_op <= y = sqrt(1+rho^2),
||R||_op <= 1.
```

For unit Frobenius directions `H,K,L`:

```text
DY[H] = [0;H],
||DY[H]||_F <= 1,
D^2Y = 0.
```

Also

```text
DS[H] = H^T X + X^T H,
D^2S[H,K] = H^T K + K^T H,
D^3S = 0.
```

Using `||H||_op <= ||H||_F = 1`, the following bounds hold:

```text
||DS[H]||_F <= 2 rho,
||D^2S[H,K]||_F <= 2.
```

## 2. Inverse Bounds

Since `R=S^{-1}`,

```text
DR[H] = -R DS[H] R.
```

Thus

```text
||DR[H]||_F <= 2 rho.
```

The second derivative is

```text
D^2R[H,K]
= R DS[K] R DS[H] R
 + R DS[H] R DS[K] R
 - R D^2S[H,K] R.
```

Therefore

```text
||D^2R[H,K]||_F <= 8 rho^2 + 2.
```

For the third derivative, `D^3S=0`. The inverse derivative has:

- six ordered terms containing three first derivatives `DS`;
- six ordered terms containing one `D^2S` and one `DS`.

Hence

```text
||D^3R[H,K,L]||_F
<= 6(2rho)^3 + 6(2)(2rho)
= 48 rho^3 + 24 rho.
```

## 3. Projector Derivative Bounds

Differentiate

```text
P = Y R Y^T.
```

### First Derivative

The first derivative has three terms:

```text
DP[H] = DY[H] R Y^T + Y DR[H] Y^T + Y R DY[H]^T.
```

Thus

```text
||DP[H]||_F
<= 2 ||DY[H]||_F ||R||_op ||Y||_op
   + ||Y||_op^2 ||DR[H]||_F
<= 2y + 2rho y^2.
```

So

```text
C1(rho)=2y+2rho y^2.
```

### Second Derivative

The second derivative has:

- two terms with both derivatives landing on the two `Y` factors;
- four terms with one derivative on a `Y` factor and one on `R`;
- one term with both derivatives on `R`.

Therefore

```text
||D^2P[H,K]||_F
<= 2
 + 4 y ||DR||_F
 + y^2 ||D^2R||_F.
```

Using the inverse bounds:

```text
||D^2P[H,K]||_F
<= 2 + 8rho y + (8rho^2+2)y^2.
```

So

```text
C2(rho)=2+8rho y+(8rho^2+2)y^2.
```

### Third Derivative

Because `Y` is linear, the third derivative has:

- one term with all three derivatives on `R`;
- six terms with two derivatives on `R` and one on a `Y` factor;
- six terms with one derivative on `R` and one derivative on each `Y` factor.

Hence

```text
||D^3P[H,K,L]||_F
<= y^2 ||D^3R||_F
 + 6 y ||D^2R||_F
 + 6 ||DR||_F.
```

Substituting the inverse bounds:

```text
C3(rho)
= (24rho + 48rho^3)y^2
 + 6(8rho^2+2)y
 + 12rho.
```

## 4. Improved Weighted-Frontier Certificate Constant

Keep the projector-functional derivative bounds from
`audit/0_3_2_weighted_frontier_certificate_note.md`:

```text
G1 = [3(1+mu)sqrt(m) + mu sqrt(n)] A2,
G2 = 6(1+mu) A2,
G3 = 6(1+mu) A2,
A2 = sum_i w_i ||A_i||_op^2.
```

Then the graph-chart Hessian-Lipschitz certificate may use

```text
L_cert(rho) =
  G3 C1(rho)^3
  + 3 G2 C2(rho) C1(rho)
  + G1 C3(rho).
```

For `rho=1`, this improves the constants from the earlier `8,64,512` bound
while preserving the same proof architecture.

## 5. Numerical Sanity Check

A randomized directional finite-difference stress check over graph-chart
points with `||X||_F <= rho <= 1` found no violations of the sharper bounds.
The largest observed directional derivatives in the previous stress pass were
approximately:

```text
max ||dP||  = 1.4143,
max ||d2P|| = 2.8282,
max ||d3P|| = 8.4843.
```

These are far below the sharpened `rho=1` constants, so the constants remain
comfortably conservative.

The more extensive brutality check in
`audit/0_3_2_certificate_brutality_check.py` tested

```text
rho in {0.02, 0.05, 0.1, 0.25, 0.5, 1.0}
```

with random dimensions and random unit directions. It again found `0`
violations with a 2 percent numerical slack. The observed maxima stayed near:

```text
max ||dP||  ~= 1.4143,
max ||d2P|| ~= 2.8285,
max ||d3P|| ~= 8.4853.
```

The sharpened constants are still conservative, especially for `D^3P`, but the
rho-dependent proof path is consistent with numerical stress tests.

## 6. Certificate Non-Vacuity

The same brutality check compared the old constants `(8,64,512)` against the
rho-dependent constants. In a synthetic exact-branch-plus-off-block family, the
old-to-new `L_cert` ratio was about:

```text
old_L_cert / new_L_cert ~= 2.06.
```

The sharpened certificate accepted a case with off-block scale

```text
delta = 3e-5
```

where the old certificate failed, while both failed by

```text
delta = 1e-4.
```

So the improvement is not just cosmetic: it roughly doubles the certified
radius in the tested family. However, the condition remains stringent because
the certificate is a sufficient local theorem. It requires stationarity drift
small compared to the Hessian margin and the global chart Lipschitz bound.

## 7. Remaining Caution

The note uses Frobenius norm for tangent directions and output matrices. The
constants should not be moved into a public theorem without keeping that norm
contract explicit.

The constants are sufficient, not sharp. They likely still make many
near-branch certificates vacuous. Their purpose is to establish a rigorous
research bridge; practical diagnostics may need local derivative estimates or
better chart constants.

## Final Position

The weighted-frontier near-branch certificate is stronger now:

- the projector derivative constants are rho-dependent;
- the proof is line-by-line rather than a sketch with rounded constants;
- the Hessian-Lipschitz constant is still conservative but less crude.

The next research frontier is practical non-vacuity: quantify when
`L_cert(rho)` is small enough relative to `lambda^2/eps` to certify a real
nearby branch in realistic families.
