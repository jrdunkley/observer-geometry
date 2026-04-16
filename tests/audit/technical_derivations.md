# Technical Derivations and Counterexamples

## 1. Why `Phi_C(H)` Is Exact Without Gaussianity

Let `H in SPD(n)` and `C:R^n -> R^m` be surjective. For fixed visible displacement `y`, minimize

```text
E(x)=x^T H x subject to Cx=y.
```

The Lagrangian `x^T H x - 2 lambda^T(Cx-y)` gives

```text
Hx = C^T lambda,
x = H^{-1} C^T lambda,
CH^{-1}C^T lambda = y.
```

Since `C` is surjective and `H>0`, `CH^{-1}C^T>0`. Thus

```text
lambda = (CH^{-1}C^T)^{-1} y = Phi_C(H)y,
x_* = H^{-1}C^T Phi_C(H)y.
```

The minimum value is

```text
x_*^T H x_* = y^T Phi_C(H) y.
```

No probability law entered. This is exact Hilbert/positive-cone quotient geometry.

If `X~N(0,H^{-1})`, then `CX` is Gaussian with covariance `CH^{-1}C^T`, so its precision is `Phi_C(H)`. This last sentence is the Gaussian full-law upgrade.

## 2. Why The Local Calculus Is Exact But Only Local

Write `K(t)=C(H+t Delta)^{-1}C^T` and `Phi(t)=K(t)^{-1}`. Using

```text
d(H^{-1}) = -H^{-1} Delta H^{-1},
d^2(H^{-1})[Delta,Delta] = 2H^{-1}Delta H^{-1}Delta H^{-1},
d(K^{-1}) = -K^{-1}(dK)K^{-1},
```

one obtains

```text
Phi_C(H+t Delta) = Phi + t L^T Delta L - t^2 L^T Delta P H^{-1} Delta L + O(t^3),
```

where

```text
L = H^{-1}C^T Phi,
P = I - LC.
```

The quadratic defect has the Gram form

```text
Q = (H^{-1/2}P^T Delta L)^T (H^{-1/2}P^T Delta L) >= 0.
```

Again this is a theorem about the SPD matrix path `H(t)`. If `H(t)` is the Hessian/Fisher path of a non-Gaussian model, this is its local quotient-Hessian calculus. It is not the Taylor expansion of the entire pushed-forward density unless a separate theorem identifies those objects.

## 3. Same Hessian Does Not Determine A Non-Gaussian Visible Law

Let

```text
p_G(x) = (2pi)^(-1/2) exp(-x^2/2),
p_T(x) = t_nu(x; scale=sqrt((nu+1)/nu)), nu=5.
```

For the Student density

```text
log p_T(x) = const - (nu+1)/2 log(1 + x^2/(nu s^2)).
```

Therefore

```text
-d^2/dx^2 log p_T(0) = (nu+1)/(nu s^2) = 1
```

by the scale choice. The Gaussian has the same local Hessian. Any one-dimensional quadratic object therefore agrees.

But the full laws differ. The audit script computes symmetric KL `0.255837` and squared Hellinger `0.0195813`.

Conclusion: local Hessian exactness does not imply full-law exactness.

## 4. Observer Ranking Reversal By A Quartic Tail

Define a product pair `(X,Y)` with baseline marginals

```text
X0 ~ N(0,1),
Y0 density proportional to exp(-y^2/2).
```

Alternative marginals:

```text
X1 ~ N(0, 1/1.08),
Y1 density proportional to exp(-y^2/2 - 0.35 y^4).
```

At the mode, the Hessian change is

```text
Delta_H = diag(0.08, 0).
```

Thus any rank-one quadratic observer selection based on local Hessian change ranks the `X` observer above the `Y` observer.

The full law does the opposite. The script computes

```text
symKL(X0,X1) = 0.00296296,
symKL(Y0,Y1) = 0.904204.
```

Conclusion: full-law observer ranking can reverse while the local quadratic calculation is exact in its own layer.

## 5. Quadratic Closure Does Not Imply Full-Law Closure

Consider the two-dimensional density

```text
p_t(x,z) ∝ exp[-0.5(x^2+z^2) - alpha(x^4+z^4) - t x z^2],
alpha = 0.15.
```

The term `t x z^2` is cubic. Its first and second derivatives at `(0,0)` vanish in the Hessian block sense relevant to the local quadratic precision. Hence `H_t(0)=H_0(0)` and the quadratic quotient layer predicts no visible response and no closure leakage for observing `x`.

The actual marginal is

```text
p_t^vis(x) ∝ exp[-0.5x^2 - alpha x^4] * integral exp[-0.5z^2 - alpha z^4 - t x z^2] dz.
```

The integral depends on `x` when `t != 0`. Therefore the visible law changes even though the local Hessian did not. Numerically, for `t=0.28`, the visible marginal symmetric KL is `0.0139792`, and the mean moves from `0` to `-0.087432`.

Conclusion: the exact quadratic closure criterion `Q=0` is not a full non-Gaussian marginal closure theorem.

## 6. Matrix Support Strata Are Not Probability Support

Let `p` be `N(0,1)` and let `q` be `N(0,1)` conditioned on `x>=-2`. Around the interior point `0`,

```text
log q(x) = log p(x) - log P(N(0,1)>=-2),
```

so the local Hessian is exactly the same as the Gaussian Hessian.

But

```text
P_p(x<-2)=0.0227501,
P_q(x<-2)=0,
KL(p||q)=infinite.
```

Conclusion: support/restart objects in the module classify ranks and kernels of supplied PSD matrix paths, not hard support boundaries of arbitrary distributions.

## 7. Gaussian Gluing Is A Law-Level Closure Convenience

For centred Gaussian variables, pair laws are determined by their covariance blocks. Therefore common Gaussian gluing on a graph reduces to covariance/correlation completion plus variance consistency.

For non-Gaussian variables, covariance completion does not determine joint law compatibility. Distinct non-Gaussian edge laws can share the same covariance/correlation data and have incompatible higher moments, tail dependence, atoms, or support constraints.

Conclusion: Bell/temporal/graph Gaussian gluing should remain labelled `G`, even though the underlying Schur/correlation geometry supplies useful shadows.

## 8. What Would Be Needed For Higher-Order Observation Geometry

A plausible order-three extension would have to start with a local log-density expansion

```text
-log p(x) = const + 1/2 H_ij x_i x_j + 1/6 T_ijk x_i x_j x_k
           + 1/24 U_ijkl x_i x_j x_k x_l + ...
```

The quotient operation at order two is the constrained minimisation over hidden fibres. At order three and above, eliminating hidden variables would require either:

- constrained minimisation of a non-quadratic jet, producing visible tensors and hidden-response tensors; or
- integration over hidden variables, producing cumulant corrections that are not the same operation as minimisation unless a Laplace/small-noise limit is specified.

This split is the real next frontier. Gaussian theory avoids it because the quadratic jet is the whole law.
