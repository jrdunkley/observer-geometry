# Non-Gaussian Survival Library

These examples show where the local quadratic layer remains useful outside exact Gaussian full-law mode. They are not promoted to universal non-Gaussian theorems; they identify conditions under which higher-order structure does not change the observational verdict in the tested family.

## S1. Elliptic Heavy-Tail Scale Perturbations

Construction: use univariate Student `t_7` laws with scale chosen so that the local Hessian is `h`. Compare a baseline `h=1` to `h=1.30` and `h=1.08`.

Quadratic prediction: the `h=1.30` perturbation is stronger than the `h=1.08` perturbation.

Full visible law: the same ranking holds. The script gives symmetric KL `0.0240971` for the strong perturbation and `0.00207309` for the weak perturbation, a ratio about `11.62`.

Why it survives: the family is fixed-shape and elliptic/heavy-tailed; changing the local scale also changes the full law in a monotone way for this comparison. Gaussian exactness is not needed to get the qualitative ordering, though it is needed for Gaussian closed-form divergences.

## S2. Small Non-Gaussian Quartic Family

Construction: use normalized densities

```text
p_h(x) ∝ exp[-0.5 h x^2 - 0.03 x^4].
```

Compare baseline `h=1` to `h=1.30` and `h=1.08`.

Quadratic prediction: the stronger Hessian perturbation should outrank the weaker perturbation.

Full visible law: the same ranking holds. The script gives symmetric KL `0.0205057` for the strong perturbation and `0.00164434` for the weak perturbation, a ratio about `12.47`.

Why it survives: the higher-order term is shared and small, so it does not dominate the Hessian perturbation.

## S3. Pure Positive-Cone Algebra

Construction: choose any SPD form `H`, surjective linear map `C`, fixed ceiling `T`, visible form `X` with `0<X<=T` on `Ran(T)`, and hidden load `Lambda=T^{1/2}X^{-1}T^{1/2}-I`.

Quadratic prediction: `Phi_C`, hidden load, rank, determinant clock, support-stratum transport, and closure scores are well-defined.

Full visible law: no full-law claim is made.

Why it survives: the statement is not probabilistic. It is exact operator geometry on a supplied local Hessian/Fisher/SPD object.

## S4. Regular Fisher Quotients

Construction: a statistical model factors locally through a regular structural quotient `q:Theta -> barTheta`; the experiment law is constant on each structural fibre.

Quadratic prediction: vertical scores vanish, and the Fisher form descends to the quotient tangent.

Full visible law: this is not Gaussian-specific, but it is local Fisher geometry, not a claim that `Phi_C(H)` computes arbitrary marginal law divergences.

Why it survives: Fisher information is itself a local quadratic object. The exactness is at order two in parameter space.

## S5. KL Chain Rule and Total Covariance

Construction: ordinary probability spaces with a visible variable and a hidden variable.

Quadratic prediction: not needed. KL marginal descent and total covariance decomposition have exact law-level forms.

Full visible law: exact beyond Gaussian when the relevant integrability and absolute-continuity assumptions hold.

Why it survives: these are broader Hilbert/probability projection identities. They support the claim that the broader programme is not merely Gaussian, but they are not the same as closed-form Gaussian precision geometry.

## S6. Fixed-Precision Affine-Hidden Sector

Construction: choose a law

```text
p(v,h) proportional to exp(-A(v)-1/2 h^T D h-J(v)^T h)
```

with `D` fixed SPD and arbitrary visible action/coupling satisfying the usual integrability assumptions.

Quadratic prediction: hidden elimination acts by the Schur/variational correction `-1/2 J(v)^T D^{-1}J(v)`.

Full visible law: exactly the same visible action up to a constant:

```text
S_vis(v)=A(v)-1/2 J(v)^T D^{-1}J(v)+const.
```

Why it survives: the hidden fibre is conditionally Gaussian with fixed precision, so the determinant contribution is independent of `v`. Staged hidden elimination is exact by Schur-complement associativity. The weakpoint numerics found a branch gap of `0.0243008` in a nonlinear fixed-precision affine-hidden example, with no extra law-level correction.

## S7. Affine-Hidden Tower Elimination

Construction: split a fixed hidden precision `D` into several hidden coordinate blocks and eliminate them in any order.

Quadratic prediction: staged Schur complements should match one-step hidden elimination.

Full visible law: exact agreement in the fixed-precision affine-hidden sector.

Why it survives: block Gaussian integration and Schur complements are associative. The weakpoint stress test eliminated five hidden coordinates in all `5!` orders with maximum residual `2.78e-16`.

## S8. White Rank-One Differential-Correlation Sector

Construction: use a covariance family

```text
Sigma_eps = I + eps f f^T, f=(a,b) in visible plus hidden coordinates.
```

Quadratic prediction: the hidden gap is one-channel:

```text
T - Phi =
  eps^2 ||b||^2 / ((1+eps||a||^2)(1+eps||f||^2)) * a a^T.
```

Full visible law: in the Gaussian/Fisher covariance reading, this is exact. The same Sherman-Morrison calculation also gives the Fisher ceiling `f^T Sigma_eps^{-1} f = ||f||^2/(1+eps||f||^2)`, tending to `1/eps`.

Why it survives: this is not a generic non-Gaussian law theorem; it is an exact rank-one covariance/Fisher perturbation sector. The Claude-notes follow-up script verified the hidden-gap formula to `5.03e-17` and rank one numerically.

## S9. Residual-Margin Branch Survival

Construction: compare two quadratic branch scores with gap `gamma`, while law-level or modelling residuals are only known to be bounded by `R`.

Quadratic prediction: if `gamma > 2R`, the branch ordering cannot reverse under those residuals.

Full visible law: the conclusion is exact as an ordering theorem whenever the residual bound is valid. It does not identify the residual; it tells how much unmodelled law structure the quadratic verdict can tolerate.

Why it survives: it is a margin argument, not a Gaussian closure claim. The follow-up script used `gamma=0.08`; with `R=0.03`, no reversal is possible, while with `R=0.05`, a constructive reversal exists and random residual draws reversed the verdict `953` times in `50000` trials.
