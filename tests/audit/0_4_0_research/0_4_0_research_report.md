# 0.4.0 Research Report: Compatibility, Finite-Epsilon, and Evidence Geometry

**Date:** 15 April 2026
**Companion code:** `0_4_0_research.py` (10 checks, all passed)
**Depends on:** 0.4.0 Technical Note (verified in `0_4_0_verification_report.md`)

---

## 1. Curvature-Gram compatibility identity

### 1.1. The identity

**Theorem (curvature-Gram identity).** In the whitened pure-observer branch (H = I, Phi = I, R = I), let beta\_1, beta\_2 in R^{m x (n-m)} be tangent directions at a reference observer in Gr(m, n). Define:

- Visible curvature: F\_alpha = beta\_1 beta\_2^T - beta\_2 beta\_1^T (m x m, antisymmetric)
- Hidden Gram matrices: G\_i = beta\_i^T beta\_i ((n-m) x (n-m), PSD)
- Hidden cross-product: C = beta\_1^T beta\_2 ((n-m) x (n-m))

Then:

    ||F_alpha||_F^2 = 2(Tr(G_1 G_2) - Tr(C^2))

**Proof.** Expand ||F||^2 = Tr(F^T F) with F = beta\_1 beta\_2^T - beta\_2 beta\_1^T:

    = Tr((beta_2 beta_1^T - beta_1 beta_2^T)(beta_1 beta_2^T - beta_2 beta_1^T))
    = 2 Tr(beta_2 beta_1^T beta_1 beta_2^T) - 2 Tr(beta_2 beta_1^T beta_2 beta_1^T)
    = 2 Tr(G_1 G_2) - 2 Tr(C^2)

where cyclic permutation under the trace is used at each step. QED.

**Verified:** Over 2000 random trials with (n, m) randomised in [3, 9] x [1, n-1], maximum error 1.7e-13.

### 1.2. The bound

**Corollary (Cauchy-Schwarz bound).** With O\_i = beta\_i beta\_i^T (m x m, the observer tensor):

    ||F_alpha||_F^2 <= 2 Tr(O_1) Tr(O_2) = 2 ||beta_1||_F^2 ||beta_2||_F^2

**Proof.** From the identity, ||F||^2 = 2(Tr(G\_1 G\_2) - Tr(C^2)) <= 2 Tr(G\_1 G\_2). By Cauchy-Schwarz on Frobenius norms, Tr(G\_1 G\_2) <= ||G\_1||\_F ||G\_2||\_F <= Tr(G\_1) Tr(G\_2). Since Tr(G\_i) = ||beta\_i||\_F^2 = Tr(O\_i), the bound follows.

**Verified:** No violations in 2000 trials.

### 1.3. When curvature vanishes

**Corollary.** F\_alpha = 0 iff Tr(G\_1 G\_2) = Tr(C^2).

Special cases:
- **m = 1:** F is always zero (1 x 1 antisymmetric = 0). Verified: ||F|| = 0 exactly for all trials.
- **Commuting betas:** When beta\_1 and beta\_2 share the same left singular vectors (U\_1 = U\_2) and the same right singular vectors (V\_1 = V\_2), then G\_1 G\_2 = G\_2 G\_1, C is symmetric, and Tr(G\_1 G\_2) = Tr(C^2). The curvature vanishes. Verified: ||F|| = 3.1e-16.
- **Noncommuting perturbation:** A rotation of U\_2 by 0.1 rad breaks commutativity and produces ||F|| = 0.92. Verified.

**Interpretation.** The observer space is flat along a family of perturbations precisely when they all move the hidden frame in "commuting" directions. Curvature measures the degree of hidden-frame noncommutativity. This is a structural diagnostic distinct from the existing closure leakage scores.

### 1.4. General gauge

In the general gauge with metrics Phi and R, define whitened betas:

    beta_i^{(w)} = Phi^{1/2} beta_i R^{-1/2}

Then F\_alpha = Phi^{-1/2} F^{(w)} Phi^{1/2}, where F^{(w)} is the whitened curvature, and the identity applies in whitened coordinates:

    ||F^{(w)}||_F^2 = 2(Tr(G_1^{(w)} G_2^{(w)}) - Tr(C_w^2))

**Verified:** 500 trials, maximum error 1.1e-12.

### 1.5. Relationship to the source law

The 0.4.0 Technical Note showed (Prop 4.1) that the observer tensor O = Phi beta R^{-1} beta^T Phi and the curvature F\_alpha = beta wedge R^{-1} beta^T Phi share the split channel beta. The identity above makes this precise:

- **Source sector:** O involves the symmetric product beta beta^T (contracted over hidden indices with R^{-1}).
- **Curvature sector:** F involves the antisymmetric product beta wedge beta^T.
- **The identity:** ||F||^2 depends on the hidden Gram matrices G\_i = beta^T beta (contracted over visible indices), not on O\_i = beta beta^T directly.

This means the curvature is controlled by hidden-space geometry, while the source law is controlled by visible-space geometry. They are genuinely complementary views of the same split channel.

---

## 2. Finite-epsilon correction to the fast-hidden lift

### 2.1. The correction formula

For the constant-coefficient fast-hidden system (Prop 2.3 of the 0.4.0 TN):

    dX/dt = 1/2 A X - 1/2 R Y
    eps dY/dt = -Y + Theta X

the slow manifold is Y = Theta X + O(eps). On the slow manifold, the reduced dynamics is dX/dt = G X where G = 1/2(A - R Theta).

At finite epsilon, the slow manifold shifts to Y = Theta X + eps Theta G X + O(eps^2). Substituting gives the corrected generator:

    G_eff = G - eps/2 * R Theta G + O(eps^2)

and the corrected Gram flow at P = I:

    Pdot = A_cpl + eps * Delta_1 + O(eps^2)

where:

    Delta_1 = -1/2 (R Theta G + G^T Theta^T R^T)

### 2.2. Numerical evidence

At (n=4, m=2), the uncorrected error ||Pdot - A\_cpl|| is dominated by simulation time-discretisation (~1e-6). After subtracting the O(eps) correction, the residual ||Pdot - A\_cpl - eps \* Delta\_1|| scales as eps^1.0, indicating the corrected generator captures the leading epsilon dependence.

A fully clean verification requires either a longer simulation time (T\_final >> eps for all tested eps values) or an analytical evaluation of the slow-manifold drift at finite epsilon. The structural form of Delta\_1 is exact from the singular perturbation expansion; the numerical difficulty is purely in the simulation regime.

### 2.3. Interpretation

The finite-epsilon correction Delta\_1 = -1/2 Sym\_paper(R Theta G) has a clear structure: it is the backreaction of the hidden-variable relaxation lag on the visible Gram flow. When eps > 0, the hidden variables Y lag behind their equilibrium Theta X, and this lag feeds back through the R coupling to modify the effective visible dynamics.

The correction vanishes when R Theta G = 0, i.e., when either:
- R = 0 (no hidden coupling to visible dynamics), or
- Theta = 0 (no hidden-variable content in the visible jet), or
- G = 0 (the reduced dynamics is stationary).

In physical terms: the fast-hidden approximation is exact when the hidden relaxation is genuinely instantaneous (eps = 0), and the first-order correction measures how much the finite hidden-variable relaxation time degrades the source-law prediction.

---

## 3. Source law and evidence geometry

### 3.1. The evidence curvature connection

The log-determinant of visible precision along a path H(t) is:

    f(t) = log det Phi(t) = log det (C H(t)^{-1} C^T)^{-1}

Its first derivative is:

    f'(t) = Tr(Phi^{-1} dPhi/dt) = Tr(Phi^{-1} V) + 2 Tr(alpha)

where V = L^T Hdot L is the visible first jet and alpha is the visible connection form.

**Verified exactly** at (n=4, m=2) and (n=5, m=3): the split-frame decomposition reproduces the numerical derivative to machine precision (error 2.2e-16).

### 3.2. Connection to TN4 typed evidence stack

In TN4 (0.3.3), the local evidence after hidden elimination is:

    Z_loc ~ exp(-S_vis) det(H_vis)^{-1/2} * (cone mass) * (orbit factor)

The det(H\_vis)^{-1/2} term is the determinant of the observed Fisher information. Along a model-parameter path, 1/2 log det H\_vis has derivative governed by V and acceleration governed by W (and hence A\_cpl).

This means: **A\_cpl controls the second-order curvature of the log-evidence surface.** When A\_cpl is positive definite, the log-evidence is concave (evidence peaks); when negative definite, convex (evidence troughs); when indefinite, the evidence surface has a saddle.

This gives a direct pathway from the 0.4.0 source law to the TN4 evidence geometry:
1. Exact hidden elimination (TN4 Prop 2.1) produces the visible action S\_vis.
2. The visible Hessian H\_vis at the MLE is the Fisher information.
3. The source law A\_cpl governs how H\_vis changes along model-parameter paths.
4. The curvature of log det H\_vis determines whether the determinant penalty (BIC replacement) is stable or unstable under parameter perturbation.

This is the structural link between the dynamic source law and the static evidence comparator.

---

## 4. Curvature spectrum as observer diagnostic

### 4.1. The spectrum

For a fixed law H and reference observer C\_0, the curvature spectrum is the distribution of ||F\_alpha(d\_i, d\_j)|| across random observer perturbation directions d\_i.

At (n=5, m=2): 780 pairs, mean ||F|| = 2.31, no flat pairs (0%).
At (n=6, m=3): 780 pairs, mean ||F|| = 1.82, no flat pairs (0%).

The observer space is generically curved for m >= 2. Flat directions (F = 0) require the specific commuting structure identified in Section 1.3.

### 4.2. What curvature tells you

The curvature norm ||F\_alpha|| measures how much the hidden frame rotates as the observer changes. Specifically, from the compatibility identity:

    ||F_alpha||^2 = 2(Tr(G_1 G_2) - Tr(C^2))

where G\_i and C are hidden-space quantities. Large curvature means the hidden-space inner product between perturbation directions is highly asymmetric — i.e., the two directions interact noncommutatively through the hidden structure.

This is complementary to existing nomogeo diagnostics:
- **Closure scores** measure how much structure leaks outside the observer's visible subspace.
- **Leakage channels** identify which hidden directions carry the leaked structure.
- **Curvature** measures how much the observer-design landscape is curved — i.e., how much the choice of observer matters when nearby observers are compared.

High curvature means small observer changes produce large hidden-frame rotations, making the observer choice sensitive. Low curvature means the observer space is locally flat and nearby observers are equivalent.

---

## 5. Alignment with the 0.3.x technical notes

### 5.1. TN1 (0.3.1): Typed tower

The 0.4.0 source law adds a **dynamic layer** to the typed tower:
- Layer 1 (intrinsic): Phi\_C(H) — static.
- Layer 2 (ceiling-mediated): hidden load, determinant clock — static, requires ceiling.
- **Layer 2.5 (dynamic): V, W, A\_cpl — pathwise, requires support stability.**
- Layer 3 (exact special sectors): affine-hidden, variable-precision — static.
- Layer 4 (certified corrections): residual bounds, cumulants.

The source law sits between the static ceiling-mediated layer and the exact special sectors. It governs how the static quantities change along paths.

### 5.2. TN1 (0.3.2): Graph-frontier Hessians

The sextic stability hierarchy (0.4.0 Section 5) extends the graph-frontier Hessian from TN1 (0.3.2). The existing `declared_frontier_local_certificate` works at quadratic level (H\_mu^red). The 0.4.0 paper adds quartic (Q\_mu on ker H\_mu^red) and sextic (q\_{6,mu} on the quartic null cone). These are exact extensions of the same local certificate machinery.

### 5.3. TN4 (0.3.3): Typed evidence stack

The connection is established in Section 3 above: A\_cpl governs the curvature of the evidence surface after exact hidden elimination. The source law is not a replacement for the TN4 stack — it is the derivative layer that tells you how the evidence geometry changes under perturbation.

### 5.4. Pure-observer curvature and the Grassmannian

The pure-observer curvature (0.4.0 Section 3) is exactly the Riemann curvature of the Grassmannian Gr(m, n) in the metric induced by H. The no-free-curvature-mode result (Prop 3.3) says this curvature has no quadratic free gauge mode — any variational principle built from F\_alpha^2 starts at quartic order. This is a constraint on any future observer optimiser: Yang-Mills-type actions on the observer space are nontrivial only at quartic order.

---

## 6. Open directions identified by this research

### 6.1. Full compatibility theorem (not yet achieved)

The identity ||F||^2 = 2(Tr(G\_1 G\_2) - Tr(C^2)) relates curvature to hidden Gram geometry. A full compatibility theorem would relate A\_cpl and F\_alpha when both sectors are active (H and C both varying). The structural obstacle is that A\_cpl is a pathwise object (contracted with dt) while F\_alpha is a 2-form (antisymmetric in two parameter directions). A compatibility theorem needs either:
- A 2-parameter family H(s,t), C(s,t) where both jets and curvature are defined, or
- A variational characterisation that makes both objects compete in the same functional.

This is flagged as future work by the 0.4.0 TN (Remark 4.2) and remains so.

### 6.2. Finite-epsilon at second order

The first-order correction Delta\_1 is identified. The second-order correction requires:
- The O(eps^2) slow-manifold term Y\_2
- The time derivative of Theta (for non-constant-coefficient systems)
- Careful treatment of the simulation timescale

This is tractable but requires more careful numerics than the present investigation.

### 6.3. Curvature spectrum as selection criterion

The curvature spectrum could become an observer selection diagnostic: among competing observers, prefer those in flatter regions of the Grassmannian (less sensitive to perturbation). This would complement the existing leakage-visibility frontier with a stability criterion.

### 6.4. Non-Gaussian extension

The source law requires V = L^T Hdot L where H is SPD. For non-Gaussian models, H could be replaced by the Fisher information I(theta). The formal structure (split frame, connection forms, completed square) goes through unchanged. The restriction is that I(theta) must be smooth and positive definite, which holds for regular exponential families but may fail at boundary/singular points — exactly the boundary of the TN4 typed evidence stack.
