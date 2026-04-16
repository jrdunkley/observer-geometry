# T5: Non-Gaussian Extension Boundary of the Source Law

**Date:** 15 April 2026
**Status:** Theoretical analysis (no code required)

---

## 1. What the source law requires

The 0.4.0 source law operates on a smooth field H: U -> SPD(n) of positive definite quadratic forms, with a rank-m observation map C: U -> R^{m x n}. The key objects are:

- Split frame M = [L Z] with M^T H M = diag(Phi, R)
- Visible first jet V = L^T Hdot L (requires H differentiable)
- Covariant visible second jet W = D_t V (requires H twice differentiable)
- Support stability: rank(V) constant along the path (requires V full rank on its support)
- V\_S > 0 on the active support (required for A\_cpl)

**The source law does NOT require:**
- H to be a covariance or precision matrix
- H to come from a Gaussian distribution
- H to be related to a likelihood function
- n to be the dimension of any particular space

**The source law DOES require:**
- H(t) smooth and SPD for all t in the path
- C smooth and rank m
- V\_S > 0 on the active support

---

## 2. Non-Gaussian models via Fisher information

For a parametric statistical model p(x | theta), the Fisher information matrix is:

    I(theta) = E[score * score^T] = -E[d^2 log p / dtheta dtheta^T]

This is SPD whenever the model is identifiable and the score exists. If we set H(theta) = I(theta) and let C be an observer that extracts visible parameter directions, the source law applies directly.

**What changes in the non-Gaussian case:**
- I(theta) is no longer the inverse of a covariance matrix
- The hidden defect Q\_hat = B R^{-1} B^T no longer has a simple interpretation as "hidden variable precision"
- The evidence connection (TN4) involves the Hessian of the log-likelihood at the MLE, which equals I(theta\_hat) only in the regular interior sector

**What does NOT change:**
- The split-frame identities (purely algebraic, hold for any SPD field)
- The connection forms (alpha, beta, theta, omega)
- The curvature F\_alpha = -beta wedge theta
- The completed square for W
- A\_cpl decomposition

The source law is mathematically valid for ANY smooth SPD field. The non-Gaussian restriction is not in the source law itself but in the **interpretation** of A\_cpl as an evidence curvature.

---

## 3. Where the source law breaks down

The source law fails at precisely the same boundary as the TN4 regular interior sector:

### 3.1. Fisher information not SPD (identifiability failure)

When the model is not identifiable at theta\_0, the Fisher information I(theta\_0) is only PSD, not SPD. The split frame is not well-defined (Phi = (C I^{-1} C^T)^{-1} requires I invertible).

**This is the quotient boundary.** TN4 handles this via exact slice reduction (Theorem 6.2): removable symmetry directions are quotiented out, and the source law applies on the transverse slice.

### 3.2. Fisher information not smooth (boundary/singular points)

At boundary points of the parameter space (e.g., a variance parameter approaching zero), I(theta) may fail to be smooth or even continuous. The derivative dI/dtheta may not exist.

**This is the cone boundary.** TN4 handles this via the tangent-cone Gaussian mass formula (Theorem 5.2). The source law does not apply at the boundary itself, but it applies along interior paths approaching the boundary.

The branch-restiffening analysis (Goal E) showed exactly this: as Phi -> 0 (approaching the boundary), A\_cpl diverges as -1/(2t), and the evidence second derivative develops a cusp. The source law DETECTS the approach to the boundary through A\_cpl divergence, even though it cannot operate AT the boundary.

### 3.3. Support instability (rank transitions)

When the rank of V = L^T Hdot L changes along the path, the support-stable stratum assumption fails. This happens when:
- An eigenvalue of V crosses zero (a visible direction "turns on" or "turns off")
- The observation map C loses rank (an observer channel dies)

**This is the support-event boundary.** The existing nomogeo observation-field layer handles support-stratum transport for the static case. The source law operates on each support-stable stratum separately.

### 3.4. The genuinely unresolved boundary

The only case where the source law provides no information is when H(theta) is not twice differentiable at theta\_0. This can happen for:
- Models with non-analytic Fisher information (e.g., mixture models at the mixing proportion boundary)
- Models where the Fisher information depends on the data (post-selection inference)
- Nonparametric models without a smooth parameter space

These are exactly the cases where the TN4 reduced kernel datum is needed: the source law's domain of validity (smooth SPD field) coincides with the TN4 regular and cone sectors, and the genuinely unresolved kernel sector lies outside both.

---

## 4. The extension map

| TN4 sector | Fisher I(theta) | Source law | Evidence curvature |
|------------|----------------|------------|-------------------|
| Regular interior | SPD, smooth | Fully valid | f'' from A\_cpl |
| Regular cone | SPD interior, PSD boundary | Valid on interior, detects boundary via A\_cpl divergence | Cusp at boundary |
| Regular quotient | PSD (kernel from symmetry) | Valid on transverse slice | Slice-reduced f'' |
| Certified singular (weighted-homogeneous) | PSD, higher-order kernel | Source law breaks down | Need kernel germ |
| Certified singular (branch-restiffening) | PSD, quadratic kernel | Source law divergence detects approach | log(n) correction |
| Unresolved kernel | Non-smooth or degenerate | Not applicable | Need reduced kernel datum |

---

## 5. The practical boundary

For practitioners: the source law applies whenever you can compute dI/dtheta and d^2I/dtheta^2 and the Fisher information is invertible. This covers:

- **All regular exponential families** (Fisher information is smooth and SPD on the interior)
- **All Gaussian models** (Fisher = function of design and variance parameters, smooth on the interior)
- **All smooth nonlinear regression models** (Fisher from the Jacobian, smooth away from singularities)
- **All GLMs** (Fisher = W^T X^T X W where W is the weight matrix, smooth at interior points)

It does NOT cover:
- Mixture models at component-death boundaries
- Singular learning machines (RLCT theory)
- Nonparametric models
- Post-selection inference

This boundary exactly matches the boundary of BIC-type asymptotics, which is reassuring: the source law extends the determinant penalty (the TN4 regular sector), and its boundary of validity is the same as the boundary of the penalty itself.

---

## 6. Conclusion

The source law has a clean non-Gaussian extension: replace H with the Fisher information I(theta). The mathematical content (split frame, connection, curvature, A\_cpl) is unchanged. The interpretive content (evidence curvature, model comparison) is valid whenever the TN4 regular sector applies. The boundary of the source law coincides exactly with the boundary of the regular sector in the TN4 typed evidence stack.

The source law does not need a separate "non-Gaussian theory." It needs the Fisher information to be smooth and SPD, which is the same condition that the TN4 evidence geometry needs to work. The two frameworks are structurally compatible because they share the same domain of validity.
