# Field Equation — Coercive Regularisation

**Date:** 16 April 2026  
**Status:** Theorem + open derivation  

---

## Notation (corrected)

- **Ξ** ∈ Sym(n): the externally prescribed perturbation (static theory)
- **Ḣ** = dH/dt: time derivative of H along a path (dynamic theory)
- These are distinct objects; the static theory fixes Ξ, the dynamic theory has Ḣ as a velocity.

---

## Lemma (Spectral Collapse)

**Statement.** Let f(H) = vis_rate(H, C, Ξ) for fixed (C, Ξ) with Ξ having a positive component in the visible sector. For any admissible class S ⊂ SPD(n) that contains a sequence Hₖ with λ_min(Hₖ) → 0 and C aligned with the collapsing eigenspace:

    sup_{H ∈ S} f(H) = +∞

**Proof.** Align C with the eigenvector of Hₖ corresponding to λ_min = εₖ. In this basis, H ≈ diag(εₖ, ...) and vis_rate ≈ Ξ_vis/εₖ → ∞. ∎

**Corollary.** Fixed det(H) = δ, fixed Tr(H) = τ, or any scalar constraint that does not control λ_min/λ_max fails. The static optimisation is ill-posed under purely scalar budgets.

---

## The coercive regularisation

### The objective

    J(H) = vis_rate(H, C, Ξ) − (γ/2) d_FR(H, H₀)²

where:
- H₀ ∈ SPD(n) is a reference metric (given)
- d_FR is the Fisher-Rao (affine-invariant) distance:

    d_FR(H, H₀)² = ||log(H₀^{-1/2} H H₀^{-1/2})||²_F = Σᵢ (log λᵢ)²

  with λᵢ the eigenvalues of H₀⁻¹H.
- γ > 0 is the regularisation strength

### Why this is well-posed

**Proposition (Well-posedness).** For any γ > 0, J(H) is bounded above on SPD(n) and attains its maximum in the interior.

**Proof sketch.**

*Boundedness above:* As any eigenvalue λᵢ of H₀⁻¹H approaches 0 or ∞:
- vis_rate grows at most as O(1/λ_min) = O(exp(|log λ_min|))
- d_FR² grows as (log λ_min)² 

For any polynomial p and any γ > 0: p(e^x) − (γ/2)x² → −∞ as |x| → ∞. So the penalty dominates and J → −∞ at the PSD boundary and at infinity.

*Attainment:* The sublevel sets {H : J(H) ≥ c} are bounded and bounded away from ∂SPD(n) for large enough c. Since J is continuous on SPD(n), the maximum is attained in the interior. ∎

**Remark.** The proof sketch above is incomplete. The vis_rate growth rate depends on the alignment of C with the collapsing eigenspace, so the bound vis_rate = O(1/λ_min) needs to be made uniform over the constraint surface. A complete proof would use the fact that d_FR controls ALL eigenvalue ratios simultaneously, not just λ_min.

### The stationarity equation

At the maximum H*, ∇_H J = 0:

    −H⁻¹ Ξ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = γ · G_FR(H, H₀)

where G_FR(H, H₀) = ∇_H d_FR(H, H₀)² is the Euclidean gradient of the squared Fisher-Rao distance.

**In the commuting case** (H and H₀ simultaneously diagonalisable):

    G_FR = 2 H⁻¹ log(H₀⁻¹H) H⁻¹

So the stationarity becomes:

    −H⁻¹ Ξ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = 2γ H⁻¹ log(H₀⁻¹H) H⁻¹

Multiplying left and right by H:

    −Ξ + H Z R⁻¹ U_h R⁻¹ Z^T H = 2γ log(H₀⁻¹H)

i.e., **S − Ξ = 2γ log(H₀⁻¹H)** where S is the hidden stress tensor.

### Structural interpretation

    S − Ξ = 2γ log(H₀⁻¹H)

The hidden stress minus the perturbation equals a logarithmic deviation from the reference metric. This has the structure:

    (geometry) − (source) = (regulariser)

At H = H₀ (no deviation from reference): S₀ − Ξ = 0, i.e., the hidden stress must exactly balance the perturbation. This is a genuine balance equation.

For small deviations H = H₀ + εδH: log(H₀⁻¹H) ≈ H₀⁻¹δH, so:

    S − Ξ ≈ 2γ H₀⁻¹ δH

which is a linear response equation. The deviation from reference is proportional to the stress-perturbation imbalance.

---

## Split-frame projection of the stationarity equation

In the adapted frame M = [L|Z], using Theorem 3 (tensor residual):

The LHS (S − Ξ) in adapted coordinates:

    M^T(S − Ξ)M = [[0,0],[0,U_h]] − [[V,B],[B^T,U_h]]
                 = [[-V, -B], [-B^T, 0]]

The RHS (2γ log(H₀⁻¹H)) in adapted coordinates depends on the relationship between H₀ and the adapted frame. For the special case H₀ = H (γ → 0 limit):

    log(H₀⁻¹H) = 0, so S = Ξ.

But we know S ≠ Ξ generically (Theorem 3). So H₀ = H is NOT a solution of the stationarity equation unless S = Ξ happens to hold (which requires the adapted-frame condition [[-V,-B],[-B^T,0]] = 0, i.e., V = 0 and B = 0).

### The B = 0 (adapted observer) sector

At B = 0:

    M^T(S − Ξ)M = [[-V, 0], [0, 0]]

So the off-diagonal and hidden-hidden blocks vanish. The stationarity equation reduces to:

    Visible block: −V = 2γ (log(H₀⁻¹H))_vis
    Hidden block:  0 = 2γ (log(H₀⁻¹H))_hid
    Coupling:      0 = 2γ (log(H₀⁻¹H))_coupling

The hidden block equation says log(H₀⁻¹H) has zero hidden-hidden component, meaning H and H₀ agree in the hidden sector (up to the adapted frame projection).

The visible block equation says V = −2γ (log deviation)_vis, which determines the visible component of H's deviation from H₀ in terms of the visible jet V.

**This is a genuine block-diagonal balance law at B = 0.** It recovers the visible source term (V) as the driver and the reference-metric deviation as the response.

---

## Four problems, clearly separated

### Problem 1: Pointwise optimisation (fixed C, Ξ, vary H)

**Result:** Ill-posed under scalar constraints. Well-posed under coercive (FR) regularisation. Interior critical points exist. Stationarity equation: ∇vis_rate = γ · G_FR.

### Problem 2: Joint static optimisation (fixed Ξ, vary H and C)

**Result:** The observer equation (vary C) gives the adapted observer (B = 0). Combined with Problem 1, the joint optimum has B = 0 and the stationarity equation simplifies to a block-diagonal law.

### Problem 3: Path optimisation (vary H(t), C(t))

**Status:** The Lagrangian formulation with Fisher-Rao kinetic term is a candidate but not yet derived. The E-L equation needs proper covariant derivation on SPD(n).

### Problem 4: Geometric compatibility (Codazzi along paths)

**Status:** Codazzi constrains pathwise evolution of the split frame. It does not regularise the static problem. Its role in Problem 3 needs investigation.

---

## Next steps

1. **Complete the well-posedness proof** for the coercive regularisation (make the vis_rate growth bound rigorous)
2. **Compute G_FR in the non-commuting case** (the general Fréchet derivative of log on SPD(n))
3. **Derive the split-frame projection** of the stationarity equation without the commuting assumption
4. **Verify numerically** on the 2×2 split model (lightweight: check that J has an interior maximum and the stationarity equation holds there)
5. **Connect to Problem 3** — does the coercive static regulariser emerge as the zero-velocity limit of a dynamic action?
