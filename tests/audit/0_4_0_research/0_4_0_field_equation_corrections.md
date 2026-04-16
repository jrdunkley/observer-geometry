# Field Equation — Corrections and Revised Analysis

**Date:** 16 April 2026  
**Prompted by:** Director review  

---

## Errors corrected

### 1. Tr(H) = τ does NOT fix the static problem

**My claim:** "Tr(H) = τ is well-posed because the feasible set is compact."

**Correction:** The set {H ≻ 0 : Tr(H) = τ} is NOT compact. Its closure includes the boundary of the PSD cone (matrices with zero eigenvalues), and vis_rate diverges there.

**Counterexample (Director):** n=2, m=1, C = e₁, Z = e₂.

    H = diag(ε, τ−ε),  Ξ = diag(h, f)

    vis_rate = h/ε → +∞  as ε ↓ 0 (for h > 0)

So Tr(H) = τ permits spectral collapse and vis_rate is unbounded. The same mechanism kills det(H) = δ.

### 2. No scalar constraint suffices

**Lemma (Spectral Collapse).** Any admissible class that allows λ_min(H) → 0 while keeping the observer aligned with a positive visible component of Ξ makes vis_rate unbounded above.

**Proof.** Let C be aligned with the eigenspace where H collapses (λ_min → 0). The visible precision Φ = (CH⁻¹C^T)⁻¹ grows as λ_min, and vis_rate = Tr(Φ⁻¹V) where V = L^TΞL also involves H⁻¹ terms. In the aligned case (C = eigenvector of H with eigenvalue ε), vis_rate = Ξ_vis/ε → ∞. ∎

**Consequence:** Fixed det, fixed trace, or any scalar budget that does not control anisotropy fails. The real fix requires a coercive regulariser that penalises spectral collapse.

### 3. The displayed ODE is not a derived field equation

**My claim:** "The field equation (final form) is α·d/dt[H⁻¹ḢH⁻¹] = ..."

**Correction:** This is a candidate E-L equation from a chosen action, not forced by the geometry. The geometry gives conservation laws and gradient formulas. A field equation requires an additional variational principle.

### 4. Four problems must be separated

| Problem | Fixed | Optimised | Type |
|---------|-------|-----------|------|
| Pointwise | (C, Ξ) | H | Static algebra |
| Joint static | Ξ | (H, C) | Static algebra |
| Path | — | H(t), C(t) | Calculus of variations |
| Compatibility | — | — | Codazzi, structure eqs |

I was drifting between these. Codazzi constrains Problem 4, not Problem 1.

### 5. Codazzi is pathwise, not pointwise

Codazzi (Dβ = Dθ = 0) governs how the split-frame data varies along a path, bundle, or field. It does not regularise a single-point optimisation in H.

### 6. Notation fix

Rename the externally prescribed perturbation from Ḣ to **Ξ** in the static theory. Reserve Ḣ = dH/dt for the dynamic theory.

---

## What remains valid

| Result | Status |
|--------|--------|
| Gradient: ∇_H vis_rate = −H⁻¹ΞH⁻¹ + ZR⁻¹U_hR⁻¹Z^T | PROVED |
| Lambda trace identity: λ = −vis_rate/n | PROVED |
| Tensor residual: T̃ = [[-V+(vis/n)Φ, -B], [-B^T, (vis/n)R]] | PROVED |
| T = 0 requires vis = 0 | PROVED |
| Visible rate projector: P_vis = LΦ⁻¹L^T | PROVED |
| Channel unification: K = βθ, three projections | VERIFIED |
| det(H) = δ is ill-posed | PROVED |
| Spectral collapse kills all scalar constraints | PROVED (lemma above) |

## What is retracted

| Claim | Status |
|-------|--------|
| "Tr(H) = τ works structurally" | RETRACTED — same spectral collapse |
| "Fisher-Rao is the only formulation" | RETRACTED — it's the most natural candidate, not unique |
| "The field equation (final form) is ..." | RETRACTED — candidate, not derived |
| "Hydrogen atom" labelling | RETRACTED — use "2×2 split model" |

---

## The real obstruction (Director's formulation)

**Spectral collapse makes the static visible-rate maximisation ill-posed unless the geometry includes a coercive penalty or a genuinely dynamical cost.**

This is a strong structural result. It tells us:

1. The problem CANNOT be solved by constraint selection
2. The problem REQUIRES either a reference-metric regulariser or a dynamical action
3. The natural candidate is a Fisher-Rao / log-det barrier relative to a reference H₀
