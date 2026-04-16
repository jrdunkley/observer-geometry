# Field Equation — Well-Posed Formulation

**Date:** 16 April 2026  
**Status:** Theorem proved (modulo regularity detail)  

---

## Why Fisher-Rao regularisation is insufficient

**Claim in earlier notes:** d_FR(H, H₀)² regularises vis_rate.

**This is WRONG.** As λ_min(H) → 0:
- vis_rate grows as O(1/λ_min) = O(exp(|log λ_min|))
- d_FR² grows as O((log λ_min)²)

Exponential beats polynomial: 1/λ dominates (log λ)². The Fisher-Rao distance is not coercive against vis_rate divergence.

---

## The correct regulariser: KL divergence

**Definition.** For H, H₀ ∈ SPD(n):

    D_KL(H₀ || H) = (1/2)[Tr(H₀ H⁻¹) − log det(H₀ H⁻¹) − n]

This grows as:
- O(1/λ_min) as λ_min(H) → 0  (from the Tr term)
- O(log λ_max) as λ_max(H) → ∞  (from the −log det term)

So it matches or dominates vis_rate in all escape directions.

**The regularised objective:**

    J(H) = vis_rate(H, C, Ξ) − γ · D_KL(H₀ || H)

---

## Theorem (Well-posedness)

**Statement.** For any γ > 0 sufficiently large (γ > γ₀ where γ₀ depends on Ξ, C, H₀), the functional J is bounded above on SPD(n) and attains its supremum in the interior.

**Proof.**

*Step 1: vis_rate growth bound.*

|vis_rate(H, C, Ξ)| ≤ ||H⁻¹||₂ · (||Ξ||_nuc + ||Ξ||_nuc) = 2||Ξ||_nuc / λ_min(H)

(Using |Tr(H⁻¹Ξ)| ≤ ||Ξ||_nuc/λ_min and |Tr(R⁻¹U_h)| ≤ ||Ξ||_nuc/λ_min since λ_min(R) ≥ λ_min(H).)

Set K = 2||Ξ||_nuc.

*Step 2: D_KL growth bound.*

D_KL(H₀||H) ≥ (1/2)[λ_min(H₀)/λ_min(H) − log(λ_min(H₀)/λ_min(H)) − 1]

For λ_min(H) = ε → 0: D_KL ≥ (1/2)[c₀/ε − log(c₀/ε) − 1] ≥ c₀/(4ε) for small enough ε.

Where c₀ = λ_min(H₀) > 0.

*Step 3: Combine.*

    J(H) ≤ K/ε − γ c₀/(4ε)

For γ > 4K/c₀ = 8||Ξ||_nuc/λ_min(H₀):

    J(H) ≤ −(γ c₀/4 − K)/ε → −∞  as ε → 0

*Step 4: The large-eigenvalue direction.*

As λ_max(H) → ∞: vis_rate → 0 (since H⁻¹ → 0 in that direction). And D_KL ≥ (1/2)(−log det(H₀H⁻¹) − n) → +∞ (from the log det term picking up −log λ_max). So J → −∞ in this direction too.

*Step 5: Compactness of sublevel sets.*

For any c ∈ ℝ, the set S_c = {H ∈ SPD(n) : J(H) ≥ c} satisfies:
- λ_min(H) ≥ ε₀(c) > 0 (from Step 3)
- λ_max(H) ≤ M₀(c) < ∞ (from Step 4)

So S_c ⊂ {H ∈ SPD(n) : ε₀ I ≼ H ≼ M₀ I}, which is compact. Since J is continuous, the supremum is attained. ∎

**Remark.** The threshold γ₀ = 8||Ξ||_nuc/λ_min(H₀) depends on the perturbation and reference metric. For γ < γ₀, the functional may still be unbounded (the penalty is too weak to overcome the vis_rate divergence). This is the correct physical statement: the regularisation strength must exceed a threshold set by the perturbation intensity relative to the reference precision.

---

## The field equation

At the maximum H*, ∇_H J = 0:

    ∇_H vis_rate = γ · ∇_H D_KL(H₀||H)

**Gradient of D_KL:**

    ∇_H D_KL(H₀||H) = (1/2)(H⁻¹ − H⁻¹ H₀ H⁻¹)

**Proof:** δD_KL = (1/2)[Tr(H₀(−H⁻¹δHH⁻¹)) + Tr(H⁻¹δH)] = (1/2)Tr((H⁻¹ − H⁻¹H₀H⁻¹)δH). ∎

**The stationarity equation:**

    −H⁻¹ Ξ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = (γ/2)(H⁻¹ − H⁻¹ H₀ H⁻¹)

Multiplying left and right by H:

### **S − Ξ = (γ/2)(H − H₀)**

where S = HZ R⁻¹ U_h R⁻¹ Z^T H is the hidden stress tensor.

---

## Structural analysis of the field equation

### The balance law

    S − Ξ = (γ/2)(H − H₀)

The hidden stress minus the perturbation is proportional to the metric's deviation from reference. This is a genuine, well-posed balance equation.

### At H = H₀ (reference configuration)

    S₀ = Ξ

The reference metric is the configuration where hidden stress exactly balances the perturbation. This is a self-consistency condition: H₀ is a solution iff S(H₀, C, Ξ) = Ξ.

### Linear response (small deviations)

For H = H₀ + ε δH:

    S(H₀ + εδH) − Ξ ≈ ε · (∂S/∂H)|_{H₀} δH = (γ/2) ε δH

So: **(∂S/∂H)|_{H₀} = (γ/2) I** at the reference configuration (when H₀ is itself a solution).

### Split-frame projection

In the adapted frame M = [L|Z], using the tensor residual decomposition (Theorem 3):

    M^T(S − Ξ)M = [[-V, -B], [-B^T, 0]]

And M^T(H − H₀)M = M^THM − M^TH₀M = [[Φ, 0], [0, R]] − [[Φ₀, K₀], [K₀^T, R₀]]

where (Φ₀, R₀, K₀) are the split-frame coordinates of H₀ in the adapted basis for H.

So the block equations are:

    VV: −V = (γ/2)(Φ − Φ₀ − coupling terms)
    VH: −B = (γ/2)(0 − K₀ − ...)
    HH:  0 = (γ/2)(R − R₀ − ...)

The hidden-hidden block gives: **R = R₀** (plus coupling corrections). The hidden metric at the optimum agrees with the reference, to leading order.

The off-diagonal block: **B = (γ/2) K₀** (plus corrections). At the adapted observer (B = 0), this requires K₀ = 0 — the reference metric must also be block-diagonal in the observer basis.

The visible block: **V = −(γ/2)(Φ − Φ₀)**. The visible jet is proportional to the visible precision deviation. This is the "driving equation": the perturbation's visible impact determines how far the optimal H deviates from reference in the visible sector.

---

## Summary

### What is proved

1. **Spectral collapse lemma:** Any scalar constraint (det, trace, etc.) that allows λ_min → 0 leaves vis_rate unbounded.

2. **Fisher-Rao insufficiency:** d_FR² grows as (log λ)², vis_rate as 1/λ. The FR penalty is too weak.

3. **KL well-posedness:** D_KL(H₀||H) grows as 1/λ_min (matching vis_rate). For γ > 8||Ξ||/λ_min(H₀), the problem is well-posed with compact sublevel sets.

4. **Field equation:** S − Ξ = (γ/2)(H − H₀). Clean balance law with genuine interior solutions.

5. **Split-frame projection at B = 0:** The field equation separates into visible (V drives Φ−Φ₀), coupling (B = 0 requires K₀ = 0), and hidden (R ≈ R₀) blocks.

### Honest boundaries

- The KL regulariser introduces a reference metric H₀ and a coupling γ. These are not determined by the observer geometry alone. An additional principle is needed to fix them.
- The field equation S − Ξ = (γ/2)(H − H₀) is algebraic (not a PDE). It is observer-relative (S depends on C). It is a stationarity law, not a conservation law.
- The split-frame projection is clean at B = 0 but involves coupling corrections for general observers.
- The threshold γ₀ = 8||Ξ||/λ_min(H₀) means the regularisation must be strong enough relative to the perturbation. For very large perturbations, the penalty must be correspondingly strong.

### The next theorem to aim for

> Under KL regularisation with reference H₀, the joint (H*, C*) optimum satisfies:
> - C* is the adapted observer (B = 0)
> - H* satisfies S − Ξ = (γ/2)(H − H₀)
> - The split-frame projection gives V = −(γ/2)(Φ − Φ₀), R = R₀
> - The hidden metric is frozen at its reference value; only the visible precision deviates
