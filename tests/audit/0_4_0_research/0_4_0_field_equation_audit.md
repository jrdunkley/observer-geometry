# Field Equation — Theorem Hardening Audit

**Date:** 16 April 2026  
**Type:** Audit (not expansion)  

---

## SECTION A. THEOREM AUDIT

| Claim | Status | Reason |
|-------|--------|--------|
| Conservation vis + hid = amb | PROVED | Algebraic identity. Verified 50/50. |
| Gradient ∇_H vis = −H⁻¹ΞH⁻¹ + ZR⁻¹U_hR⁻¹Z^T | PROVED | Standard matrix calculus on vis = Tr(H⁻¹Ξ) − Tr(R⁻¹U_h). Verified 40/40. |
| Lambda trace identity λ = −vis/n | PROVED | Cyclic trace Tr(H⁻¹S) = Tr(R⁻¹U_h) = hid. Verified 50/50. |
| Tensor residual decomposition (Thm 3) | PROVED | Block computation, L^THZ = 0 annihilates S's vis-vis and vis-hid blocks. |
| T = 0 ⟹ vis = 0 | PROVED | T̃₂₂ = (vis/n)R ≻ 0. |
| P_vis = LΦ⁻¹L^T | PROVED | H⁻¹ = LΦ⁻¹L^T + ZR⁻¹Z^T, subtract hidden part. |
| Spectral collapse lemma | PROVED | vis = Ξ_vis/ε → ∞ when aligned. Kills det, Tr, all scalar constraints. |
| FR distance insufficient | PROVED | (log λ)² vs 1/λ. |
| Channel unification K = βθ | PROVED | F antisymm 60/60, Q̂ PSD 60/60, Q̂ = BR⁻¹B^T 60/60. |
| ∇D_KL = (1/2)(H⁻¹ − H⁻¹H₀H⁻¹) | PROVED | Direct variational computation. |
| KL well-posedness, fixed C | PROVED (sketch) | vis ≤ K/λ_min, D_KL ≥ c₀/(2λ_min). For γ > 2K/c₀ the penalty dominates in all escape directions. Main argument correct; constant tracking informal. |
| S − Ξ = (γ/2)(H − H₀), fixed C | PROVED | Stationarity of J, multiply by H on both sides. Holds at any interior critical point of J. No closed form for H* claimed. |
| Adapted observer at joint (H*,C*) | PROVED (corrected) | ∂J/∂C = ∂vis/∂C since D_KL is C-independent. Gives B(H*,C*) = 0. The adaptedness is w.r.t. H*, not H₀. The notes' claim "adapted for (H₀,Ξ)" is **wrong**. |
| K₀* = 0 at joint optimum | PROVED | At B(H*,C*) = 0, field equation VH block gives 0 = −(γ/2)K₀* where K₀* = L(H*)^TH₀Z. Forces K₀* = 0. Derived, not assumed. |
| R* = Z^TH₀Z at joint optimum | PROVED | HH block gives 0 = (γ/2)(R* − R₀*) where R₀* = Z^TH₀Z (Z depends on C*, not H). |
| Closed-form Φ* = Φ₀ − (2/γ)V | CONDITIONAL | VV block gives Φ* = Φ₀* − (2/γ)V* where Φ₀* = L(H*)^TH₀L(H*) and V* = L(H*)^TΞL(H*). Both depend on L(H*) which depends on H*. This is an implicit fixed-point equation, not a closed form, except in the compatibility sector where L(H*) = L(H₀). |
| L(H*) = L(H₀) in compatibility sector | PROVED | When H* = M₀⁻ᵀ diag(Φ*,R₀)M₀⁻¹, compute (H*)⁻¹C^T = L₀Φ*⁻¹, hence L(H*) = L₀. |
| Rank-m deviation formula | CONDITIONAL | Requires compatibility sector (L(H*) = L(H₀)). |
| "First equation that determines H" | CONDITIONAL | Determines H given external (H₀, γ). These are not geometric. |
| Thm 5 (Lagrangian Ḣ_opt) | FALSE-OLD | Problem 3, not Problem 1. Lagrangian is a choice. Well-posedness unaddressed. |
| Thm 6 (Tr(H) boundary optima) | FALSE-OLD | Tr(H)=τ is unbounded. Retracted. |
| FR well-posedness | FALSE-OLD | Retracted. |

---

## SECTION B. CLEAN THEOREM STATEMENT

### Theorem 1 (Fixed-observer KL-regularised stationarity)

**Hypotheses.** H₀ ∈ SPD(n), C ∈ Gr(m,n), Ξ ∈ Sym(n), γ > 2K/c₀ where K = 2||Ξ||_nuc and c₀ = λ_min(H₀). Z = ker(C).

**Define** J(H) = vis_rate(H,C,Ξ) − γ D_KL(H₀||H) on SPD(n).

**Conclusion.** J is bounded above, attains its supremum at an interior H* ∈ SPD(n), and H* satisfies:

    S(H*) − Ξ = (γ/2)(H* − H₀)

where S(H) = HZ(Z^THZ)⁻¹(Z^TΞZ)(Z^THZ)⁻¹Z^TH.

**Proof.** (1) |vis_rate| ≤ K/λ_min(H). (2) D_KL ≥ c₀/(2λ_min(H)) as λ_min → 0 and D_KL → ∞ as λ_max → ∞ (log det term). (3) J → −∞ in all escape directions for γ > 2K/c₀. (4) Sublevel sets compact. (5) J continuous ⟹ sup attained. (6) Interior maximum ⟹ ∇J = 0. (7) Multiply stationarity by H on both sides. ∎

**Remark.** No claim about C. No closed form for H*. The maximiser is implicitly defined by the stationarity equation, which is a nonlinear system in H.

---

### Proposition 2 (Joint optimum block structure)

**Hypotheses.** Same as Theorem 1, but now optimise J over (H, C) ∈ SPD(n) × Gr(m,n). Assume the joint supremum is attained at an interior point (H*, C*) and that H* is a nondegenerate critical point of J(·, C*).

**Conclusion.** Let Z = ker(C*), M* = [L(H*)|Z], and define:

    Φ₀* = L(H*)^T H₀ L(H*),   V* = L(H*)^T Ξ L(H*),   R₀* = Z^T H₀ Z

Then:

(a) B(H*,C*) = L(H*)^T Ξ Z = 0 (adapted observer for H*).

(b) K₀* = L(H*)^T H₀ Z = 0 (reference block-diagonal in H*'s adapted frame).

(c) R* = R₀* (hidden metric of H* equals hidden metric of H₀, both through Z).

(d) Φ* = Φ₀* − (2/γ)V* (implicit: LHS and RHS both depend on H* through L(H*)).

**Proof.** (a) ∂J/∂C = ∂vis/∂C since D_KL is C-independent. Standard adapted-observer theory gives B = 0. (b,c,d) Project the field equation S − Ξ = (γ/2)(H* − H₀) through the adapted frame M*. At B = 0: M*^T(S−Ξ)M* = [[-V*,0],[0,0]] and M*^T(H*−H₀)M* = [[Φ*−Φ₀*, −K₀*],[−K₀*^T, R*−R₀*]]. Equate blocks. ∎

**Remark.** Part (d) is a fixed-point condition, not a closed-form solution. The equation Φ* = Φ₀* − (2/γ)V* involves L(H*) which depends on Φ* through the reconstruction H* = M*⁻ᵀ diag(Φ*, R₀*) M*⁻¹.

---

### Proposition 3 (Compatibility-sector closed form)

**Additional hypotheses beyond Proposition 2:**

(H-compat) L(H₀)^T Ξ Z = 0 (perturbation block-diagonal in H₀'s adapted frame).

**Conclusion.** L(H*) = L(H₀), and the fixed point in Proposition 2(d) collapses to the explicit formula:

    Φ* = Φ₀ − (2/γ)V

where Φ₀ = L(H₀)^T H₀ L(H₀) and V = L(H₀)^T Ξ L(H₀), provided Φ₀ − (2/γ)V ≻ 0.

**Proof.** Under (H-compat), construct H* = M₀⁻ᵀ diag(Φ₀−(2/γ)V, R₀) M₀⁻¹ where M₀ = [L(H₀)|Z]. Direct computation: (H*)⁻¹C^T = M₀ diag((Φ₀−(2/γ)V)⁻¹, R₀⁻¹) M₀^TC^T = L₀(Φ₀−(2/γ)V)⁻¹, hence L(H*) = L₀. Therefore Φ₀* = Φ₀, V* = V, and the fixed point is solved. Verified 63/63 numerically. ∎

---

## SECTION C. GAP ANALYSIS

### Gap 1: First variation ∂vis_rate/∂C

The claim that ∂vis/∂C = 0 ⟹ B(H,C) = 0 is standard in the nomogeo kernel (closure_adapted_observer). But the current audit files do not contain a self-contained derivation. To make Proposition 2(a) fully proved rather than citing the kernel, an explicit computation of the C-variation of vis_rate = Tr(H⁻¹Ξ) − Tr(R⁻¹U_h) through the Grassmannian tangent space is needed.

**Status:** Believed correct, but a reference or derivation should be included.

### Gap 2: Existence and uniqueness of the joint optimum

Proposition 2 assumes the joint supremum is attained. The Grassmannian is compact, so the C-supremum at fixed H exists. By Theorem 1, the H-supremum at fixed C exists. But the joint supremum requires a separate argument (e.g., J(H,C) is upper semicontinuous on the product, and the joint sublevel sets are compact). This is plausible but not proved.

**Status:** Open. Likely closable with a compactness argument.

### Gap 3: Uniqueness of the fixed point in Proposition 2(d)

The equation Φ* = Φ₀* − (2/γ)V* with Φ₀* and V* depending on Φ* through L(H*) is a nonlinear fixed-point equation. Existence is guaranteed by Theorem 1 + Proposition 2. Uniqueness is open. Multiple solutions would correspond to multiple local maxima of J.

**Status:** Open.

### Gap 4: Generic non-compatibility

For generic Ξ (with B(H₀, C₀, Ξ) ≠ 0 at the initial observer), the compatibility sector is not entered. The joint optimum still exists (Gap 2 aside), and Proposition 2 still gives the block structure. But the explicit formula of Proposition 3 does not apply. The maximiser H* is implicitly defined.

**Question:** For generic Ξ, is the implicit fixed point of Proposition 2(d) well-behaved (smooth in the data, computable by iteration, etc.)? Or can it exhibit pathologies (bifurcation, non-smoothness)?

**Status:** Open.

### Gap 5: SPD condition on Φ₀ − (2/γ)V

Proposition 3 requires Φ₀ − (2/γ)V ≻ 0. For the generic Proposition 2, the analogous condition is Φ₀* − (2/γ)V* ≻ 0 at the fixed point. This is guaranteed by the well-posedness of J (the maximiser has H* ≻ 0, hence Φ* ≻ 0). But the explicit threshold γ₀ ensuring this may differ from the well-posedness threshold γ > 2K/c₀.

**Status:** Partially addressed. The well-posedness threshold guarantees H* ≻ 0, which implies Φ* ≻ 0. But the explicit relationship between γ₀ and the SPD condition on the closed-form formula needs clarification.

### Gap 6: Maximiser vs. stationary point

Theorem 1 proves the field equation at the maximiser. The maximiser is a stationary point but J may have other stationary points (saddle points, local minima). The field equation holds at all of them. No second-order condition (Hessian check) is provided to confirm the maximiser is a nondegenerate maximum.

**Status:** Open. Not required for the theorem but relevant for numerical methods.

---

## SECTION D. VERIFICATION REWRITE PLAN

### Test 1: Fixed-C generic optimiser

**Tests:** Theorem 1.

**Method:** For random (H₀, C, Ξ) with no B=0 projection, numerically optimise J(H) over SPD(n) using Cholesky parametrisation + Nelder-Mead (or L-BFGS). At the numerical maximiser H*, check:

- ||∇J(H*)|| < tol (stationarity)
- ||S(H*) − Ξ − (γ/2)(H* − H₀)|| < tol (field equation)
- H* ≻ 0 (interiority)
- J(H*) > J(H₀) (improvement)

**Label:** "Testing Theorem 1 (fixed-C stationarity)."

**Dimensions:** n ∈ {3,4}, m ∈ {1,2}, 10 restarts, γ = 50. Should be lightweight since the KL penalty makes the landscape well-behaved.

### Test 2: Compatibility-sector closed form

**Tests:** Proposition 3.

**Method:** Exactly the current `0_4_0_field_equation_verify.py`. Random (H₀, C), then construct adapted Ξ via `construct_adapted_Xi`. Evaluate the closed-form H* and check field equation + stationarity.

**Label:** "Testing Proposition 3 (compatibility sector)." Each check must state it is conditional on (H-compat).

**No changes needed** to the existing script except labelling.

### Test 3: Joint (H,C) optimiser

**Tests:** Proposition 2 (open claim).

**Method:** For random (H₀, Ξ) with generic Ξ (NOT projected to B=0), optimise J(H,C) over SPD(n) × Gr(m,n). At the numerical maximiser (H*, C*), check:

- B(H*, C*) ≈ 0 (adapted observer — Proposition 2a)
- K₀* = L(H*)^T H₀ Z ≈ 0 (reference compatibility — Proposition 2b)
- R* ≈ R₀* (hidden frozen — Proposition 2c)
- Φ* ≈ Φ₀* − (2/γ)V* (implicit equation — Proposition 2d)
- Whether L(H*) ≈ L(H₀) or not (compatibility sector vs. generic)

**Label:** "Testing Proposition 2 (joint optimum block structure) — exploratory."

**Dimensions:** n ∈ {3,4}, m = 1, γ = 50. Joint optimisation over Cholesky(H) × Stiefel(C). Heavier than Tests 1-2 but feasible for small n.

### Test 4: Counterexample search

**Tests:** Whether the compatibility sector is generic or special.

**Method:** For generic Ξ, compare:

- L(H*) vs L(H₀)
- ||Ξ − construct_adapted_Xi(H*, C*, Ξ)|| (how far Ξ is from the compatibility sector of H*)

Report statistics on whether generic Ξ approximately satisfies the compatibility hypothesis at the joint optimum.

**Label:** "Probing Gap 4 (generic non-compatibility)."
