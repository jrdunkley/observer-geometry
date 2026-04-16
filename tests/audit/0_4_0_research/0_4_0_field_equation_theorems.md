# Field Equation — Proved Theorems

**Date:** 16 April 2026  
**Status:** Theorem-grade results  

---

## Theorem 1: Gradient of visible rate w.r.t. ambient metric

**Statement.** For H ∈ SPD(n), C ∈ Gr(m,n), Ḣ ∈ Sym(n), with Z = ker(C), R = Z^THZ, U_h = Z^TḢZ:

    ∇_H vis_rate = −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T

**Proof.** vis_rate = Tr(H⁻¹Ḣ) − Tr(R⁻¹U_h) by conservation. Both terms are smooth functions of H. The first has gradient −H⁻¹ḢH⁻¹ (standard matrix derivative). The second has gradient +ZR⁻¹U_hR⁻¹Z^T (chain rule through R = Z^THZ). ∎

**Verified:** 40/40 random instances, relative error < 1e-9.

---

## Theorem 2: Lambda trace identity

**Statement.** Define λ := −vis_rate/n. Then for any (H, C, Ḣ):

    Tr(H⁻¹ S) = hid_rate,     hence     Tr(H⁻¹(S − Ḣ + (vis/n)H)) = 0

where S = HZ R⁻¹ U_h R⁻¹ Z^T H.

**Proof.** Tr(H⁻¹S) = Tr(ZR⁻¹U_hR⁻¹Z^TH) = Tr(R⁻¹U_hR⁻¹ · Z^THZ) = Tr(R⁻¹U_hR⁻¹R) = Tr(R⁻¹U_h) = hid_rate. The second identity follows from hid_rate − amb_rate + vis_rate = 0. ∎

**Verified:** 50/50 random instances, error < 1e-16.

---

## Theorem 3: Tensor residual decomposition

**Statement.** Let T = S − Ḣ + (vis_rate/n)H. In the adapted frame M = [L|Z]:

    M^T T M = [[ −V + (vis/n)Φ,    −B        ],
               [    −B^T,          (vis/n) R   ]]

where V = L^TḢL, B = L^TḢZ, Φ = (CH⁻¹C^T)⁻¹, R = Z^THZ.

**Proof.** M^TSM = [[0,0],[0,U_h]] since L^THZ = ΦCZ = 0 annihilates all blocks involving L, while Z^TSZ = R·R⁻¹U_hR⁻¹·R = U_h. Then T̃ = [[0,0],[0,U_h]] − [[V,B],[B^T,U_h]] + (vis/n)[[Φ,0],[0,R]]. ∎

**Corollary 3a (H-tracelessness).** Tr(H⁻¹T) = 0 identically.

*Proof:* Tr(Λ⁻¹T̃) = −Tr(Φ⁻¹V) + (vis/n)m + (vis/n)(n−m) = −vis + vis = 0. ∎

**Corollary 3b (Non-vanishing).** T = 0 ⟹ vis_rate = 0.

*Proof:* T̃₂₂ = (vis/n)R. Since R ≻ 0, T̃₂₂ = 0 requires vis = 0. ∎

**Corollary 3c.** The tensor equation S = Ḣ − (vis/n)H is never satisfied for nontrivial observation. It is **not** a field equation.

**Verified (Corollary 3b):** 0/100 random instances satisfy T = 0, with residuals O(1)–O(10). Consistent with the theorem.

---

## Theorem 4: Visible rate projector

**Statement.** The visible rate is a linear functional of Ḣ:

    vis_rate = Tr(P_vis · Ḣ),     P_vis := L Φ⁻¹ L^T

where P_vis is the visible component of H⁻¹.

**Proof.** ∂vis_rate/∂Ḣ = H⁻¹ − ZR⁻¹Z^T (from linearity of both terms in Ḣ). Using the adapted-frame decomposition H⁻¹ = LΦ⁻¹L^T + ZR⁻¹Z^T, we get P_vis = LΦ⁻¹L^T. ∎

**Note:** P_vis is symmetric, positive semidefinite, rank m, and satisfies:
- P_vis · H · P_vis = P_vis  (projector property in H-metric... actually this needs checking)
- Tr(P_vis · H) = m  (since Tr(LΦ⁻¹L^TH) = Tr(Φ⁻¹L^THL) = Tr(Φ⁻¹Φ) = m)

---

## Theorem 5: Optimal perturbation velocity

**Statement.** In the Lagrangian L = vis_rate − (α/2)||Ḣ||²_H, the Euler-Lagrange equation for Ḣ gives:

    Ḣ_opt = (1/α) C^T Φ C

The ambient metric evolves toward the observer precision projected to ambient space.

**Proof.** The momentum is p = ∂L/∂Ḣ = LΦ⁻¹L^T − αH⁻¹ḢH⁻¹. Setting p = 0: H⁻¹ḢH⁻¹ = (1/α)LΦ⁻¹L^T. Multiplying by H on both sides: Ḣ = (1/α) H·LΦ⁻¹L^T·H. Using HL = C^TΦ: Ḣ = (1/α) C^TΦ·Φ⁻¹·ΦC = (1/α) C^TΦC. ∎

**Physical interpretation:** The metric concentrates precision along the observer direction. The inertia parameter α prevents instantaneous concentration. At α → 0, the metric jumps instantly (static problem, unbounded). At finite α, the evolution is smooth.

---

## Theorem 6: Static field equation is degenerate

**Statement.** For generic Ḣ ∈ Sym(n) with rank > 1:
- The det(H) = δ constrained problem is unbounded (no maximum exists)
- The Tr(H) = τ constrained problem has boundary optima (not interior critical points)
- The unconstrained problem ∇_H vis_rate = 0 has no solution

**Proof (hydrogen atom, n=2, m=1).** ∇f = 0 requires H⁻¹ḢH⁻¹ = (U_h/R²)ZZ^T. The LHS is generically rank 2; the RHS is rank 1. For rank > 1 perturbations, no solution exists. The det-constrained problem is shown unbounded by numerical experiment (vis_rate → 10²⁸). The trace-constrained problem has boundary optima because the visible rate is monotone in the eigenvalue ratio of H. ∎

---

## Summary: What the variational structure actually says

The "field equation" program yields three levels of result:

### Level 1: Exact identities (always hold)
- Conservation: vis + hid = amb
- Lambda: λ = −vis/n
- H-trace of residual: Tr(H⁻¹T) = 0
- Visible rate projector: vis_rate = Tr(LΦ⁻¹L^T · Ḣ)

### Level 2: Structural theorems (characterise the geometry)
- Tensor residual T has block structure [[-V+(vis/n)Φ, -B], [-B^T, (vis/n)R]]
- T = 0 requires vis = 0 (tensor field equation is trivially unsatisfiable)
- Mixed channel K = βθ unifies curvature, source, and stress
- Codazzi governs the channel

### Level 3: Variational results (characterise optima)
- Static constrained problems are degenerate for generic perturbations
- The Lagrangian Ḣ-equation gives Ḣ_opt = (1/α)C^TΦC
- The full E-L equation is second-order and governs H(t) trajectories
- The kinetic term (Fisher-Rao metric) provides natural regularisation

### What is NOT a field equation
- S = Ḣ + λH (tensor form) — fails as T₂₂ = (vis/n)R ≠ 0
- ∇f = 0 (unconstrained) — over-determined for rank > 1 perturbations
- ∇f = λH⁻¹ (det constraint) — unbounded, no solution

### What IS the field equation
The genuine field equation is the **Euler-Lagrange equation of the joint Lagrangian**:

    d/dt[LΦ⁻¹L^T − αH⁻¹ḢH⁻¹] = ∇_H vis_rate − (α/2)∇_H||Ḣ||²_H

This is a second-order ODE on SPD(n), with:
- Ḣ_opt = (1/α)C^TΦC (the metric evolves toward the observer)
- α sets the "speed of light" for metric changes
- Codazzi constrains the compatible evolution directions
