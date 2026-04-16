# Field Equation Investigation — Findings

**Date:** 16 April 2026  
**Code:** `0_4_0_field_equation.py`, `0_4_0_field_equation_part6.py`  
**Total checks:** 300+ across 7 verification suites  

---

## 1. The variational problem

Given observer C ∈ Gr(m,n) and perturbation Ḣ ∈ Sym(n), extremise:

    f(H) = vis_rate(H, C, Ḣ) = Tr(H⁻¹Ḣ) − Tr(R⁻¹U_h)

where R = Z^THZ, U_h = Z^TḢZ, Z = ker(C).

### Gradient (PROVED, 40/40 at 1e-10)

    ∇_H f = −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T

Verified against corrected finite differences (factor-of-2 fix for Sym(n) off-diagonal entries).

### Hidden stress tensor

    S := H Z R⁻¹ U_h R⁻¹ Z^T H

This is a symmetric ambient tensor of rank ≤ n−m, built from the hidden sector.

### Lambda trace identity (PROVED, 50/50 at 1e-17)

From S = Ḣ + λH, taking Tr(H⁻¹·):

    Tr(H⁻¹S) = Tr(Z R⁻¹ U_h R⁻¹ Z^T H)
              = Tr(R⁻¹ U_h R⁻¹ · Z^THZ)    [cyclic]
              = Tr(R⁻¹ U_h)
              = hid_rate

So: hid_rate = amb_rate + λn, hence **λ = −vis_rate / n**.

This is an exact scalar identity verified to machine precision. It determines the Lagrange multiplier from the visible rate alone.

### S = Ḣ − (vis/n)H is NOT an identity (0/100, residuals O(1)–O(10))

Despite the scalar trace identity, the full tensor equation fails for random (H, C, Ḣ). Residuals of order 1–10 across 100 trials. This means the tensor equation is genuinely selective — it holds only at configurations satisfying additional structure.

---

## 2. The constraint problem

### det(H) = δ is ILL-POSED

Numerical optimisation (10 trials, Part 6) shows vis_rate diverging to 10¹⁷–10²⁸. The optimiser drives H toward extreme condition numbers (one eigenvalue huge, others tiny) while maintaining fixed determinant.

**Mechanism:** With det(H) = δ, eigenvalues can be redistributed arbitrarily. Concentrating precision in the visible direction while collapsing the hidden metric R makes vis_rate unbounded.

**Conclusion:** The determinant constraint allows arbitrary degeneracy. It is the wrong constraint for a well-posed field equation.

### Tr(H) = τ is well-posed (theoretical argument)

With Tr(H) = τ, eigenvalues are bounded above (by τ) and below (by positivity). The feasible set {H ∈ SPD(n) : Tr(H) = τ} is compact. A continuous function on a compact set attains its maximum. So the problem has a finite solution.

The stationarity condition changes from:

    ∇f = λ H⁻¹       (det constraint — allows anisotropic H⁻¹ response)

to:

    ∇f = μ I           (trace constraint — forces isotropic response)

Explicitly:

    −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = μ I

This says the optimal ambient metric distributes its gradient response isotropically. The hidden-sector backreaction (second term) compensates the ambient perturbation response (first term) up to a uniform scalar.

### Which constraint is forced by the geometry?

**Open question.** Candidates:

| Constraint | Stationarity | Interpretation | Status |
|-----------|-------------|----------------|--------|
| det(H) = δ | ∇f = λ H⁻¹ | Volume normalisation | Unbounded — rejected |
| Tr(H) = τ | ∇f = μ I | Total precision budget | Well-posed (compact) — candidate |
| Tr(H⁻¹) = σ | ∇f = μ H⁻² | Total uncertainty budget | Well-posed — candidate |
| Tr(H) + Tr(H⁻¹) = κ | ∇f = μ(I − H⁻²) | Balanced budget | Well-posed — candidate |

The trace constraint Tr(H) = τ is the most natural: it fixes the total available precision, which is the physical resource the observer geometry distributes.

---

## 3. Mixed channel unification (VERIFIED, 60/60 on all sub-checks)

The three candidate field-equation objects are **not** the same tensor. They are three projections of one Codazzi-constrained mixed coupling channel:

    K(u,v) := β(v) θ(u) ~ B_v R⁻¹ B_u^T

### Curvature (antisymmetric off-diagonal)

    F_α(u,v) = K(v,u) − K(u,v) = −Alt(K)

Verified antisymmetric: 60/60.

### Source / hidden defect (symmetric diagonal)

    Q̂(u) = K(u,u) = B_u R⁻¹ B_u^T

Verified PSD: 60/60. Verified Q̂ = BR⁻¹B^T from information budget: 60/60.

### Variational stress (ambient lift)

    S = H Z R⁻¹ U_h R⁻¹ Z^T H = ambient lift of K(u,u)

Consistent with hidden_stress function: 60/60.

### The structural hierarchy

1. **Primary object:** mixed channel (β, θ)
2. **Compatibility:** Codazzi (Dβ = Dθ = 0)
3. **Field strength:** F_α = −Alt(K) — antisymmetric wedge
4. **Source law:** A_cpl = V⁻¹/² Diag(K) V⁻¹/² — symmetric diagonal contraction
5. **Variational balance:** S = ambient lift(K) + constraint term

---

## 4. Summary table

| Claim | Status | Evidence |
|-------|--------|----------|
| Conservation vis + hid = amb | PROVED | 50/50, errors < 1e-15 |
| Gradient formula | PROVED | 40/40, errors < 1e-9 |
| λ = −vis_rate/n (scalar) | PROVED | 50/50, errors < 1e-16 |
| S = Ḣ − (vis/n)H (tensor) | NOT an identity | 0/100, residuals O(1) |
| F_α antisymmetric from K | VERIFIED | 60/60 |
| Q̂ PSD from K | VERIFIED | 60/60 |
| Q̂ = BR⁻¹B^T | VERIFIED | 60/60 |
| det(H) field eq well-posed | REJECTED | vis_rate → ∞ (10/10 diverge) |
| Tr(H) field eq well-posed | THEORETICAL | Compact feasible set |

## 5. Honest boundaries

- The "field equation" S = Ḣ + λH has the structural shape of Einstein's equation but is **algebraic in H** (not a PDE), **observer-relative** (not intrinsic), and requires a well-posed constraint to be meaningful
- The determinant constraint produces an unbounded problem — no finite optimal H exists
- The trace constraint is the natural well-posed candidate but has not yet been verified numerically
- No propagator, conservation law for S, or observer-independent formulation has been shown
- The Einstein analogy is structural shape only: hidden stress = perturbation + constraint term
- The genuine variational content lives in the gradient equation ∇f = μ · (constraint gradient), not in the tensor rewriting S = Ḣ + λH
