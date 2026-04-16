# Field Equation — Hydrogen Atom (n=2, m=1)

**Date:** 16 April 2026  
**Status:** Analytical solution  

---

## Setup

Observer C = [cos φ, sin φ], hidden Z = [−sin φ, cos φ]^T.

In the rotated frame (observer-adapted):

    H̃ = [[a, b], [b, d]]     with Tr(H) = a + d = τ, Δ = ad − b² > 0
    Ḣ̃ = [[h, g], [g, f]]     (given perturbation, rotated)

Split:
    Φ = a − b²/d (visible precision)
    R = d (hidden metric)
    U_h = f (hidden first jet)

## The three stationarity equations

From ∇_H vis_rate = μ I on Tr(H) = τ:

### (12): Off-diagonal

    b(dh + af) = g(ad + b²)

This constrains the coupling b as a function of (a, d) and the perturbation (h, g, f).

**Key consequence:** If the perturbation is diagonal in the observer basis (g = 0), the optimal H is block-diagonal (b = 0). If the perturbation has off-diagonal structure (g ≠ 0), the optimal H must have nonzero coupling.

### (11) = (22): Diagonal balance

    −(d²h − 2bdg + b²f)/Δ² = −(b²h − 2abg + a²f)/Δ² + f/d²

Rearranging:

    [(b² − d²)h − 2bg(a − d) + (a² − b²)f] / Δ² = f/d²

### Plus: a + d = τ (trace constraint)

Three equations in three unknowns (a, b, d). Generically this gives a discrete set of solutions (finitely many optimal H*).

## Special case: diagonal perturbation (g = 0)

When g = 0, the (12) equation gives b = 0. The system reduces to:

    [(−d²h + b²f) − (−b²h + a²f)] / Δ² = f/d²

With b = 0, Δ = ad:

    [−d²h + a²f] / (ad)² = f/d²
    [a²f − d²h] / a²d² = f/d²
    f/d² − h/a² = f/d²
    −h/a² = 0

So h = 0. This means: if the perturbation is diagonal (g = 0) AND has h = 0 (no visible component), the field equation is satisfied for any (a, d) with a + d = τ.

If h ≠ 0 and g = 0: no solution with b = 0. We need to recheck — but (12) forces b = 0 when g = 0. So h = 0 is required. This means:

> **For a purely hidden perturbation (Ḣ = f · e₂e₂^T), any H with Tr(H) = τ and b = 0 is a critical point.** The entire τ-hyperplane of diagonal H is critical. The observer sees nothing (vis_rate = 0 from the hidden perturbation at the adapted observer).

> **For a perturbation with visible component (h ≠ 0) and no coupling (g = 0), there is no critical point.** The vis_rate is monotone in (a, d).

This makes physical sense: a purely visible perturbation (h ≠ 0, f = 0, g = 0) has vis_rate = h/a (in the diagonal case), which increases monotonically as a → 0. The trace constraint bounds a > 0 but the maximum is at the boundary a → 0, not at an interior critical point.

## Special case: coupled spring

k₁ = 3, k₂ = 2, kc = 5. In the lab frame:

    H = [[8, −5], [−5, 7]]
    Ḣ = [[1, −1], [−1, 1]]

The optimal observer is at φ = 135° (from Part 2 of the numerics). In the rotated frame:

    c = cos(135°) = −1/√2, s = sin(135°) = 1/√2
    Q = [[-1/√2, -1/√2], [1/√2, -1/√2]]

    H̃ = Q^T H Q = [[12.5, 0.5], [0.5, 2.5]]
    Ḣ̃ = Q^T Ḣ Q = [[2, 0], [0, 0]]

So h = 2, g = 0, f = 0, τ = 15.

With g = 0 and f = 0: the (12) equation gives b(2d) = 0, so b = 0 or d = 0. Since d > 0, we need b = 0.

The diagonal balance becomes: −d²·2/Δ² = −b²·0/Δ² + 0, with b = 0, Δ = ad:
    −2d²/(ad)² = 0
    −2/(a²) = 0

This has no solution. Consistent with our analysis: a purely visible perturbation has no interior critical point. The vis_rate increases as a decreases (more precision allocated to visible direction).

**This confirms:** the trace-constrained problem for the coupled spring with a purely visible perturbation has its optimum at the BOUNDARY of the feasible set, not at an interior critical point.

## What this means

The hydrogen atom analysis reveals:

1. **Interior critical points exist only for perturbations with specific structure** — not generic perturbations
2. **Generic perturbations have boundary optima** — the maximum vis_rate occurs at the most anisotropic H allowed by the constraint
3. **The static field equation ∇f = μI has solutions only in special cases**

This strengthens the case for the **dynamic Lagrangian formulation** as the right framework. The kinetic term ||Ḣ||² prevents H from jumping to boundary configurations, and the field equation becomes a second-order ODE rather than an algebraic equation.

## The definitive structural picture

| Formulation | Well-posed? | Solutions? | Type |
|-------------|-------------|------------|------|
| det(H) = δ | No (unbounded) | None | — |
| Tr(H) = τ | Yes (compact) | Boundary optima generically | Algebraic, degenerate |
| ∇f = 0 (unconstrained) | — | Only for rank-1 perturbations | Over-determined |
| Lagrangian with kinetic term | Yes | Interior solutions | 2nd-order ODE |

**The Lagrangian formulation is the only one that gives genuine interior solutions for generic perturbations.** This is the correct field equation.
