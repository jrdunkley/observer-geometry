# Field Equation — Theoretical Analysis

**Date:** 16 April 2026  
**Status:** Active research  

---

## 1. Why the static variational problem is ill-posed

The problem "maximise vis_rate(H, C, Ḣ) over H at fixed det(H)" has no solution because vis_rate is unbounded. The optimiser drives H toward degenerate configurations (extreme condition number) that concentrate all precision in the visible direction.

This is not a numerical artefact. It reflects a structural fact:

> **Without a cost for changing H, the observer can extract unlimited information by reshaping the ambient metric.**

This is analogous to saying "without a speed-of-light constraint, you can transmit unlimited information." The constraint must come from the geometry.

## 2. Three candidate regularisations

### 2a. Static constraint: Tr(H) = τ

The trace constraint fixes the total precision budget. The feasible set is compact, so the problem is well-posed. The field equation becomes:

    −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = μ I

This says the gradient of vis_rate must be isotropic at the optimum. Interpretation: the optimal ambient metric distributes its sensitivity uniformly across all directions.

**Limitation:** This is a static snapshot. It doesn't explain *why* H has a particular trace, or how H evolves.

### 2b. Kinetic regularisation: the joint Lagrangian

The coupling-budget Lagrangian from the 0.3.3 session governed observer paths C(t). Generalising to joint (H(t), C(t)) evolution:

    L[H, C] = ∫₀^T [ vis_rate(H, C, Ḣ) − (α/2)||Ḣ||²_H − (β/2)||Ċ||²_Gr ] dt

where:
- ||Ḣ||²_H = Tr(H⁻¹ Ḣ H⁻¹ Ḣ) is the natural Riemannian metric on SPD(n)
- ||Ċ||²_Gr is the Grassmannian metric on Gr(m,n)
- α, β > 0 are inertia parameters

The kinetic term ||Ḣ||² prevents instantaneous reshaping. The field equation is the Euler-Lagrange equation:

    d/dt[∂L/∂Ḣ] − ∂L/∂H = 0

The ∂L/∂Ḣ term gives a "momentum" conjugate to H, and the d/dt derivative gives acceleration. This produces a **second-order equation for H** — unlike the static case which is algebraic.

**Key structural point:** The kinetic term on SPD(n) is the Fisher-Rao metric. This is the unique Riemannian metric on SPD(n) invariant under congruence H → AHA^T. So the regularisation is geometrically forced, not arbitrary.

### 2c. Codazzi constraint

The Codazzi identity Dβ = 0 removes m(n−m) directions from Sym(n). If H evolves along a path, the Codazzi identity constrains which Ḣ are compatible with the current observer C.

This is an evolution constraint, not a static one. It says:

> Not every ambient perturbation Ḣ is compatible with a given observer motion Ċ.

Codazzi-compatible Ḣ live in a subspace of dimension:
    dim Sym(n) − m(n−m) = n(n+1)/2 − m(n−m)

This removes ~40% of directions for balanced m ~ n/2.

## 3. The natural formulation

The three regularisations combine into one picture:

> **H evolves along Codazzi-compatible directions to extremise the joint Lagrangian, whose kinetic term (Fisher-Rao metric on SPD(n)) provides the natural regularisation.**

The field equation is then:

    α · ∇_Ḣ ||Ḣ||²_H + ∇_H vis_rate = 0
    subject to: Codazzi compatibility Dβ = Dθ = 0

Expanding:
    α · d/dt[H⁻¹ Ḣ H⁻¹] + (−H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T) = 0

This is a **second-order ODE on SPD(n)** constrained to the Codazzi-compatible subspace.

### Zero-velocity limit

When Ḣ = 0 (static H), the kinetic term vanishes and we recover the static field equation ∇_H vis_rate = 0, i.e., the gradient vanishes without any constraint. This gives:

    H⁻¹ Ḣ H⁻¹ = Z R⁻¹ U_h R⁻¹ Z^T

which says: **the ambient perturbation response exactly equals the hidden-sector backreaction**. This is the genuine static field equation — no constraint needed, just vanishing gradient.

But wait: this is the equation ∇f = 0 (unconstrained), not ∇f = λ·(something). Is this well-posed?

On SPD(n) without constraint, ∇f = 0 could have solutions if f has critical points. The unboundedness on det(H) = const doesn't imply unboundedness on all of SPD(n) — it implies unboundedness along certain directions. The gradient could vanish at saddle points.

## 4. Connection to existing results

### The adapted observer already solves the C-equation

At the adapted observer C*, B = L^TḢZ = 0. This kills the hidden coupling entirely. The source law reduces to A_cpl = A_direct. The observer field equation is solved exactly by the adapted observer.

### The joint (H*, C*) system

At the joint optimum:
- C* is the adapted observer (B = 0)
- H* satisfies the ambient field equation

When B = 0, the gradient simplifies:
    ∇_H f = −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T

But at B = 0, the coupling between visible and hidden sectors vanishes. The vis_rate becomes:

    vis_rate = Tr(Φ⁻¹V) where V = L^TḢL is the visible first jet

And hid_rate = Tr(R⁻¹U_h) is determined by the hidden sector alone.

The gradient at B = 0 tells us how to adjust H to further improve vis_rate. The hidden-sector term Z R⁻¹ U_h R⁻¹ Z^T acts as a "source" that the ambient perturbation response must balance.

### The hidden stress at the adapted observer

When C = C* (adapted), the hidden stress S = HZ R⁻¹ U_h R⁻¹ Z^T H still depends on H. The tensor equation S = Ḣ + λH would become a constraint on H if we could determine λ from the geometry.

But we showed λ = −vis_rate/n as a trace identity. At B = 0, vis_rate simplifies, so λ is determined by the visible sector alone.

## 5. The key open question

**Does the unconstrained gradient equation ∇_H vis_rate = 0 have solutions on SPD(n)?**

If yes: these are the natural "equilibrium" ambient metrics for a given (C, Ḣ). No additional constraint is needed.

If no: the kinetic regularisation is essential, and the field equation is dynamical (second-order ODE on SPD(n)).

This can be checked analytically for n=2, m=1 (the hydrogen atom) without heavy computation.

## 6. The hydrogen atom: n=2, m=1

For C = [cos φ, sin φ], Z = [−sin φ, cos φ]^T:

    R = Z^THZ    (scalar)
    U_h = Z^TḢZ  (scalar)

The gradient is:
    ∇_H f = −H⁻¹ Ḣ H⁻¹ + (U_h/R²) Z Z^T

Setting ∇f = 0:
    H⁻¹ Ḣ H⁻¹ = (U_h/R²) Z Z^T

The LHS is a rank-n matrix (generically full rank). The RHS is rank 1. So equality requires the LHS to also have rank 1, i.e., H⁻¹Ḣ has rank 1.

For generic Ḣ with rank > 1, **this equation has no solution**. The unconstrained critical point does not exist for generic perturbations.

**Conclusion:** The unconstrained problem has no critical points generically. The field equation requires either:
1. A constraint (Tr(H) = τ is the natural choice), or
2. A kinetic term (the Lagrangian formulation)

This is the definitive answer to the constraint question: **the constraint is forced by the non-existence of unconstrained critical points.**

## 7. The field equation (final form)

### Static (with trace constraint)

    −H⁻¹ Ḣ H⁻¹ + Z R⁻¹ U_h R⁻¹ Z^T = μ I

where μ = (Tr(gradient)) / n. This selects H* on the trace hypersurface.

### Dynamic (Lagrangian)

    α · d/dt[H⁻¹ Ḣ H⁻¹] = H⁻¹ Ḣ H⁻¹ − Z R⁻¹ U_h R⁻¹ Z^T
    subject to Codazzi: Dβ = Dθ = 0

This is the second-order evolution equation on SPD(n).

### Structural comparison

| GR | Observer geometry |
|----|-------------------|
| G_μν = 8πG T_μν + Λg_μν | ∇f = μ I (static) |
| Einstein-Hilbert action | Joint Lagrangian with Fisher-Rao kinetic term |
| Bianchi identity ∇·G = 0 | Codazzi Dβ = Dθ = 0 |
| Geodesic equation | Euler-Lagrange on Gr(m,n) |
| Metric g_μν | Ambient precision H |
| Stress-energy T_μν | Hidden stress S |
| Cosmological constant Λ | λ = −vis_rate/n |

---

## Next steps

1. **Verify Tr(H) = τ field equation numerically** (lightweight: n=3, m=1, should converge)
2. **Derive the dynamic Lagrangian EL equation explicitly** for n=2, m=1
3. **Check whether Codazzi compatibility is automatic** for the Tr(H)-constrained solution
4. **Classify the solution space** for the static field equation
