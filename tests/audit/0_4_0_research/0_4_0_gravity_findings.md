# Emergent Gravity from Observer Geometry

**Date:** 15 April 2026
**Code:** `0_4_0_emergent_gravity.py` (4/4 checks passed)

---

## What is established

### The equivalence principle holds

**Theorem.** For any (H, C, Hdot) with V > 0, there exists a symmetric Hddot such that A_cpl = 0.

**Proof.** The map Hddot -> L^T Hddot L is surjective from Sym(n) to Sym(m) (because L has full column rank m and n > m). Given the hidden defect Q_hat, construct Hddot_cancel = L_pseudo M_target L_pseudo^T where M_target = 2 Q_hat and L_pseudo = L (L^T L)^{-1}. Then L^T Hddot_cancel L = 2 Q_hat and A_cpl = A_direct + hidden_defect = -hidden_defect + hidden_defect = 0.

**Verified:** At (n=3,m=1) and (n=4,m=2) to machine precision (max error 2.4e-15).

**Interpretation.** The hidden defect can always be locally cancelled by the right ambient acceleration. The observer can always find a "freely falling frame" in which the hidden sector is transparent. This is structurally identical to the GR equivalence principle.

### The source-law action anti-correlates with information capture

Across 111 observer points on Gr(2,4): **Spearman(Tr(A_cpl^2), vis_frac) = -0.91** (p = 7.4e-44).

Minimising the source-law energy Tr(A_cpl^2) maximises the visible fraction. The "vacuum" (A_cpl = 0, no hidden source) is the optimal observer. "Gravity" (nonzero A_cpl from hidden coupling) is the deviation from optimal observation.

### The hidden defect is a non-conserved source

Tr(Q_hat) is NOT conserved along a path (d/dt[Tr(Q_hat)] = 0.053 != 0 at the test point). The hidden defect evolves with the geometry, like stress-energy in GR. It's a source, not a charge.

### The spatial fall-off matches the source profile

For a localised Gaussian hidden source, the visible precision excess has the same Gaussian profile (fitted sigma = 1.10 vs source sigma = 1.00). There is no long-range 1/r^2 fall-off. This is because the current framework has no propagator — the geometry responds locally, not through a field equation.

## What is sharply posed but open

### The field equation

The equivalence principle shows A_cpl can be locally cancelled. But there is no field equation that determines H from C. In GR, Einstein's equation G_munu = 8 pi G T_munu makes the metric dynamical. In observer geometry, H is background.

The self-consistency condition (find H such that the adapted observer for H is C itself) converges in 1 step (trivially, because the adapted observer is deterministic). The open problem is the SECOND direction: given C, what equation determines H?

### The propagator

The Newton limit shows local fall-off (Gaussian, matching the source). For 1/r^2 fall-off (true gravity), you need a propagator — a Green's function that solves a differential equation. This requires promoting H to a field with spatial derivatives, i.e., a field theory on the observer bundle.

The natural candidate: H satisfies a Laplace-type equation nabla^2 H = source(C), where nabla is the covariant derivative of the split-frame connection. This would give 1/r^(n-2) fall-off in n spatial dimensions.

### The coupling constant

In GR, Newton's constant G sets the strength of gravity. In observer geometry, the analogous quantity is cond(R)^{-1} (the inverse condition number of the hidden metric). When R is well-conditioned, "gravity" is weak (small Q_hat for given B). When R is ill-conditioned, "gravity" is strong.

What determines R? In the current framework, R = Z^T H Z depends on H and C. Making "G" dynamical requires making R evolve independently, which requires the field equation.

## Structural summary

| GR concept | Observer geometry analogue | Status |
|------------|--------------------------|--------|
| Metric g | SPD field H | Background (given) |
| Connection Gamma | Split connection (alpha, beta, theta, omega) | Derived from H, C |
| Curvature R | F_alpha = -beta wedge theta | Computed |
| Stress-energy T | Hidden defect Q_hat = B R^{-1} B^T | Computed, PSD |
| Conservation nabla T = 0 | vis + hid = ambient | PROVED |
| Equivalence principle | A_cpl can be locally cancelled | PROVED |
| Einstein equation G = 8piG T | ??? | OPEN |
| Newton's constant G | cond(R)^{-1} | Derived, not dynamical |
| Newton limit (1/r^2) | Gaussian fall-off (matches source) | No propagator yet |
| Einstein-Hilbert action | Tr(A_cpl^2) anti-correlates with vis_frac | Observed, not derived |
