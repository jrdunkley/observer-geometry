# T1 + T2 Findings: Source Law on Physical Systems + Curvature-Evidence Inequality

**Date:** 15 April 2026
**Code:** `0_4_0_T1_T2.py`

---

## T1: Source law on coupled physical systems

### Setup

Coupled spring: H(kc) = [[k1+kc, -kc], [-kc, k2+kc]], C = [1,0] (observe mass 1).
Coupled LC: H(gc) = [[g1+gc, -gc], [-gc, g2+gc]], C = [1,0] (observe node 1).
Both: n=2, m=1, parameterised by coupling strength.

### Results

**Analytical formulas verified:**
- Phi(kc) = (k1\*k2 + k1\*kc + k2\*kc) / (k2 + kc) (error < 4.4e-15)
- dPhi/dkc = k2^2 / (k2+kc)^2 > 0 always (error < 1.1e-16)

**Physical predictions from the source law:**

1. **Coupling always increases visible precision.** V > 0 across the entire range. Adding a coupling spring to mass 1 always makes it stiffer (more informative).

2. **A_cpl > 0 everywhere** (no sign changes). The evidence surface is always concave — the system is in the "near peak" regime at every coupling strength.

3. **Source term decays with coupling.** At kc=0.01: source = -0.33; at kc=30: source = -0.00005. The evidence becomes less sensitive to coupling changes at high coupling (diminishing returns).

4. **Hidden load increases monotonically.** From 0 at kc=0 to 5.77 at kc=30. More coupling means more hidden structure (mass 2 increasingly affects what mass 1 sees).

5. **Substrate-agnostic.** At matched parameters (k1=g1, k2=g2, kc=gc), A_cpl is identical for springs and LC circuits to machine precision. The source law sees only the mathematical structure (H, C, dH), confirming the cross-substrate prediction from the operationalise layer.

### Asymptotic behaviour

| Regime | Phi | V | A_cpl | Hidden load |
|--------|-----|---|-------|-------------|
| kc -> 0 | k1 (isolated) | 1 | 0.50 | 0 |
| kc -> inf | k1+k2 (rigid block) | 0 | 0.03 | large |

As coupling increases, the system transitions from "isolated spring" to "rigid block." The source law tracks this transition through A_cpl, which decreases monotonically from 0.5 to 0.03.

---

## T2: Curvature-evidence inequality

### The coupling identity

In the mixed 2-parameter case (H varies in s, C varies in t):

    ||F(ds,dt)||^2 = Tr(G_beta * R^{-1} B^T B R^{-1})

where G\_beta = beta\_t^T beta\_t (hidden Gram of the C-perturbation) and B = L^T H\_s Z (source coupling).

### Connection to A_cpl

The hidden defect Q\_hat = B R^{-1} B^T appears in both:
- **Curvature:** ||F||^2 involves Q\_hat through R^{-1} B^T B R^{-1}
- **Source law:** A\_cpl = A + V^{-1/2} Q\_hat V^{-1/2} (when beta = 0)

For m=1 (scalar case):

    ||F||^2 = ||beta_t||^2 * V * (A_cpl - A) / R

This is the **curvature-source coupling identity**: the mixed curvature equals the observer coupling strength times the source law's hidden correction, scaled by V/R.

### The inequality

    ||F||^2 <= ||beta_t||^2 * cond(R) * m * ||A_cpl - A|| * ||V||

**Interpretation:** High curvature requires BOTH strong observer coupling (large ||beta||) AND strong source-hidden coupling (large hidden defect A\_cpl - A). If either is zero, the curvature vanishes. This is the structural reason why the source and curvature sectors interact only through the shared split channel.
