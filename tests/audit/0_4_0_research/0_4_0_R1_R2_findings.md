# R1 + R2 Findings: Information Destruction and Determinant Conservation

**Date:** 15 April 2026
**Code:** `0_4_0_R1_R2.py` (15/15 checks passed)

---

## R1: Information-destruction locus

### Definition

For a path H(theta) with fixed C, the **information-destruction locus** is the set where V = L^T dH L drops rank — where a visible direction goes dark.

### Coupled spring (n=2, m=1)

V = k2^2 / (k2 + kc)^2 > 0 always. The information-destruction locus is **empty**. The observer (mass 1) always receives information about the coupling change because it is directly coupled to mass 2 through the coupling spring. V -> 0 only as kc -> infinity (rigid limit where internal dynamics disappear).

### Two-mode squeezed state (n=4, m=2)

Both eigenvalues of V are negative throughout the squeezing range r in [0.05, 2.0]. This means increasing squeezing **always destroys** visible precision in every direction.

| Squeezing r | V eigenvalues | Phi min eigenvalue |
|-------------|--------------|-------------------|
| 0.05 | -0.20, -0.20 | 0.995 |
| 0.54 | -0.96, -0.96 | 0.607 |
| 1.04 | -0.48, -0.48 | 0.249 |
| 2.00 | -0.07, -0.07 | 0.037 |

**Physical interpretation:** The price of entanglement is visible information. To entangle two modes, you must sacrifice local precision. The information-destruction locus is at r = 0 (no squeezing). Beyond that point, every increase in squeezing destroys visible Fisher information. The source law cannot operate (V < 0, not on the support-stable stratum), but this is itself informative: **the failure of the support condition IS the physical content**. It says "you are moving in a direction that costs information."

---

## R2: Split determinant conservation law

### The identity

**Theorem (split determinant conservation).** For M = [L Z] with M^T H M = diag(Phi, R):

    det(Phi) * det(R) = det(H) * det(M)^2

Taking logs and differentiating along a path H(t):

    d/dt[log det Phi] + d/dt[log det R] = d/dt[log det H] + 2 d/dt[log det M]

Or in words:

    (visible evidence rate) + (hidden evidence rate) = (ambient rate) + (frame transport rate)

### Verification

Verified to machine precision (max error 2.7e-15) along paths at three dimensions: (3,1), (4,2), (5,3).

### The fixed-C simplification

When C is constant (Z fixed, frame transport rate = 0):

    visible rate + hidden rate = ambient rate

This says: **the total information rate Tr(H^{-1} Hdot) splits exactly into what the observer sees and what's hidden.** No information is created or destroyed by the observation — only redistributed.

### Quantitative examples

| (n,m) | Visible rate | Hidden rate | Sum | Ambient rate |
|-------|-------------|-------------|-----|-------------|
| (3,1) | -0.432 | -0.101 | -0.533 | -0.533 |
| (4,2) | -0.927 | -0.495 | -1.422 | -1.422 |
| (5,3) | -0.087 | -0.617 | -0.704 | -0.704 |

Note the varying split ratios: at (5,3) the observer sees 3 of 5 dimensions but captures only 12% of the total rate. This quantifies observer inefficiency.

### Connection to source law

The visible rate is:

    d/dt[log det Phi] = Tr(Phi^{-1} V) + 2 Tr(alpha) = f'(0)

which is governed by V (the visible first jet from the source law). The conservation law therefore links the source law to the total ambient information change.

### Connection to evidence geometry

In the TN4 evidence framework, the log-evidence involves 1/2 log det(Fisher information). The conservation law says that along a model-parameter path, the evidence contribution from visible parameters plus the evidence contribution from hidden parameters equals the total Fisher information change. This is the structural reason why exact hidden elimination (TN4 Prop 2.1) preserves the evidence: nothing is lost, only split.

---

## Interpretation

The conservation law and the information-destruction locus together give the **information budget** of an observation:

1. The total information rate Tr(H^{-1} Hdot) is the "budget" — fixed by the ambient geometry.
2. The conservation law splits it: visible + hidden = total.
3. The information-destruction locus marks where the visible share drops to zero.
4. The source law A\_cpl governs the acceleration of the visible share.

This is the precise version of "every optimisation problem is an alignment problem with the underlying geometry": the geometry determines the total information budget, the observer determines the split, and the source law governs how the split evolves.
