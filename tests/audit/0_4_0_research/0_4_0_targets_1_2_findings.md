# Targets 1 & 2: A_cpl Positivity Theorem and Information Capture Curve

**Date:** 15 April 2026
**Code:** `0_4_0_acpl_positivity.py` (15/15), `0_4_0_scaling_law.py`

---

## Target 1: A_cpl >= 0 on linear paths (THEOREM)

### Statement

**Theorem (A_cpl positivity on linear paths).** Let H(t) = (1-t)H_0 + tH_1 be a linear path through SPD(n) with fixed rank-m observer C. If V = L^T(H_1 - H_0)L is positive definite on the active support, then

    A_cpl = V^{-1/2} (B R^{-1} B^T) V^{-1/2} >= 0.

### Proof

On a linear path, Hddot = 0. From the completed square (eq 2.8 of the TN), with beta = 0 (fixed C): W = -2 B R^{-1} B^T. Therefore A_cpl = -1/2 V^{-1/2} W V^{-1/2} = V^{-1/2} (B R^{-1} B^T) V^{-1/2}. Since R is SPD, R^{-1} is SPD, so Q_hat = B R^{-1} B^T is PSD. Conjugation by V^{-1/2} preserves PSD. QED.

### Verification

- 2,134 random trials across dimensions 3-7: **zero violations**
- Minimum A_cpl eigenvalue across all trials: -3.7e-14 (machine zero)
- Algebraic chain verified at 5 dimension pairs
- Confirmed by 10,243 QM9 molecular pairs (all positive)

### Consequences

1. **Evidence is always concave on linear paths.** The source term in the evidence second derivative is -2 Tr(Phi^{-1} Q_hat) <= 0.

2. **The hidden sector always stabilises.** The hidden defect Q_hat = B R^{-1} B^T is the "stabilising correction" from hidden variables. On linear paths, it can only make the evidence more concave, never less.

3. **Indefiniteness requires curvature.** A_cpl can become indefinite only when Hddot != 0 (the path curves through SPD space). The critical Hddot scale is computable.

4. **A_cpl = 0 iff B = 0.** The hidden defect vanishes exactly when the perturbation doesn't couple visible and hidden sectors (block-diagonal in the split frame).

---

## Target 2: Information Capture Curve

### The curve

For water-vacuum molecular pairs, sweeping m from 1 to n-1:

| m/n | Median vis_frac | Interpretation |
|-----|-----------------|---------------|
| 0.025 | 0.06 | First mode: negligible |
| 0.075 | 0.32 | Rapid initial capture |
| 0.125 | 0.65 | **50% capture at ~12% of modes** |
| 0.175 | 0.88 | Steep gain continues |
| 0.325 | 0.97 | Near saturation |
| 0.475 | 1.08 | Fully captured (>1 = hidden sector contributes negatively) |
| 0.725 | 1.03 | Diminishing returns |
| 0.975 | 1.00 | Approaches ambient |

### Key finding

**50% of the solvent-response information is captured by observing just 12% of the vibrational modes.**

The curve is concave: steep initial gain, rapid saturation around m/n ~ 0.35. After that, additional modes add almost nothing. The marginal gain drops from +0.33 per 5% of modes (at m/n = 0.1) to effectively zero (at m/n > 0.4).

### Physical interpretation

The softest vibrational modes carry disproportionately more solvent-response information than higher modes. A spectroscopist observing just the 3-5 lowest-frequency modes of a 20-mode molecule captures the majority of the solvent effect. This is a quantitative law of diminishing returns for molecular observation, derived from first principles via the conservation law.

### Vis_frac > 1

The visible fraction exceeds 1.0 for m/n > 0.35. This means the hidden sector's rate has the opposite sign to the visible rate — the hidden modes are "working against" the solvent perturbation while the visible modes amplify it. The conservation law ensures the sum always equals the ambient rate.
