# Anti-Correlation Investigation: Tr(A_cpl^2) vs vis_frac

**Date:** 15 April 2026
**Code:** `0_4_0_anticorrelation_proof.py`

---

## Result: NOT a universal theorem

The anti-correlation between Tr(A_cpl^2) and vis_frac is strong in some configurations but **reverses sign** in others. It is not a theorem.

### Dimension sweep

| (n,m) | n_valid | Spearman(A^2, vf) | Median vf |
|-------|---------|-------------------|-----------|
| (3,1) | 320 | **-0.94** | 0.85 |
| (4,1) | 261 | **+0.98** | -2.46 |
| (4,2) | 109 | -0.76 | 1.21 |
| (5,1) | 321 | **-0.94** | 0.53 |
| (5,2) | 46 | +0.62 | -1.96 |
| (6,2) | 51 | -0.64 | 4.93 |

**Pattern:** When median vis_frac > 0 (observer aligned), the anti-correlation holds. When median vis_frac < 0 (observer anti-aligned), the correlation reverses.

### Perturbation sweep

| Hdot | n_valid | Spearman | Median vf |
|------|---------|----------|-----------|
| 0 | 231 | -0.72 | 0.63 |
| 2 | 143 | -0.60 | 0.77 |
| 6 | 54 | -0.74 | 1.60 |
| 9 | 20 | +0.55 | -2.13 |

Same pattern: positive vis_frac → anti-correlation; negative vis_frac → positive correlation.

### Why the sign flips

The ratio vis_frac = vis_rate / amb_rate has a sign pathology. When both vis_rate and amb_rate are negative (which happens when Phi is decreasing), the ratio is positive. When A_cpl is large (strong hidden coupling), it pushes vis_rate further negative. If amb_rate is also negative, vis_frac = (negative)/(negative) = positive, and it INCREASES with A_cpl.

The anti-correlation is real when measured in the aligned regime (vis_frac > 0, amb_rate > 0). In the anti-aligned regime, the ratio flips.

### The underlying mechanism is real but weak

||B|| (the coupling norm) correlates with hid_rate at Spearman 0.26 and anti-correlates with vis_frac at -0.26. The conservation law (large hid => small vis) provides the mechanism, but it's only a moderate effect — not the dominant one.

### What IS true

The following are theorems (proved earlier in the session):
1. **A_cpl >= 0 on linear paths** (positivity theorem)
2. **A_cpl = 0 iff B = 0** (zero coupling = zero hidden defect)
3. **The conservation law** vis + hid = ambient
4. **The equivalence principle** (A_cpl can always be locally cancelled)

The anti-correlation is an empirical pattern in the aligned regime, not a structural theorem.
