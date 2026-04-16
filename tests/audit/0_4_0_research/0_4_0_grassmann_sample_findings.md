# Grassmannian Optimisation: 5-Molecule Sample Results

**Date:** 15 April 2026
**Code:** `0_4_0_grassmann_sample.py` (3131s on 20 cores)
**Data:** 5 matched water-vacuum molecules, m = 2, 3, 5

---

## The adapted observer does NOT match the Grassmannian optimum

The earlier Spearman 0.83 (from 50-trial random search) was misleading. With proper L-BFGS-B Grassmannian optimisation (15 restarts, 150 iterations each), the true dynamic optimum is substantially better than the adapted observer on most molecules.

| Molecule | n_vib | Adapted vf (m=2) | Grassmann(vf) (m=2) | Gap |
|----------|-------|-------------------|---------------------|-----|
| 011106 | 36 | 4.79 | **5.85** | -1.06 |
| 011777 | 39 | 1.23 | **1.39** | -0.17 |
| 030178 | 42 | **1.54** | 1.52 | +0.02 |
| 028438 | 48 | **-18.94** | **6.15** | **-25.09** |
| 013023 | 48 | **1.65** | 1.63 | +0.02 |

On 2/5 molecules the adapted observer is near-optimal (gap < 0.05). On 2/5 it's moderately worse. On 1/5 (028438) it's catastrophically anti-aligned (vis_frac = -19 vs +6).

## V > 0 rates by observer type

At m=2:

| Observer | V > 0 count | Rate |
|----------|-------------|------|
| Canonical | 0/5 | 0% |
| Adapted | 2/5 | 40% |
| Grassmann(vf) | 4/5 | **80%** |
| Grassmann(v_min) | 4/5 | **80%** |
| Grassmann(combined) | 3/5 | 60% |

Proper optimisation achieves V > 0 on 80% of molecules — far higher than the 8.6% canonical rate or the 38% random-search rate.

## The adapted observer has high leakage

| Molecule | Adapted leakage (m=2) | Grassmann(combined) leakage | Ratio |
|----------|----------------------|----------------------------|-------|
| 030178 | 61.5 | 7.9 | 7.8x |
| 013023 | 37.4 | 11.5 | 3.3x |
| 011777 | 23.3 | 10.2 | 2.3x |

The adapted observer consistently has the highest leakage of all five methods. This is because it was designed for static visibility of a single family member, not for leakage minimisation or dynamic capture.

## The combined objective is best overall

Grassmann(combined) = vis_frac + 0.3 * max(v_min, 0) - 0.1 * leakage. It consistently achieves:
- Reasonable vis_frac (not the highest, but positive)
- Low leakage (often the lowest of all methods)
- Moderate V > 0 rate

This is the objective that nomosteer should use by default.

## m-dependence

V > 0 rates decrease with m across all observers:
- At m=2: Grassmann(vf) achieves V > 0 on 4/5
- At m=3: 2/5
- At m=5: 1/5

Higher observer rank makes it harder to enter the exact sector. This is expected: more visible dimensions means more eigenvalues of V must be positive simultaneously.

## Key correction to earlier results

The static-dynamic Spearman of 0.83 was computed with 50-trial random search as the "dynamic optimum." The true Grassmannian optimum (from proper L-BFGS-B) is substantially better. The adapted observer is:
- Near-optimal on 40% of molecules
- Moderately suboptimal on 40%
- Catastrophically wrong on 20%

The unification is partial, not complete. The static design criterion captures the right general direction but misses the fine structure. **For nomosteer, proper Grassmannian optimisation is essential — the adapted observer is a warm start, not the answer.**
