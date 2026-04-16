# R5: Optimal Observer Alignment

**Date:** 15 April 2026
**Code:** `0_4_0_R5_optimal_observer.py` (1 check passed, landscape + 5 alignment trials)

---

## The optimisation problem

Given fixed ambient geometry (H, Hdot, Hddot), find the observer C that maximises the visible evidence acceleration f''_vis. The second-order conservation law constrains:

    f''_vis = f''_amb - f''_hid

So maximising f''_vis is equivalent to minimising f''_hid.

## Results

### Landscape on S^2 (n=3, m=1)

The full landscape of f''_vis was computed at 5000 points on S^2. Conservation holds at every point (error 0.0).

- Range of f''_vis: [-0.015, +0.068]
- Ambient acceleration: 0.095
- Optimal observer captures 71% of the ambient acceleration
- Worst observer captures -16% (negative = the observer is working against the geometry)

### Optimal observer structure

The optimal observer does NOT simply align with the largest eigenvector of H or Hdot. Instead, it aligns with a specific eigenvector of the **effective acceleration matrix**:

    A_eff = H^{-1} Hddot - (H^{-1} Hdot)^2

which is the matrix whose trace gives the ambient acceleration f''_amb.

Across 5 random trials, the optimal observer consistently shows high alignment (0.71 to 0.998) with an eigenvector of A_eff, but which eigenvector varies. The selection depends on the interplay between the eigenvalue magnitudes and the hidden-sector geometry.

### Key finding: the observer is not a source

The landscape confirms quantitatively what the conservation law says structurally: the observer redistributes acceleration between visible and hidden sectors. The total is fixed at f''_amb = 0.095. The best observer can capture 71% of this. The worst observer can make the visible acceleration negative (working against the geometry by pushing information into the hidden sector).

This is the precise sense in which "every optimisation problem is an alignment problem with the underlying geometry." The geometry (A_eff) determines the acceleration budget. The observer's job is to align with it.

### Alignment is not eigenvalue selection

The alignment principle does NOT reduce to "pick the largest eigenvalue of A_eff." Trial 0 aligns with the middle eigenvalue (-0.22). Trial 2 aligns with the largest (0.64). Trial 3 aligns with the smallest (-11.43). The optimal observer balances the eigenvalue magnitude against the cost of the hidden-sector projection, which depends on R = Z^T H Z.

This means observer optimisation requires solving the full conservation-constrained problem on the Grassmannian, not just an eigendecomposition. The TN1 weighted-family stationarity equation is the right framework — the conservation law provides the additional constraint that A_cpl must respect.

---

## Connection to TN1 observer design

The TN1 weighted-family frontier stationarity equation (Theorem 2.2) finds observers that are stationary for a visibility-leakage trade-off:

    F_mu(P) = S(P) - mu * L(P)

The second-order conservation adds a dynamic constraint: among observers that satisfy the TN1 stationarity condition, the one that best aligns with A_eff will have the most stable evidence geometry (largest f''_vis in the peak regime, or most controlled f''_vis in the saddle regime).

This suggests a two-stage observer design:
1. **Static stage:** Use TN1 weighted-family frontier to find observers that maximise visibility for a given leakage tolerance.
2. **Dynamic stage:** Among the frontier observers, use the conservation law to select the one whose evidence acceleration is best aligned with the ambient geometry.

The conservation law provides the missing second-stage criterion.
