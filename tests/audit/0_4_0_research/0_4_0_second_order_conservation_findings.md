# Second-Order Conservation Law

**Date:** 15 April 2026
**Code:** `0_4_0_second_order_conservation.py` (21/21 checks passed)

---

## The theorem

**Second-order information conservation.** For a path H(t) through SPD(n) with fixed observer C:

    d^2/dt^2 [log det Phi] + d^2/dt^2 [log det R] = d^2/dt^2 [log det H]

Or in words:

    (visible evidence acceleration) + (hidden evidence acceleration) = (ambient acceleration)

This holds exactly (machine precision, residual < 5e-16) at all tested dimensions (n,m) in {(3,1), (4,2), (5,3), (6,2)}.

## What each term is

**Ambient acceleration** (fixed by the geometry of H):

    f''_amb = -Tr((H^{-1} Hdot)^2) + Tr(H^{-1} Hddot)

This depends only on H and its derivatives. No observer involved.

**Hidden acceleration** (fixed by the hidden sector):

    f''_hid = -Tr((R^{-1} U_h)^2) + Tr(R^{-1} Z^T Hddot Z)

where U_h = Z^T Hdot Z. This depends on H, its derivatives, and the hidden frame Z — but NOT on the observer C.

**Visible acceleration** (involves A_cpl):

    f''_vis = f''_amb - f''_hid

This is determined by the conservation law. It involves A_cpl through:

    Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})

## The conservation constraint on observer design

Given fixed ambient geometry (H, Hdot, Hddot) and hidden frame (Z), the visible acceleration is constrained:

    f''_vis = f''_amb - f''_hid

The observer C affects f''_vis by changing the split (Phi vs R, V vs U_h, etc.), but the total is fixed.

**Verified:** 20 random observers on the same (H, Hdot, Hddot):
- Conservation exact to 8.3e-17 for every observer
- Visible acceleration ranges from -0.012 to +0.068
- Hidden acceleration compensates exactly each time
- Sum always equals 0.0953 (the ambient acceleration)

## The acceleration budget

| (n,m) | Visible fraction | Hidden fraction |
|-------|-----------------|----------------|
| (3,1) | 68% | 32% |
| (4,2) | 98% | 2% |
| (5,3) | 75% | 25% |
| (6,2) | 45% | 55% |

The split varies with dimension and observer, but always sums to 100%.

## Interpretation

**The observer is not a source of information — only a splitter.**

The total information acceleration is fixed by the ambient geometry. The observer can only choose how to distribute it between visible and hidden sectors. An observer that maximises visible acceleration must minimise hidden acceleration, and vice versa.

This is the information-theoretic analogue of energy conservation in mechanics. The ambient geometry plays the role of total energy. The observer plays the role of a coordinate choice. You can make one component large by making another small, but you cannot create information from nothing.

This is also the structural reason why "every optimisation problem is an alignment problem with the underlying geometry": the geometry determines the total budget, and the observer's job is to align with it — to place the visible directions where the acceleration is most useful, not to create acceleration that isn't there.

## Connection to the full theory

The first-order conservation (R2) says: **the observer splits information.**
The second-order conservation says: **the observer splits the rate of change of information.**
A_cpl governs the visible share of the second-order split.

Together with the curvature-source coupling identity (T2), this gives the complete picture:
- The conservation law constrains A_cpl (it can't be arbitrary)
- The curvature F_alpha measures the cost of changing the split
- The hidden defect Q_hat measures how much the hidden sector affects the visible share
- The source law A_cpl governs the dynamics of the split

The observer's job is to navigate this budget optimally.
