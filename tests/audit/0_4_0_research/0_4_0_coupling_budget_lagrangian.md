# The Coupling Budget Lagrangian

**Date:** 15 April 2026
**Status:** Theoretical sketch (to be formalised)

---

## The setup

The conservation law constrains the observer:

    f''_vis + f''_hid = f''_amb    (fixed by geometry)

The curvature F_alpha measures the cost of changing the observer:

    ||F||^2 = 2(Tr(G_1 G_2) - Tr(C^2))    (Grassmannian curvature)

The source law A_cpl governs the visible evidence acceleration:

    Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})

The question: what is the optimal observer trajectory C(t) as the law H(t) changes?

## The variational problem

Given a path H(t) through SPD(n), find the observer path C(t) through Gr(m, n) that optimises a trade-off between:

1. **Information capture:** maximise the visible rate d/dt[log det Phi]
2. **Observer stability:** minimise the curvature ||F_alpha||^2 of the observer path
3. **Exact-sector preference:** maximise the minimum eigenvalue of V (support stability)

Subject to:
- The conservation constraint: vis_rate + hid_rate = amb_rate at all t
- The split-frame identities: M^T H M = diag(Phi, R)

## The Lagrangian

Define the functional on observer paths C: [0,T] -> Gr(m,n):

    L[C] = integral_0^T [ lambda_1 * vis_rate(t) - lambda_2 * ||dC/dt||^2_Gr - lambda_3 * max(0, -V_min(t)) ] dt

where:
- lambda_1 weights the information capture objective
- lambda_2 weights the observer motion cost (curvature penalty)
- lambda_3 weights the exact-sector preference (V > 0 penalty)
- ||dC/dt||^2_Gr is the Grassmannian speed (related to ||beta||^2)

The conservation constraint is automatic (it holds for any C).

## The Euler-Lagrange equation

For the simplified case lambda_3 = 0 (no V > 0 penalty):

    lambda_1 * d/dC [vis_rate] = lambda_2 * d/dC [||dC/dt||^2_Gr]

The LHS is the gradient of the visible rate with respect to the observer. From the conservation law, this equals -d/dC [hid_rate] (since amb_rate doesn't depend on C).

The RHS is the Riemannian acceleration on the Grassmannian. When this is zero (geodesic), the observer moves at constant speed.

**At a stationary observer (dC/dt = 0):** the Euler-Lagrange reduces to d/dC [vis_rate] = 0, which is the TN1 weighted-family stationarity equation evaluated for the perturbation Hdot instead of a static family.

This confirms the unification: **the TN1 frontier is the zero-velocity limit of the coupling budget Lagrangian.**

## Connection to the empirical results

The Spearman correlation of 0.83 between static and dynamic optima means the TN1 solution (zero-velocity limit) is already near the dynamic optimum for most molecules. The Lagrangian explains why: when the observer doesn't need to move (the law change is captured by the static design), the geodesic solution is to stay put. The static optimum IS the dynamic optimum in the slow-observer limit.

The 17% of molecules where adapted and optimised diverge are the cases where the observer needs to move — where the Grassmannian curvature term matters and the static design is no longer sufficient.

## What needs to be formalised

1. The Grassmannian metric on Gr(m, n) in terms of the split-frame connection forms
2. The gradient d/dC [vis_rate] in terms of beta and theta
3. The Euler-Lagrange equation with the exact-sector penalty
4. The geodesic equation on Gr(m, n) with the H-induced metric
5. The stability analysis around the static solution

This is a well-posed variational problem. The ingredients are all in the 0.4.0 TN. The Lagrangian connects the TN1 static frontier (zero-velocity limit) to the full dynamic observer trajectory.

## Sketch of the gradient

From the conservation law: vis_rate = amb_rate - hid_rate.
amb_rate = Tr(H^{-1} Hdot) does not depend on C.
hid_rate = Tr(R^{-1} U_h) where R = Z^T H Z and U_h = Z^T Hdot Z.

So: d/dC [vis_rate] = -d/dC [hid_rate].

Under a variation delta C (which changes Z to Z + delta Z):
- R changes by delta R = (delta Z)^T H Z + Z^T H (delta Z)
- U_h changes by delta U_h = (delta Z)^T Hdot Z + Z^T Hdot (delta Z)
- hid_rate changes by Tr(R^{-1} delta U_h) - Tr(R^{-1} delta R R^{-1} U_h)

This is computable from the existing split-frame data. The gradient of the visible rate IS the gradient of the Lagrangian in the zero-velocity limit.
