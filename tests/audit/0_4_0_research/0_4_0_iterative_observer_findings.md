# Iterative Observer Refinement: Findings

**Date:** 15 April 2026

---

## Key discovery: the adapted observer already achieves B = 0

The adapted observer from `closure_adapted_observer(H, [Hdot], m)` finds the m-dimensional subspace that commutes with the whitened Hdot. In this subspace, the visible-hidden coupling B = L^T Hdot Z = 0 exactly.

Verified on n = 20 (synthetic) and n = 48 (synthetic): ||B|| = 0 to machine precision in both cases.

**Consequence:** The hidden defect Q_hat = B R^{-1} B^T = 0. The source law reduces to A_cpl = A_direct (no hidden correction). The observer is in the zero-coupling sector.

## Why iteration doesn't help

Gradient descent on ||B||^2 doesn't move because the adapted observer already sits at ||B|| = 0, a global minimum. There's nothing to iterate toward — the adapted observer found it in one shot.

## Why the Grassmannian optimiser finds something different

The Grassmannian optimiser maximises vis_frac = vis_rate / amb_rate, not ||B||^2. These are different objectives:

- **Adapted:** minimises ||B|| → achieves B = 0 → Q_hat = 0 → source law simplest
- **Grassmann(vf):** maximises vis/amb → may accept nonzero B if it improves the ratio

When amb_rate is large and positive, B = 0 also gives high vis_frac (because vis_rate ≈ amb_rate when there's no hidden coupling). But when amb_rate is small or negative, the B = 0 observer can give pathological vis_frac even though the coupling is minimal.

## The resolution

The adapted observer is the right answer for a different question: "find the observer with minimal hidden coupling." This is the exact-sector answer — A_cpl = A_direct, no hidden mediation, cleanest source law.

The Grassmannian optimiser answers a different question: "find the observer that captures the most information." This involves the full nonlinear chain and can accept hidden coupling if it helps the ratio.

For nomosteer:
- **If the goal is exact-sector entry:** the adapted observer is already optimal (B = 0)
- **If the goal is maximum vis_frac:** Grassmannian optimisation is needed, but only when the adapted observer gives pathological vis_frac (e.g., vis_frac < 0)
- **The adapted observer is the natural Tier 1** because it achieves the cleanest structural result (zero coupling) in milliseconds

## The right iteration

If iteration is needed, it should optimise vis_rate (not ||B||^2), using the adapted observer as warm start. But this is equivalent to Grassmannian optimisation with a warm start — there's no faster iterative scheme because the adapted observer has already found the structural optimum.
