# Static-Dynamic Unification

**Date:** 15 April 2026
**Code:** `0_4_0_static_dynamic_unification.py`
**Data:** 139 matched water-vacuum molecular pairs

---

## The question

Does the static observer optimum (from `closure_adapted_observer`, which maximises visibility at fixed H) coincide with the dynamic optimum (maximum visible rate under solvent perturbation)?

## Results

| Observer | Median vis_frac | V > 0 rate | Beats canonical |
|----------|-----------------|------------|-----------------|
| Canonical (soft modes) | 0.384 | 8.6% | -- |
| **Adapted (static)** | **1.349** | 21.6% | 87.1% |
| Optimised (dynamic) | 1.155 | 38.1% | 99.3% |

**Spearman rank correlation between adapted and optimised: 0.83 (p = 1.6e-36).**

## Interpretation

### The static observer is competitive with the dynamic optimum

The closure-adapted observer, designed purely from static visibility considerations, achieves a higher median visible fraction (1.35) than the random-search dynamic optimum (1.16). This means the static geometry (TN1) already encodes most of the dynamic information (0.4.0).

### The same geometry governs both

The Spearman correlation of 0.83 means the molecular ranking is nearly the same under both criteria. Molecules that are good for static observation are also good for dynamic observation. This is the unification result: **the static visibility frontier and the dynamic information budget are controlled by the same underlying split geometry.**

### The dynamic optimum has better V > 0 activation

The random-search optimised observer puts 38% of molecules into the exact sector (V > 0), vs 22% for the adapted observer. This makes sense: the adapted observer was designed for visibility, not for V > 0 specifically. An observer designed to maximise V's minimum eigenvalue would likely do even better.

### Both massively beat the canonical observer

Both the static and dynamic observers beat the canonical soft-mode observer on 87-99% of molecules. The canonical observer is a poor default for solvent-response analysis.

## What this means for the theory

The adapted observer from TN1 was designed using the weighted-family stationarity equation, which optimises a leakage-visibility trade-off. The dynamic visible fraction from 0.4.0 measures how much solvent-response information the observer captures. The fact that these two quantities have rank correlation 0.83 means:

1. The TN1 frontier is dynamically relevant (not just a static property)
2. The 0.4.0 conservation law is captured by the TN1 design criterion
3. A single observer design framework serves both the static and dynamic layers

This is the structural unification across the stack.
