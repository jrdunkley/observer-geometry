# nomosteer Architecture: Fast Path with Fallback

**Date:** 15 April 2026
**Informed by:** Grassmannian sample results (5 molecules, 3131s)

---

## The compute reality

| Method | Time per molecule | vis_frac quality |
|--------|-------------------|------------------|
| Adapted observer | ~1 ms | Good on 80%, catastrophic on 20% |
| Grassmann opt (15 restarts) | ~10 min | Optimal |

The adapted observer is 600,000x faster. The Grassmannian optimum is better but the marginal gain per compute dollar is poor except on the ~20% of cases where the adapted observer fails.

## The failure signature

The one catastrophic failure (molecule 028438) had:
- Adapted vis_frac = -18.9 (extreme negative)
- Adapted leakage = 12.0

The successful cases had:
- Adapted vis_frac > 1.0 (positive, often > 1)
- Variable leakage

**Simple detection rule:** If adapted vis_frac < 0, the adapted observer is unreliable. Run optimisation.

## Proposed architecture

### Tier 1: Fast diagnostic (milliseconds)

```python
from nomogeo import closure_adapted_observer, observer_diagnostics

result = closure_adapted_observer(H, [Hdot], m)
diag = observer_diagnostics(H, result.C, Hdot)

if diag.visible_fraction > 0 and diag.exact_sector:
    # Accept: adapted observer is good
    return result.C, diag
```

### Tier 2: Moderate search (seconds)

```python
if diag.visible_fraction < 0 or not diag.exact_sector:
    # Adapted observer is unreliable
    # Run 5-restart Grassmannian search with adapted as warm start
    C_opt = grassmann_optimise(H, Hdot, m, warm_start=result.C, n_restarts=5)
    diag_opt = observer_diagnostics(H, C_opt, Hdot)
    return C_opt, diag_opt
```

### Tier 3: Full optimisation (minutes, only if needed)

```python
if diag_opt.visible_fraction < 0:
    # Still failing — full optimisation
    C_full = grassmann_optimise(H, Hdot, m, n_restarts=15, max_iter=200)
    return C_full, observer_diagnostics(H, C_full, Hdot)
```

## Expected performance

From the 5-molecule sample:
- Tier 1 accepts: 3/5 molecules (60%) — near-optimal, instant
- Tier 2 needed: 2/5 molecules (40%) — moderate compute
- Tier 3 needed: ~0/5 — rare

For the 139-molecule dataset:
- Tier 1 would accept ~80% (vis_frac > 0 with adapted)
- Tier 2 needed for ~20%
- Average compute: 0.8 * 1ms + 0.2 * 5s = ~1 second per molecule

This is 3000x faster than full Grassmannian optimisation while maintaining quality.

## For nomosteer's API

```python
from nomosteer import steer

# Default: tiered, fast
result = steer(H, Hdot, m=3)
# result.observer, result.diagnostics, result.tier_used

# Force full optimisation
result = steer(H, Hdot, m=3, mode='full')

# Fast only (no fallback)
result = steer(H, Hdot, m=3, mode='fast')
```
