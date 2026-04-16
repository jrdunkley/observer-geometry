# Source Law on HessianQM9 Full Hessian Data

**Date:** 15 April 2026
**Code:** `0_4_0_qm9_hessian_source_law.py`
**Data:** HessianQM9 Arrow (vacuum + THF, 200 molecules per shard, 5 shards)

---

## Results

### Conservation on 17 matched vacuum-THF pairs

The first-order conservation law was verified on 17 molecules appearing in both vacuum and THF environments, with full vibrational Hessians up to 57x57 (19-atom molecules).

**Maximum conservation error: 6.4e-14.** The law holds at machine precision on real molecular Hessian data with full off-diagonal structure.

### Visible fraction distribution

| Statistic | Value |
|-----------|-------|
| Molecules analysed | 17 |
| V > 0 (source law active) | 0/17 |
| Mean visible rate | 5.91 |
| Mean hidden rate | 2.36 |
| Mean visible fraction | 94% |
| Visible fraction range | -53% to 99.5% |

### Honest boundary: V > 0 is restrictive

The source law's support condition V > 0 is never satisfied for vacuum-to-THF solvent perturbation on real molecules. Solvent changes generically mix positive and negative eigenvalue changes in the visible jet V = L^T (H_THF - H_vac) L.

This means A_cpl cannot be computed in the standard way for these paths. However:
- The conservation law vis_rate + hid_rate = ambient_rate operates unconditionally
- The visible/hidden rate split is the primary diagnostic
- A_cpl would become accessible on restricted paths where V > 0 (e.g., perturbation within a single solvent by varying concentration)

### Physical interpretation

The visible fraction varies from -53% to 99.5% across the 17 molecules. This means:
- Some molecules have nearly all their solvent response in the soft (observed) modes — the observer is well-aligned with the perturbation
- Others have their response mostly in higher modes — the observer misses the action
- Some have negative visible fraction — the observer is actively misaligned, seeing the opposite of what the perturbation does

This is the information budget of molecular spectroscopy: which vibrational modes to observe depends on which modes carry the solvent response, and the conservation law quantifies the trade-off exactly.
