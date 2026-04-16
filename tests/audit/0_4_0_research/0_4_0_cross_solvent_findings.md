# N4: Cross-Solvent Stability

**Date:** 15 April 2026
**Data:** 37,798 molecules appearing in all 3 solvents (from 118K sweep)

---

## The question

Is the optimal observer for water-perturbation also good for THF and toluene?

## Answer: No. The optimal observer is solvent-specific.

### Cross-solvent correlations

| Pair | Spearman correlation |
|------|---------------------|
| THF vs Water | 0.257 |
| THF vs Toluene | 0.253 |
| Water vs Toluene | 0.252 |

Compare with the static-dynamic correlation of 0.83 (within the same perturbation). The cross-solvent correlation is 3x weaker.

### Universal "good" molecules are rare

| Category | Count | Fraction |
|----------|-------|----------|
| Good in all 3 (vf > 0.5) | 7,581 | 20.1% |
| Bad in all 3 (vf < 0.2) | 3,372 | 9.0% |

Only 20% of molecules have vf > 0.5 for all three solvents. The canonical soft-mode observer is a poor universal choice.

### Mean absolute difference

- |vf_THF - vf_water| = 1.20
- |vf_THF - vf_toluene| = 1.28

The visible fractions differ by more than 1.0 on average between solvents. The information budget is strongly perturbation-dependent.

## Significance

**The observer design must account for the expected perturbation.** The adapted observer (designed for Hdot = H_solvent - H_vacuum) will only be optimal for that specific solvent. For a different solvent, a different observer may be needed.

This has direct implications:
1. **nomoselect:** the task family should include the expected perturbation direction, not just the static structure
2. **nomocomp:** model comparison observers should be tailored to the specific model pair being compared
3. **Molecular spectroscopy:** the modes to observe depend on which solvent effect you're trying to measure

The static-dynamic unification (Spearman 0.83) holds within a single perturbation. Across perturbations, the geometry changes and the observer must adapt.
