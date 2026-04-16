# nomocomp Integration: Conservation Law as Model Comparison Diagnostic

**Date:** 15 April 2026
**Status:** Specification (no code yet)

---

## What nomocomp currently does

nomocomp extracts the observed Fisher information from fitted statsmodels OLS results and computes the exact geometric score: -2 log L + log det(I) - k log(2pi). It compares models by their geometric scores and detects ranking reversals vs AIC/BIC.

## What the conservation law adds

For two competing models M_1 and M_2 at the same data, the Fisher information changes: H_1 = I(theta_hat_1) and H_2 = I(theta_hat_2). The "perturbation" Hdot = H_2 - H_1 measures how the information geometry differs between models.

The conservation law says:

    vis_rate + hid_rate = amb_rate

where vis_rate is the information change captured by the observer's visible parameters, hid_rate is the change in the hidden (nuisance) parameters, and amb_rate is the total.

**New diagnostic:** the visible fraction vis_rate / amb_rate tells you how much of the model difference your comparison captures. If vis_frac is near 1, the comparison is seeing the full difference. If near 0, the difference is in parameters you can't see.

## Integration points

### 1. In `GeometricModelComparator.compare()`

After computing geometric scores for each model, also compute:

```python
from nomogeo import information_budget

# For each pair of models (M_i, M_j):
H_i = score_i.information_matrix
H_j = score_j.information_matrix
Hdot = H_j - H_i

# C = observation map (identity for full-parameter comparison)
# For subset comparison: C selects the shared parameters
budget = information_budget(H_i, C, Hdot)

result.information_budget = {
    'visible_rate': budget.visible_rate,
    'hidden_rate': budget.hidden_rate,
    'ambient_rate': budget.ambient_rate,
    'visible_fraction': budget.visible_fraction,
    'conservation_residual': budget.conservation_residual,
}
```

### 2. In `ComparisonResult`

Add a new field `information_budgets` that stores the conservation law diagnostics for each model pair. This tells the user:
- How much of the model difference is in the visible parameters
- Whether the comparison is "seeing" the full difference or missing hidden structure
- Whether the model difference is V > 0 (information-coherent) or V indefinite (information-mixed)

### 3. New diagnostic: model comparison stability

The positivity theorem says A_cpl >= 0 on linear paths (Hddot = 0). The model comparison H(t) = (1-t)H_1 + tH_2 IS a linear path. So:

- A_cpl >= 0 always holds for model comparisons
- The evidence is always concave along the model interpolation
- The geometric score ranking is stable within the hidden-defect radius

This is the "stability radius" of the model ranking: how much can the parameters change before the ranking flips? For nomocomp, this is:

    stability_radius = min_eigenvalue(A_cpl) / ||V^{-1/2}||^2

which is computable from the existing data.

## Implementation effort

Minimal. The `information_budget` function already exists in nomogeo. nomocomp already imports nomogeo. The integration is:
1. Call `information_budget(H_i, C, H_j - H_i)` for each model pair
2. Store the result in `ComparisonResult`
3. Add a `stability_radius` diagnostic

No new theory needed. The conservation law does the work.

## Expected impact

The information budget diagnostic answers a question AIC/BIC cannot: "is my model comparison seeing the full picture?" A visible fraction near 1 means yes. Near 0 means the models differ in hidden structure that the comparison is blind to. This is the same insight as the molecular result (45% of solvent response in soft modes) applied to model selection.
