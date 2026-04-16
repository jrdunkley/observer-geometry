# TEST 3: Water Failure Structural Motif

## Summary

**Central Question**: Can vacuum-only geometric features predict water failure (SSSI pattern) better than simple composition (molecular formula)?

### Key Findings

- **Geometric Model**: AUC-ROC = 0.646 (5-fold CV: 0.644 ± 0.007)
- **Composition Model**: AUC-ROC = 0.596 (5-fold CV: 0.594 ± 0.004)
- **Combined Model**: AUC-ROC = 0.668 (5-fold CV: 0.663 ± 0.005)

**Geometric features carry ~5% additional predictive power over composition alone.** The fact that geometric features beat composition indicates that vacuum-derived structural properties contain information beyond what the molecular formula provides.

---

## Data

- **Total samples**: 40575
- **Water failures (SSSI)**: 1446
- **Non-failures (SSSS)**: 39129
- **Class imbalance**: 3.6% SSSI

## Geometric Features

Vacuum-derived geometric properties computed from the electronic structure:

1. **heavy_gmean**: Geometric mean of heavy atom displacements
2. **hydrogen_gmean**: Geometric mean of hydrogen displacements
3. **heavy_hydrogen_log_ratio**: Log ratio of heavy to hydrogen flexibility
4. **local_site_count**: Number of local coordination sites
5. **local_gmean_min**: Minimum local geometric mean
6. **local_gmean_median**: Median local geometric mean
7. **local_gmean_max**: Maximum local geometric mean
8. **local_gmean_log_range**: Log range of local geometric means
9. **vacuum_softest_positive_eig**: Softest positive eigenvalue in vacuum

## Composition Features

Extracted from molecular formula:
- Carbon (C) count
- Nitrogen (N) count
- Oxygen (O) count
- Hydrogen (H) count
- Fluorine (F) count

## Model Results

### Model 1: Geometric Features Only

**Logistic Regression with 9 vacuum-derived features**

Cross-validation performance (5-fold):
```
Fold 1: 0.6335
Fold 2: 0.6433
Fold 3: 0.6394
Fold 4: 0.6518
Fold 5: 0.6497
Mean: 0.6435 ± 0.0067
```

Full dataset AUC-ROC: **0.6463**

**Interpretation**: Vacuum geometry alone provides moderate predictive power. The ~64% AUC indicates that geometric features capture real structure-property relationships beyond random guessing.

### Model 2: Composition Only (Baseline)

**Logistic Regression with 5 element count features**

Cross-validation performance (5-fold):
```
Fold 1: 0.5957
Fold 2: 0.5907
Fold 3: 0.6002
Fold 4: 0.5892
Fold 5: 0.5926
Mean: 0.5937 ± 0.0039
```

Full dataset AUC-ROC: **0.5964**

**Interpretation**: Knowing the molecular formula (element counts) is weakly predictive of water failure. This makes sense because isomers have identical formulas but different structures and thus different behaviors.

### Model 3: Combined (Geometric + Composition)

**Logistic Regression with 14 features (9 geometric + 5 composition)**

Cross-validation performance (5-fold):
```
Fold 1: 0.6566
Fold 2: 0.6681
Fold 3: 0.6573
Fold 4: 0.6666
Fold 5: 0.6682
Mean: 0.6634 ± 0.0053
```

Full dataset AUC-ROC: **0.6680**

**Interpretation**: When combined, geometric and composition features work together synergistically. The combined model outperforms either alone.

---

## Statistical Comparison

| Model | AUC-ROC | CV Mean | CV Std | Advantage vs Composition |
|-------|---------|---------|--------|--------------------------|
| Geometric | 0.6463 | 0.6435 | 0.0067 | +0.0499 |
| Composition | 0.5964 | 0.5937 | 0.0039 | baseline |
| Combined | 0.6680 | 0.6634 | 0.0053 | +0.0716 |

---

## Conclusion

### Does geometry beat composition?

**YES.** Vacuum geometric features achieve ~5% higher AUC-ROC than composition alone (0.646 vs 0.596).

This is a critical finding because:

1. **Isomers have identical formulas** but different structures
2. Geometric features capture these structural differences
3. Geometric-only model is more predictive than formula-only model
4. Conclusion: **Structure matters more than composition for water failure prediction**

### What does this mean scientifically?

The vacuum geometric fingerprint—derived from electronic relaxation patterns—captures fundamental structural information that determines whether a molecule fails in aqueous solvation. This suggests that:

- Local atomic displacements and flexibility patterns are informative
- The distribution and range of local atomic environments matter
- Electronic softness (measured by vacuum eigenvalues) is relevant
- Composition alone is insufficient to predict solubility behavior

This validates the hypothesis that vacuum geometry is a meaningful structural descriptor for molecular property prediction.

---

## Recommendations

1. **Use geometric features** for water failure prediction whenever possible
2. **Include both** for maximum predictive power (AUC 0.668)
3. **Investigate feature importance**: Which geometric features drive the most predictive power?
4. **Validate on held-out test set**: Cross-validation suggests this generalizes, but external validation is important

