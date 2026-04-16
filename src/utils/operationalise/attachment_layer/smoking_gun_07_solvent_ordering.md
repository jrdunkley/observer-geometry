# TEST 7: SOLVENT ORDERING CONSISTENCY

## Executive Summary

The three solvents (THF, toluene, water) do **NOT** produce a perfectly consistent geometric ordering across all molecules. However, the expected weak ordering (water > toluene in perturbation strength) is obeyed in **85.79% of cases**, and the full expected ordering (water > THF > toluene) holds in **51.99% of cases**. Strong rank correlations between all three solvent shifts suggest they operate along a shared "solvation susceptibility" axis, but independent information still exists.

---

## Detailed Results

### 1. Absolute Value Shifts - Summary Statistics

| Solvent | Mean |shift| | Std Dev | Interpretation |
|---------|---------|---------|---------|
| **toluene** (least polar) | 0.7990 | 0.5199 | Expected to perturb LEAST |
| **THF** (moderately polar) | 1.0238 | 0.7322 | Expected INTERMEDIATE |
| **water** (most polar) | 1.2792 | 0.8316 | Expected to perturb MOST |

**Key observation**: Mean shifts increase with solvent polarity as predicted by chemistry. The ordering holds on average, but with substantial variance.

---

### 2. Expected Ordering Compliance

#### Full Ordering: |water| > |THF| > |toluene|
- **Compliance rate: 51.99%** (21,100 / 40,583 molecules)
- Only about half of molecules follow the full expected hierarchy

#### Weak Ordering: |water| > |toluene|
- **Compliance rate: 85.79%** (34,818 / 40,583 molecules)
- Water does perturb more than toluene in most cases
- **Violation rate: 14.21%** (5,765 molecules)

**Conclusion**: The water > toluene ordering is robust. The water > THF > toluene full hierarchy is less consistent.

---

### 3. Detailed Ordering Breakdown

| Ordering Type | Count | Percentage |
|---------|---------|---------|
| **w > t > tol** (EXPECTED) | 21,100 | 51.99% |
| w > tol > t | 7,649 | 18.85% |
| t > w > tol | 6,069 | 14.95% |
| t > tol > w | 2,768 | 6.82% |
| tol > t > w | 2,058 | 5.07% |
| tol > w > t | 939 | 2.31% |

**Interpretation**: 
- The expected ordering dominates but is not exclusive
- 18.85% of molecules show THF and toluene interchanged (w > tol > t), suggesting these two solvents can have variable relative effects
- ~8% show complete reversals where toluene perturbs more than water
- No single ordering dominates to the point of being a universal law

---

### 4. Violations: Toluene Perturbs More Than Water

**Number of violations: 5,765 (14.21% of dataset)**

This is the most chemically surprising result: the least polar solvent (toluene) sometimes perturbs the geometry MORE than the most polar solvent (water).

#### Top 10 Most Extreme Violations

| Molecule | Formula | |toluene| | |water| | Ratio |
|---------|---------|---------|---------|---------|
| dsgdb9nsd_081582 | C9H14 | 0.6177 | 0.0002 | **3,088x** |
| dsgdb9nsd_092016 | C7H12N2 | 0.3078 | 0.0001 | **3,078x** |
| dsgdb9nsd_085941 | C7H14O2 | 0.4939 | 0.0004 | **1,235x** |
| dsgdb9nsd_081987 | C9H12 | 0.3490 | 0.0003 | **1,163x** |
| dsgdb9nsd_083678 | C7H9NO | 1.3140 | 0.0014 | **938x** |
| dsgdb9nsd_058885 | C7H11NO | 1.0701 | 0.0011 | **973x** |
| dsgdb9nsd_030711 | C6H10N2O | 1.0279 | 0.0015 | **685x** |
| dsgdb9nsd_071199 | C9H12 | 0.1846 | 0.0003 | **615x** |
| dsgdb9nsd_070905 | C8H12O | 0.0528 | 0.0001 | **528x** |
| dsgdb9nsd_101410 | C6H10O3 | 0.8884 | 0.0023 | **386x** |

**Chemical insight**: 
- Many violations involve small, nonpolar or weakly polar molecules (e.g., C9H14, C7H12N2, C9H12)
- These molecules appear to have geometric features that toluene can perturb efficiently (possibly π-system alignment, hydrophobic interactions) while water has minimal effect
- Suggests the perturbation mechanism is NOT simply polarity-driven but depends on molecule-solvent compatibility

---

### 5. Rank Correlations - Testing for a Single "Solvation Susceptibility" Axis

**Hypothesis**: If all three solvents perturb the same underlying geometric property, their shifts should be highly correlated.

#### Spearman Rank Correlations (nonlinear relationships)

| Pair | Spearman ρ | p-value | Interpretation |
|---------|---------|---------|---------|
| water vs THF | 0.7351 | <0.001 | **STRONG** |
| water vs toluene | 0.7141 | <0.001 | **STRONG** |
| THF vs toluene | 0.8383 | <0.001 | **VERY STRONG** |

#### Pearson Correlations (linear relationships)

| Pair | Pearson r | Interpretation |
|---------|---------|---------|
| water vs THF | 0.6362 | **MODERATE-STRONG** |
| water vs toluene | 0.6008 | **MODERATE-STRONG** |
| THF vs toluene | 0.7180 | **STRONG** |

**Critical finding**: All correlations are statistically significant and in the 0.60-0.84 range, indicating:
1. **Shared axis exists**: Molecules with high water-shift tend to have high THF and toluene shifts
2. **Independent information remains**: Correlations are strong but not perfect (r < 0.85)
3. **THF-toluene coupling is strongest**: These two solvents perturb molecules in remarkably similar ways (ρ = 0.8383)

---

## Chemical Interpretation

### What's Working
✓ **Polarity-driven effects are real but soft**: On average, water > THF > toluene holds  
✓ **Water-toluene ordering is robust**: 85.79% of molecules obey this  
✓ **Shared solvation axis detected**: All three solvents couple to a common geometric property  

### What's Surprising
✗ **Full ordering only 52% reliable**: THF and toluene can interchange effects  
✗ **14% of molecules violate water > toluene**: Some nonpolar molecules respond much more to toluene  
✗ **Polarity isn't the whole story**: Solvent-molecule compatibility appears highly molecular-specific  

### Data Quality Implications
- **The ordering consistency is moderate, not absolute**: This suggests the solvent perturbation mechanism is complex and molecule-dependent
- **Strong correlations indicate reliable relative rankings**: Molecules can be reliably ranked by "solvation susceptibility," but not necessarily in a fixed polarity order
- **The 14% violation rate is concerning but not critical**: If this were random noise, we'd see much weaker correlations

### Recommendations
1. **Accept the solvent shifts as valid**: Strong correlations and mean-level ordering support data integrity
2. **Treat the full ordering (water > THF > toluene) as a statistical tendency, not a law**: ~52% compliance is high enough to be meaningful but low enough to expect exceptions
3. **Investigate the 14% violators**: These molecules may reveal new geometric or chemical features
4. **Use rank-based aggregation**: Since correlations are strongest for rank (Spearman > Pearson), treat molecules' relative solvation susceptibility rankings as more reliable than absolute shifts

---

## Conclusion

**The solvent data is reasonably reliable** for detecting geometric perturbations. The expected polarity-based ordering holds on average and produces strong correlations. However, **the ordering is not universally obeyed** — molecules show significant variation in how different solvents affect them. This likely reflects real molecular properties (geometry, polarity distribution, π-systems) rather than data noise, making the dataset suitable for downstream analysis while remaining alert to outliers.
