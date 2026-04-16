# The Extraction Layer: Turning Data into (H, C, Hdot)

**Date:** 15 April 2026
**Status:** Architecture specification

---

## The gap

Every kernel function takes (H, C, Hdot). Users have data. The extraction step is currently:
- **nomocomp:** `extract_information(statsmodels_result)` → H. Hdot = H2 - H1.
- **nomoselect:** implicit. H = whitened within-class precision. Hdot = task family aggregate.
- **Molecules:** custom pipeline per data format.

The gap is: no unified extraction from raw data to (H, C, Hdot).

## The opportunity

A universal extraction function that takes any of:

### Supervised learning (X, y)
```python
from nomogeo import extract_geometry

# Fisher discrimination
H, Hdot = extract_geometry(X, y, task="fisher")
# H = regularised within-class precision
# Hdot = between-class scatter (= the task)

# Equal-weight discrimination  
H, Hdot = extract_geometry(X, y, task="equal_weight")

# Minority focus
H, Hdot = extract_geometry(X, y, task="minority")
```

### Model comparison (statsmodels)
```python
# Two fitted models
H, Hdot = extract_geometry(result_1, result_2)
# H = Fisher(model_1)
# Hdot = Fisher(model_2) - Fisher(model_1)
```

### Raw matrices
```python
# Direct
H, Hdot = extract_geometry(H_matrix, perturbation_matrix)
```

### Then everything follows
```python
from nomogeo import information_budget, observer_diagnostics, capture_curve

# Conservation law
budget = information_budget(H, C, Hdot)

# Observer quality
diag = observer_diagnostics(H, C, Hdot)

# How many dimensions?
curve = capture_curve(H, Hdot)

# The adapted observer
from nomogeo import closure_adapted_observer
result = closure_adapted_observer(H, [Hdot], m)
```

## What this enables

### One-line diagnosis
```python
from nomosteer import diagnose

# Complete diagnostic for any supervised learning problem
report = diagnose(X, y, task="fisher", m=3)
# report.budget: conservation law
# report.capture_curve: how many dimensions
# report.observer: adapted observer
# report.diagnostics: vis_frac, exact_sector, leakage
# report.comparison: vs PCA, vs LDA
```

### The kill shot against sklearn
```python
# sklearn way (no diagnostics, no certification)
from sklearn.decomposition import PCA
pca = PCA(n_components=3).fit_transform(X)

# nomosteer way (full diagnostics, certified)
from nomosteer import diagnose
report = diagnose(X, y, task="fisher", m=3)
X_proj = report.observer @ X.T  # project
print(report.summary())
# "Captures 95% of task structure. PCA captures 12%.
#  In exact sector. Conservation verified."
```

## Implementation plan

1. `extract_geometry(X, y, task)` — supervised data → (H, Hdot)
2. `extract_geometry(model_1, model_2)` — statsmodels → (H, Hdot)
3. `extract_geometry(H, Hdot)` — passthrough for raw matrices
4. `diagnose(data, task, m)` — full pipeline: extract → budget → curve → observer → diagnostics

The extraction layer is ~100 lines of code. The diagnosis pipeline is ~50 lines connecting existing functions. The value is in making 24 proved theorems accessible in one function call.
