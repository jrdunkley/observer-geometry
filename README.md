# nomogeo

Exact observer-geometry kernel for visible precision, hidden-load calculus, information conservation, and observer steering.

Current release: **nomogeo 0.4.0**.

## What nomogeo does

Given a symmetric positive-definite matrix H (the ambient precision / Fisher information / Hessian) and a rank-m observer C, nomogeo computes exact split-frame geometry:

- **Visible precision** Phi = (C H^{-1} C^T)^{-1} and canonical lift L
- **Information budget**: vis_rate + hid_rate = amb_rate (exact conservation)
- **Source law**: A_cpl on support-stable strata, with hidden-defect decomposition
- **Observer steering**: adapted observer (B = 0), capture curves, Grassmannian diagnostics
- **Geometry extraction**: from supervised data (X, y) or covariance matrices directly to (H, Hdot)

## Install

```bash
pip install nomogeo
```

Runtime dependencies: `numpy >= 2.0`, `scipy >= 1.15`. No other requirements.

## Quick start

```python
import numpy as np
from nomogeo import visible_precision, information_budget, steer

# Exact visible precision
H = np.array([[3.0, 1.0], [1.0, 2.0]])
C = np.array([[1.0, 0.0]])
phi = visible_precision(H, C)

# Information conservation: vis + hid = amb
Hdot = np.array([[0.1, -0.2], [-0.2, 0.3]])
budget = information_budget(H, C, Hdot)
print(f"vis={budget.visible_rate:.4f}, hid={budget.hidden_rate:.4f}, "
      f"amb={budget.ambient_rate:.4f}")
print(f"conservation residual: {budget.conservation_residual:.1e}")
```

### Observer steering from data

```python
from sklearn.datasets import load_iris
from nomogeo import steer

X, y = load_iris(return_X_y=True)
result = steer(X=X, y=y, rank=2, task="fisher")
print(f"vis_frac={result.visible_fraction:.3f}, "
      f"advantage over PCA: {result.advantage_over_pca:.3f}")
print(f"exact sector: {result.exact_sector}")
```

## Core API

### Geometry

| Function | Description |
|----------|-------------|
| `visible_precision(H, C)` | Phi_C(H) = (C H^{-1} C^T)^{-1} |
| `canonical_lift(H, C)` | L = H^{-1} C^T Phi |
| `hidden_projector(H, C)` | P = I - L C |
| `fixed_observer_coordinates(H, C)` | Exact (Phi, R, K) chart |
| `observer_transition(H, C1, C2)` | Exact transition between observers |

### Information conservation and source law

| Function | Description |
|----------|-------------|
| `information_budget(H, C, Hdot)` | vis + hid = amb (first order) |
| `source_law(H, C, Hdot, Hddot)` | A_cpl on support-stable strata |
| `evidence_decomposition(H, C, Hdot, Hddot)` | Second-order: source + kinematic + connection |
| `observer_diagnostics(H, C, Hdot)` | vis_frac, exact sector, hidden defect, leakage |
| `capture_curve(H, Hdot)` | Information capture vs observer rank |

### Extraction and steering

| Function | Description |
|----------|-------------|
| `extract_supervised(X, y, task)` | (X, y) -> (H, Hdot) for fisher/equal_weight/minority |
| `extract_covariance(cov, perturbation)` | Raw matrices -> (H, Hdot) |
| `steer(X, y, rank, task)` | One-call: extract + optimise + diagnose |
| `score_observer(H, C, Hdot)` | Score an observer against baselines |
| `optimize_observer(H, Hdot, rank)` | Adapted observer + PCA fallback |

### Observation fields and regime classification

| Function | Description |
|----------|-------------|
| `support_stratum_transport(Lambda, A_cpl)` | Hidden-load transport diagnostics |
| `local_coupled_birth(H, Hdot, Hddot, C, Cdot)` | Birth/death events at support transitions |
| `classify_regime(datum)` | Quadratic regime classification |

## 0.4.0 Technical Note

The accompanying technical note contains 30+ formal statements:

- **Conservation laws**: vis_rate + hid_rate = amb_rate (first and second order)
- **Source law**: A_cpl = A_direct + hidden_defect, with completed-square decomposition
- **Positivity theorem**: A_cpl >= 0 on linear paths (hidden sector stabilises)
- **Equivalence principle**: A_cpl can always be locally cancelled
- **Curvature-Gram identity**: ||F_alpha||^2 = 2(Tr(G1 G2) - Tr(C^2))
- **Mixed factorisation**: F_alpha(ds, dt) = -beta_t R^{-1} B_s^T
- **Gauge invariance**: vis_rate is a gauge-invariant Noether current
- **Spectral collapse lemma**: no scalar constraint makes vis_rate maximisation well-posed
- **KL field equation**: S - Xi = (gamma/2)(H - H0) under KL-divergence regularisation
- **Joint optimum theorem**: closed-form solution with frozen hidden sector

Validated on 118,539 molecular pairs from HessianQM9 and 108/108 field equation checks.

## Companion packages

- **nomocomp** (0.1.0): geometric model comparison replacing AIC/BIC with exact fibre-volume correction
- **nomoselect** (0.1.0): task-aware subspace selection replacing PCA with certified observer design

Both depend on nomogeo >= 0.4.0.

## Architecture

```
nomogeo/
  src/nomogeo/
    core.py          -- visible precision, lift, projector
    connection.py    -- fixed-observer charts, transitions, currents
    source.py        -- conservation, source law, diagnostics, capture curves
    extract.py       -- data -> geometry extraction
    steer.py         -- observer steering engine
    field.py         -- observation fields, birth/death, support transport
    regime.py        -- quadratic regime classification
    singular.py      -- singularity analysis
    kernel_reduction.py -- affine-hidden Gaussian-fibre reduction
    frontier.py      -- weighted-family frontier evaluation
    ensemble.py      -- local quadratic ensemble diagnostics
    adapted.py       -- closure-adapted observer synthesis
    affine.py        -- thin affine layers
  tests/             -- 340+ tests
    audit/           -- complete research audit trail
  examples/          -- 6+ runnable demonstrations
  papers/            -- LaTeX source for technical notes
```

## Verification

```bash
PYTHONPATH=src python -m pytest tests/ -q
```

340+ tests, all passing. 2 skipped (require external datasets).

## Exact domain

nomogeo is exact for:
- Linear observers on finite-dimensional SPD matrices
- Gaussian/quadratic visible objects
- Matrix identities where theorems apply
- Explicitly supplied special-law sectors

It is not a generic non-Gaussian engine. Outside the exact Gaussian sector, it provides local quadratic Hessian/Fisher geometry plus separately proved special sectors.

## License

BSD 3-Clause. See LICENSE
