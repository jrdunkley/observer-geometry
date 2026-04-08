# LLM Audit Brief

## Mathematical Object

Compute the `nomogeo` hidden-load clock beneath the local precision ceiling of
a bipartite Gaussian state and compare it to twice the standard Gaussian mutual
information.

## Commands

```bash
pip install -e .
python -m examples.entanglement_hidden_load.run_all
```

## Expected Headline Outputs

- `outputs/tmsv_identity.csv`
- `outputs/thermal_identity.csv`
- `outputs/headline_summary.csv`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- `max_abs_residual_tau_minus_2mi <= 1e-10`
- `max_local_symplectic_invariance_residual <= 1e-10`
- `composition_clock_residual <= 1e-10`
- `contraction_recovery_residual <= 1e-10`

## Narrow Audit Question

Does the output support the claim that `tau = 2 I(A:B)` across the stated
Gaussian sweeps within floating-point residual, while also passing the explicit
invariance and additivity checks?
