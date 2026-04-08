# Claims

## Exact Setup

- covariance:
  `Sigma = [[A, C], [C^T, B]]`
- latent precision:
  `H = Sigma^{-1}`
- observer:
  `C_A = [[1,0,0,0],[0,1,0,0]]`
- visible precision:
  `Phi = (C_A H^{-1} C_A^T)^(-1)`
- ceiling:
  `T = H_AA`
- hidden load:
  `Lambda = hidden_load(T, Phi, support_mode="ambient")`
- clock:
  `tau = log det(I + Lambda)`

Independent comparison:

`I(A:B) = 0.5 * (log det A + log det B - log det Sigma)`

## Executable Statement

For the two-mode squeezed vacuum family and the tested two-mode squeezed thermal
family, the saved outputs support:

`tau = 2 I(A:B)`

within floating-point residual.

Additional executable checks:

- local symplectic invariance under the explicit transformation family used in
  `validate.py`:
  `local_rotation(0.23, -0.41) @ local_squeezer(0.37, -0.18)`
- additivity under independent composition through
  `hidden_contraction` and `load_from_hidden_contraction`
- sanity edge cases:
  - `r = 0` gives zero clock and zero mutual information
  - small `r` stays on the identity curve
  - the thermal family still satisfies the identity

## Limitations

- Gaussian families only
- no claim about non-Gaussian entanglement
- no claim about the full quantum-information landscape

## Reproduction

```bash
pip install -e .
python -m examples.entanglement_hidden_load.run_all
```

Inspect:

- `outputs/headline_summary.csv`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- max identity residual across validation sweep: `<= 1e-10`
- max local symplectic invariance residual: `<= 1e-10`
- composition clock residual: `<= 1e-10`
- contraction recovery residual: `<= 1e-10`

The release object passes if all four checks pass.
