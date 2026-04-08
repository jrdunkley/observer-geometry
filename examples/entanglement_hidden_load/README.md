# Entanglement As Hidden Load

Reproduce from a clean installed environment:

```bash
pip install -e .
python -m examples.entanglement_hidden_load.run_all
```

Then inspect:

- [summary.json](outputs/summary.json)
- [audit.json](outputs/audit.json)

## What Is Shown

For the explicit Gaussian families in this example, the `nomogeo`
hidden-load clock satisfies

```text
tau = 2 I(A:B)
```

where

```text
I(A:B) = 0.5 * (log det A + log det B - log det Sigma)
```

is computed independently from the Gaussian covariance matrix.

The release object includes:

- a two-mode squeezed vacuum sweep
- a two-mode squeezed thermal sweep
- a local symplectic invariance check under
  `local_rotation(0.23, -0.41) @ local_squeezer(0.37, -0.18)`
- an additive composition check carried through the associative contraction
  object, not raw load coordinates

## What Is Not Claimed

- this is not a universal theorem about all quantum states
- this does not classify non-Gaussian entanglement
- this does not use non-Gaussian measurements or entanglement witnesses

## Why This Is Difficult To Dismiss

- the identity is checked against an independent standard Gaussian mutual
  information formula
- the result is a sweep, not a single point
- the local symplectic change is explicit and the residual is reported
- additivity is checked through the correct associative contraction object

## Exact Setup

- latent precision: `H = Sigma^{-1}`
- observer:
  `C_A = [[1,0,0,0],[0,1,0,0]]`
- visible object:
  `Phi = visible_precision(H, C_A)`
- ceiling:
  `T = H_AA`
- hidden object:
  `Lambda = hidden_load(T, Phi, support_mode="ambient")`
- long composition is always carried through
  `hidden_contraction(Lambda)` and `load_from_hidden_contraction(K)`, not raw
  repeated use of visible load coordinates

## Public Outputs

- headline figure:
  [tmsv_tau_vs_2mi.svg](outputs/tmsv_tau_vs_2mi.svg)
- residual figure:
  [tmsv_residual.svg](outputs/tmsv_residual.svg)
- summary table:
  [headline_summary.csv](outputs/headline_summary.csv)

See also:

- [CLAIMS.md](CLAIMS.md)
- [LLM_AUDIT_BRIEF.md](LLM_AUDIT_BRIEF.md)
- [RESULT.md](RESULT.md)

