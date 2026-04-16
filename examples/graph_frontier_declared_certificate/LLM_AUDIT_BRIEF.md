# LLM Audit Brief

## Mathematical Object

Audit a two-matrix weighted frontier where the declared observer is stationary
without being an exact branch.

## Commands

```bash
pip install -e .
python -m examples.graph_frontier_declared_certificate.run_all
```

## Expected Headline Outputs

- `outputs/frontier_sweep.csv`
- `outputs/hessian_breakpoint.svg`
- `outputs/summary.json`
- `outputs/validation.json`
- `outputs/audit.json`

## Acceptance Thresholds

- graph Hessian matches `-24 + 8 e^2` to `1e-10`
- stationarity residual is at most `1e-10`
- max certificate passes at `e = 1`
- both certificates are vacuous at `e = sqrt(3)`
- min certificate passes at `e = 2`
- `exact_branch_hessian` rejects the nonzero off-block cases

## Narrow Audit Question

Does the output support the narrow claim that the graph-frontier Hessian is the
right local quadratic object for stationary non-exact declared observers, while
the exact-branch Hessian must remain strict?
