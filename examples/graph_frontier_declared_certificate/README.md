# Graph Frontier Declared Certificate

Reproduce from a clean installed environment:

```bash
pip install -e .
python -m examples.graph_frontier_declared_certificate.run_all
```

Then inspect:

- [summary.json](outputs/summary.json)
- [audit.json](outputs/audit.json)
- [frontier_sweep.csv](outputs/frontier_sweep.csv)
- [hessian_breakpoint.svg](outputs/hessian_breakpoint.svg)

## Scope Limitation

This is a synthetic local quadratic observer-geometry example. It is not a
global observer optimizer, not a probability-support event classifier, and not
a generic non-Gaussian full-law branch selector.

## What Is Shown

For the symmetric pair

```text
A_+ = [[3, e], [ e, 1]]
A_- = [[3,-e], [-e, 1]]
```

with equal weights and observer `span(e_1)`, the observer is stationary for
every `e` even though it is not an exact branch when `e != 0`.

The graph-frontier Hessian is

```text
D^2 F(0) = -24 + 8 e^2
```

so the local sign changes at `|e| = sqrt(3)`.

The narrow message is:

use `general_graph_frontier_hessian` for declared non-exact stationary
observers; do not use an exact-branch proxy outside the exact-branch sector.

## What Is Not Claimed

- no global Grassmannian search
- no full-law branch probability
- no support/restart probability event
- no softening of `exact_branch_hessian`

## Public Outputs

- sweep table:
  [frontier_sweep.csv](outputs/frontier_sweep.csv)
- headline figure:
  [hessian_breakpoint.svg](outputs/hessian_breakpoint.svg)
- validation:
  [validation.json](outputs/validation.json)
- audit:
  [audit.json](outputs/audit.json)

See also:

- [CLAIMS.md](CLAIMS.md)
- [LLM_AUDIT_BRIEF.md](LLM_AUDIT_BRIEF.md)
- [RESULT.md](RESULT.md)
