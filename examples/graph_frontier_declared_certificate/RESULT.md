# Result

Headline finding:

- the declared observer is stationary for the symmetric pair even when each
  family member has nonzero off-block coupling
- the graph-frontier Hessian matches `-24 + 8 e^2`
- the local sign changes at `|e| = sqrt(3)`
- the max certificate passes before the threshold, the degenerate case is
  vacuous, and the min certificate passes after the threshold
- `exact_branch_hessian` rejects the nonzero off-block cases, as it should
- the invalid exact-branch proxy would remain negative and miss the sign change

See:

- `outputs/frontier_sweep.csv`
- `outputs/hessian_breakpoint.svg`
- `outputs/summary.json`
- `outputs/validation.json`
- `outputs/audit.json`
