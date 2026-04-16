# Claims

## Exact Setup

Use the family

```text
A_+ = [[3, e], [ e, 1]]
A_- = [[3,-e], [-e, 1]]
```

with equal weights, `mu = 0`, and the declared one-dimensional observer

```text
B = [[1],
     [0]].
```

## Executable Statement

The declared observer is stationary for all `e` by cancellation of first
variations, but it is not an exact branch for `e != 0`.

The exact graph-chart second variation is

```text
D^2 F(0) = -24 + 8 e^2.
```

Therefore:

- at `e = 1`, the declared observer is a strict local maximum;
- at `e = sqrt(3)`, the certificate is degenerate/vacuous;
- at `e = 2`, the declared observer is a strict local minimum;
- the invalid exact-branch proxy would remain `-24` and miss the sign change.

## Limitations

- local quadratic geometry only
- no global observer optimization
- no full non-Gaussian law selection
- no probability-support event semantics

## Acceptance Thresholds

- closed-form Hessian residual: `<= 1e-10`
- stationarity residual: `<= 1e-10`
- `e = 1` max certificate passes
- `e = sqrt(3)` max and min certificates are vacuous
- `e = 2` min certificate passes
- nonzero cases reject `exact_branch_hessian`
