# Claims

## Exact Setup

The Bell square contains four centered Gaussian pair laws:

`(A_0, B_0), (A_0, B_1), (A_1, B_0), (A_1, B_1)`.

For parameters `rho in [0, 1)` and `delta >= 0`, this release uses:

```text
Sigma_00 = [[1, rho], [rho, 1]]
Sigma_01 = [[1 + delta, rho sqrt(1 + delta)], [rho sqrt(1 + delta), 1]]
Sigma_10 = [[1, rho], [rho, 1]]
Sigma_11 = [[1, -rho], [-rho, 1]]
```

The normalized correlator matrix is:

```text
K = [[ rho,  rho],
     [ rho, -rho]]
```

and therefore does not depend on `delta`.

The correlator-level criterion is the arcsine margin:

`pi - |asin(k00) + asin(k01) + asin(k10) - asin(k11)|`.

The second obstruction is repeated-marginal variance consistency across context.

## Executable Statement

The saved outputs support the following narrow claim:

There exists an explicit Gaussian family where the normalized correlator summary
is unchanged, the correlator-level Bell obstruction is absent, and common
Gaussian gluing still fails because of a second obstruction invisible to the
correlator summary.

The release object certifies this with:

- a compatible sample
- a variance-only incompatible sample
- a correlator-only incompatible sample
- a phase sweep
- a confirmed variance-only neighborhood on a `5 x 5` grid over:
  - `delta in [0.18, 0.26]`
  - `rho in [0.56, 0.66]`

## Limitations

- Gaussian common-gluing only
- no non-Gaussian Bell claims
- no claim beyond the explicit family and validations stored here

## Reproduction

```bash
pip install -e .
python -m examples.bell_common_gluing.run_all
```

Inspect:

- `outputs/comparison_table.csv`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- compatible vs variance-only correlator matrix match residual: `<= 1e-12`
- compatible pair re-observation residual through `nomogeo`: `<= 1e-10`
- compatible linear completion residual: `<= 1e-10`
- variance-only linear completion residual: `>= 1e-3`
- correlator-only best PSD gap: `<= -1e-3`
- confirmed variance-only neighborhood flag: `true`

The release object passes if all six checks pass.
