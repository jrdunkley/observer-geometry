# LLM Audit Brief

## Mathematical Object

Compare correlator-level Gaussian Bell compatibility with full-law Gaussian
common-gluing compatibility on an explicit Bell-square family.

## Commands

```bash
pip install -e .
python -m examples.bell_common_gluing.run_all
```

## Expected Headline Outputs

- `outputs/phase_diagram.svg`
- `outputs/comparison_table.csv`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- compatible and variance-only correlator matrices match to `<= 1e-12`
- compatible re-observation residual through `nomogeo` is `<= 1e-10`
- variance-only linear completion residual is `>= 1e-3`
- correlator-only best PSD gap is negative by at least `1e-3`
- confirmed variance-only neighborhood flag is `true`

## Narrow Audit Question

Does the output support the claim that there is a strict gap between
correlator-level compatibility and full-law Gaussian compatibility in the
explicit family used here?
