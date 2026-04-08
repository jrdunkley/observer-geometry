# LLM Audit Brief

## Mathematical Object

Compute visible precisions and hidden loads for three explicit observer maps on
a synthetic latent preference geometry.

## Commands

```bash
pip install -e .
python -m examples.arrow_rank_deficiency.run_all
```

## Expected Headline Outputs

- `outputs/observer_maps.csv`
- `outputs/observer_summary.csv`
- `outputs/clock_comparison.svg`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- composition residual `<= 1e-10`
- determinant-gap cross-check residual `<= 1e-10`
- ranked-choice clock absolute value `<= 1e-10`
- clock ordering flag is `true`

## Narrow Audit Question

Does the output support the narrow synthetic claim that these three explicit
observer maps erase different latent preference directions, with exact
determinant-gap agreement and exact nested descent?
