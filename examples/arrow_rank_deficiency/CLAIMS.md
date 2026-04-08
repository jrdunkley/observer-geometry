# Claims

## Exact Setup

This is a synthetic latent preference geometry with latent precision:

```text
H = [[4.0, 1.8, 1.4],
     [1.8, 3.0, 0.1],
     [1.4, 0.1, 2.2]]
```

Latent feature directions:

- first-place dominance
- broad approval support
- transfer / runoff structure

Observer maps:

```text
plurality     = [[1, 0, 0]]
approval      = [[1, 0, 0],
                 [0, 1, 0]]
ranked_choice = I_3
```

For each observer `C`, the release computes:

- `Phi_C(H) = (C H^{-1} C^T)^(-1)`
- the hidden load beneath the natural observer ceiling block
- the determinant-gap clock

## Executable Statement

In this synthetic latent model:

- plurality erases the most latent preference geometry
- approval erases less
- ranked choice erases none, up to floating-point residual
- ranked-choice -> approval -> plurality descent matches direct plurality via
  the explicit coarsening maps
  `[[1, 0, 0], [0, 1, 0]]` and `[[1, 0]]`

## Limitations

- synthetic example only
- no empirical election claim
- no normative ranking of voting systems

## Reproduction

```bash
pip install -e .
python -m examples.arrow_rank_deficiency.run_all
```

Inspect:

- `outputs/observer_maps.csv`
- `outputs/observer_summary.csv`
- `outputs/summary.json`
- `outputs/audit.json`

## Acceptance Thresholds

- composition residual: `<= 1e-10`
- max clock-vs-direct-determinant-gap residual: `<= 1e-10`
- ranked-choice clock absolute value: `<= 1e-10`
- ordering flag `plurality >= approval >= ranked_choice`: `true`

The release object passes if all four checks pass.
