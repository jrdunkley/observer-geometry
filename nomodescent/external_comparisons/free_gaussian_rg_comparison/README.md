# Free Gaussian RG Comparison

This pack compares the QD core against standard Schur-complement elimination in
the free Gaussian regime.

External method:

- block Gaussian elimination / Schur-complement elimination of hidden modes

QD target:

- exact quotient descent under retained-mode observation
- exact tower agreement under staged elimination
- exact common refinement / non-nestedness distinctions for alternative coarse
  observers

Headline:

- QD reproduces the standard Schur-complement elimination exactly
- QD also makes the observer relation structure explicit, including when block
  averaging is non-nested with retained-mode observation

Run:

```bash
python -m external_comparisons.free_gaussian_rg_comparison.run_main
python -m external_comparisons.free_gaussian_rg_comparison.validate
```
