# External Comparisons

These packs validate the narrow Quotient Descent operational core against
external methods, tiny real bundles, or standard external calculations.

Current targets:

- `bell_chsh_comparison`
  - external method: correlator-level CHSH / arcsine compatibility logic
  - data: tiny real Bell count table
- `protocol_cca_comparison`
  - external method: canonical correlation analysis
  - data: tiny real Iris measurement slice
- `free_gaussian_rg_comparison`
  - external method: standard Schur-complement Gaussian elimination
  - data: tiny explicit free-Gaussian precision matrix

Run from the `nomodescent` workspace:

```bash
python -m external_comparisons.bell_chsh_comparison.run_main
python -m external_comparisons.protocol_cca_comparison.run_main
python -m external_comparisons.free_gaussian_rg_comparison.run_main
```
