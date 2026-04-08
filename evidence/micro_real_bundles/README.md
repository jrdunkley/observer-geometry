# Micro-Real Bundles

These bundles are intentionally tiny.

Each one is:

- legally clean to redistribute or reconstruct
- small enough to audit by inspection
- routed through `evidence` curated ingestion before `ProblemSpec` assembly
- paired with a downstream `nomodescent` result and explicit provenance

## Claim Hierarchy

Each micro-real bundle contains at least three layers and they should be read
separately:

1. Exact extracted or manually reconstructed factual material
2. Encoded inference used to turn that material into visible objects or
   observer hypotheses
3. Downstream `nomodescent` conclusion, which may be exact or audited
   approximate depending on the bundle

Current status by bundle:

- Bell counts bundle
  - exact extracted counts
  - encoded Gaussian surrogate
  - raw path gives exact incompatibility
  - projected no-signalling path gives audited approximate common descent
- Iris protocol mismatch
  - exact tiny measurement table
  - encoded covariance objects
  - exact non-nestedness plus exact underdetermination
- leaderboard benchmark slice
  - exact tiny public score table
  - encoded benchmark-observer hypotheses
  - exact initial underdetermination plus exact non-nestedness after explicit
    selection

Current bundles:

- `bell_counts_bundle`
  - real Bell-test counts from a CC BY 4.0 SciPost article
  - shows the difference between raw finite-sample evidence and an explicitly encoded no-signalling projection
- `iris_protocol_mismatch`
  - 12-row Iris slice from the UCI repository
  - shows real panel mismatch and honest underdetermination
- `leaderboard_benchmark_slice`
  - tiny manually reconstructed public benchmark score slice
  - keeps exact scores separate from encoded benchmark-observer hypotheses

Run from the `evidence` workspace:

```bash
python -m micro_real_bundles.bell_counts_bundle.run_main
python -m micro_real_bundles.iris_protocol_mismatch.run_main
python -m micro_real_bundles.leaderboard_benchmark_slice.run_main
```

