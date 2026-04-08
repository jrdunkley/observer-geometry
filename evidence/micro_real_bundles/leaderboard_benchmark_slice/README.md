# Leaderboard Benchmark Slice

This bundle is a tiny manually reconstructed public benchmark slice.

What is shown:

- exact task-by-model scores are kept distinct from the encoded observer family
- benchmark suites remain underdetermined until an observer hypothesis is selected
- after explicit selection, the two benchmark suites are non-nested observers and the common Gaussian completion remains underdetermined

Primary run:

```bash
python -m micro_real_bundles.leaderboard_benchmark_slice.run_main
python -m micro_real_bundles.leaderboard_benchmark_slice.validate
```

What is not claimed:

- this is not a theorem about LLM capability structure
- it is only a disciplined encoding of a tiny public benchmark score slice
