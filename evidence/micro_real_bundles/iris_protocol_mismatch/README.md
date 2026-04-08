# Iris Protocol Mismatch

This is a tiny real protocol-mismatch bundle built from a 12-row Iris slice.

What is shown:

- the sepal panel and petal panel are exact non-nested observers of the same 4-feature measurement space
- the two visible covariances are exactly compatible with many common latent completions, so the honest downstream result is underdetermination
- the minimal obvious common refinement is the full 4-feature panel

Primary run:

```bash
python -m micro_real_bundles.iris_protocol_mismatch.run_main
python -m micro_real_bundles.iris_protocol_mismatch.validate
```

What is not claimed:

- this is not a claim about the replication crisis
- it is a tiny real measurement-panel mismatch example only
