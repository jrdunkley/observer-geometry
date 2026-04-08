# Bell Counts Bundle

This bundle is a tiny real Bell-family evidence object.

What is included:

- exact extracted coincidence counts for four Bell contexts
- exact context selectors as observer maps
- exact arithmetic pair probabilities derived from the counts
- encoded Gaussian-surrogate covariance matrices for downstream `nomodescent`
- an explicit no-signalling projection used only as a documented modelling layer

Primary run:

```bash
python -m micro_real_bundles.bell_counts_bundle.run_main
python -m micro_real_bundles.bell_counts_bundle.validate
```

What is shown:

- raw exact counts lead to a linear repeated-marginal inconsistency certificate under Gaussian common descent
- an explicitly encoded no-signalling projection removes the finite-sample signalling mismatch and yields an audited approximate common descent

What is not claimed:

- this bundle does not claim that the underlying Bell experiment is Gaussian
- the projected problem is an encoded Gaussian surrogate, not a theorem about the original non-Gaussian experiment
