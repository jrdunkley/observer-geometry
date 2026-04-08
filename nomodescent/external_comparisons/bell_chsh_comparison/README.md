# Bell CHSH Comparison

This pack compares the Quotient Descent core against a standard correlator-level
Bell method on a tiny real Bell count table.

External method:

- CHSH / arcsine compatibility logic on the normalized correlator summary

QD target:

- full-law Gaussian-surrogate common descent on the pair covariances derived
  from the same counts

Headline:

- the CHSH-side coarse method is non-obstructive on the real count table
- raw full-law QD still detects exact incompatibility from repeated-marginal /
  linear inconsistency
- after an explicit no-signalling projection, QD returns audited approximate
  compatibility

Run:

```bash
python -m external_comparisons.bell_chsh_comparison.run_main
python -m external_comparisons.bell_chsh_comparison.validate
```
