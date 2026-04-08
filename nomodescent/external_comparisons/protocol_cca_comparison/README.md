# Protocol Mismatch Versus CCA

This pack compares the QD core against canonical correlation analysis on a tiny
real two-panel measurement bundle.

External method:

- canonical correlation analysis between the sepal and petal panels of a tiny
  Iris slice

QD target:

- exact observer-relation classification
- exact common-completion classification
- exact minimal common refinement inside a finite candidate family

Headline:

- CCA reports strong cross-panel dependence
- QD says the panels are still exactly non-nested observers
- QD also says the current evidence leaves the common completion honestly
  underdetermined

Run:

```bash
python -m external_comparisons.protocol_cca_comparison.run_main
python -m external_comparisons.protocol_cca_comparison.validate
```
