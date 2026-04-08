# Arrow As Observer Rank Deficiency

Reproduce from a clean installed environment:

```bash
pip install -e .
python -m examples.arrow_rank_deficiency.run_all
```

Then inspect:

- [summary.json](outputs/summary.json)
- [audit.json](outputs/audit.json)

## Scope Limitation

This is a synthetic observer-geometry example. It is not an empirical claim
about real elections and not a new theorem in political science.

## What Is Shown

Plurality, approval, and ranked choice are represented as explicit observer
maps on one latent preference geometry. Their hidden loads then quantify which
latent preference directions remain erased by each observer.

The narrow message is:

which latent preference directions are erased by each observer?

The explicit descent checked in the validator is:

`ranked_choice -> approval -> plurality`

## What Is Not Claimed

- no claim about which voting system is normatively best
- no claim about real electorates
- no claim beyond the explicit synthetic latent model used here

## Why This Is Difficult To Dismiss

- the observer maps are literal matrices
- the hidden-load clocks are cross-checked against direct determinant gaps
- richer-to-coarser observer descent is checked exactly
- the example is explicit enough for independent inspection without political
  interpretation

## Observer Maps

```text
plurality    = [[1, 0, 0]]
approval     = [[1, 0, 0],
                [0, 1, 0]]
ranked_choice= [[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]
```

## Public Outputs

- headline figure:
  [clock_comparison.svg](outputs/clock_comparison.svg)
- residual / structure figure:
  [hidden_spectra.svg](outputs/hidden_spectra.svg)
- summary table:
  [observer_summary.csv](outputs/observer_summary.csv)

See also:

- [CLAIMS.md](CLAIMS.md)
- [LLM_AUDIT_BRIEF.md](LLM_AUDIT_BRIEF.md)
- [RESULT.md](RESULT.md)

