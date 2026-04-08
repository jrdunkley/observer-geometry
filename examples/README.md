# Examples

This folder contains two example classes:

- small kernel examples
- public demonstration examples with frozen outputs and audits

These examples do not all have the same epistemic status. See
[docs/claim_hierarchy.md](../docs/claim_hierarchy.md)
before reusing any result.

## Quick Run

From repo root:

```bash
python -m examples.entanglement_hidden_load.run_all
python -m examples.bell_common_gluing.run_all
python -m examples.arrow_rank_deficiency.run_all
```

Then inspect `outputs/summary.json` and `outputs/audit.json` in the relevant example folder.

## Public Demonstrations

### Entanglement As Hidden Load

- type: exact Gaussian identity demonstration
- command: `python -m examples.entanglement_hidden_load.run_all`
- surface:
  [README.md](entanglement_hidden_load/README.md),
  [CLAIMS.md](entanglement_hidden_load/CLAIMS.md),
  [LLM_AUDIT_BRIEF.md](entanglement_hidden_load/LLM_AUDIT_BRIEF.md),
  [audit.json](entanglement_hidden_load/outputs/audit.json)

### Bell / Common-Gluing Failure Beyond Correlators

- type: structural separation demonstration
- command: `python -m examples.bell_common_gluing.run_all`
- surface:
  [README.md](bell_common_gluing/README.md),
  [CLAIMS.md](bell_common_gluing/CLAIMS.md),
  [LLM_AUDIT_BRIEF.md](bell_common_gluing/LLM_AUDIT_BRIEF.md),
  [audit.json](bell_common_gluing/outputs/audit.json)

### Arrow As Observer Rank Deficiency

- type: synthetic worked example
- command: `python -m examples.arrow_rank_deficiency.run_all`
- surface:
  [README.md](arrow_rank_deficiency/README.md),
  [CLAIMS.md](arrow_rank_deficiency/CLAIMS.md),
  [LLM_AUDIT_BRIEF.md](arrow_rank_deficiency/LLM_AUDIT_BRIEF.md),
  [audit.json](arrow_rank_deficiency/outputs/audit.json)

## Small Kernel Examples

- `python examples/minimal_visible_precision.py`
- `python examples/local_calculus.py`
- `python examples/hidden_load_transport.py`
- `python examples/finite_dv_bridge.py`
- `python examples/quotient_coarsening.py`

## Index Files

- [RELEASE_NOTE.md](RELEASE_NOTE.md)
- [release_manifest.json](release_manifest.json)


