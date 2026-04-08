# Audit Layer

Every `AssemblyResult` and `SuggestionResult` carries an `EvidenceAudit`.

The audit records:

- exact items
- inferred items
- ambiguous items
- authoritative items
- load-bearing items
- load-bearing exact items
- load-bearing inferred items
- load-bearing ambiguous items
- load-bearing assumptions
- unresolved items
- residuals
- theorem layer
- falsification route
- ambiguity collapse route

This is the main honesty mechanism in the package. The conclusion should never
be read without the audit.

Important distinction:

- `exact_items` / `inferred_items` / `ambiguous_items` describe what is present
  in the evidence object
- `load_bearing_*` fields describe what actually matters for the assembled path
  or suggestion result being reported

This matters especially when a bundle contains several candidate observer
hypotheses but only one is selected for downstream assembly.
