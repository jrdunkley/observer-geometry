# Curated Ingestion

The curated-ingestion layer is intentionally narrow.

Supported excerpt classes:

- `QuotedTextExcerpt`
- `ManualTableExcerpt`
- `MatrixExcerpt`
- `ProtocolExcerpt`
- `BenchmarkExcerpt`
- `NumericClaimExcerpt`

Each excerpt carries:

- source provenance
- human selection provenance
- extraction mode
- epistemic status
- downstream authority flag

The ingestion layer preserves three modes and does not let them blur:

- `exact_extraction`
  - use only when the source explicitly states the item
- `encoded_inference`
  - use when the item is constructed from the source through a documented modelling or encoding choice
- `open_ambiguity`
  - use when the excerpt does not determine a clean downstream object

`build_evidence_bundle_from_excerpts(...)` produces:

- an `EvidenceBundle`
- an `EvidenceAudit`
- unresolved ambiguity list
- required human decisions

It does not assemble a `ProblemSpec` by itself. That remains a separate step.

This separation is deliberate:

- curated ingestion preserves what the source gives you
- `assemble_problem_spec(...)` decides whether that preserved evidence is sufficient for downstream descent
