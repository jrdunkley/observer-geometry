# evidence

`evidence` is the disciplined evidence-to-spec layer above `nomodescent`.

Its job is narrow:

- preserve what a source states exactly
- separate that from encoded modelling choices
- preserve unresolved ambiguity instead of flattening it
- assemble an auditable `nomodescent.ProblemSpec` only when the evidence is
  sufficient

It is not a generic paper chatbot, not a PDF parser, and not a broad autonomous
search system.

## Claim Hierarchy

`evidence` is not itself a theorem engine. Its job is to preserve the
boundary between:

- exact extracted facts
- encoded modelling inferences
- open ambiguity

Its outputs therefore sit in two different layers:

- evidence-preservation claims
  - exact about provenance, extraction mode, and authority status
- downstream problem-assembly claims
  - exact only when the evidence really determines a clean `ProblemSpec`
  - otherwise underdetermined by design

Worked examples and micro-real bundles are also distinct:

- worked examples can be synthetic and pedagogical
- micro-real bundles contain tiny real or reconstructed evidence with explicit
  provenance

## Supported v0.3 Inputs

- explicit matrices from papers or notes
- structured numeric tables
- pairwise marginal families
- finite protocol descriptions that map into a finite observer family
- small curated text claims or excerpts

## Install

`nomodescent v0.2.0` must already be installed or available in the environment.

```bash
pip install -e nomodescent
pip install -e evidence
```

## Workflow

1. Collect curated source material into an `EvidenceBundle`.
2. Mark each item as `exact_extraction`, `encoded_inference`, or
   `open_ambiguity`.
3. Assemble to a `ProblemSpec` only if the evidence is sufficient.
4. Pass the result to `nomodescent`.
5. Read the audit trail before reading the conclusion.

Finite-family observer suggestion is intentionally light:

- `infer_observer_candidates(...)` is a deterministic ranking aid
- it is not an automatic scientific selector
- the top-ranked observer hypothesis is never authoritative by itself
- downstream assembly should still be treated as conditional on the chosen observer encoding

## Curated Ingestion

`evidence` now also exposes a narrow curated-ingestion surface for manually
selected excerpts, tables, matrices, protocol notes, and numeric claims.

The ingestion layer is designed to preserve:

- source provenance
- human selection provenance
- exact vs inferred vs ambiguous status
- required human decisions before assembly

See [docs/CURATED_INGESTION.md](docs/CURATED_INGESTION.md).

## Micro-Real Bundles

- [micro_real_bundles/README.md](micro_real_bundles/README.md)

## Worked Examples

- [worked_examples/bell_evidence_encoding](worked_examples/bell_evidence_encoding)
- [worked_examples/replication_protocol_encoding](worked_examples/replication_protocol_encoding)
- [worked_examples/benchmark_blindness_encoding](worked_examples/benchmark_blindness_encoding)

## Docs

- [docs/DOMAIN.md](docs/DOMAIN.md)
- [docs/EPISTEMIC_STATUS.md](docs/EPISTEMIC_STATUS.md)
- [docs/PROBLEM_ASSEMBLY.md](docs/PROBLEM_ASSEMBLY.md)
- [docs/AUDIT_LAYER.md](docs/AUDIT_LAYER.md)
- [docs/CURATED_INGESTION.md](docs/CURATED_INGESTION.md)
- [docs/MESSY_MATERIAL_ENCODING_GUIDE.md](docs/MESSY_MATERIAL_ENCODING_GUIDE.md)


