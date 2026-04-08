from __future__ import annotations

from evidence import (
    ExtractionRecord,
    ManualTableExcerpt,
    MatrixExcerpt,
    ProtocolExcerpt,
    QuotedTextExcerpt,
    SourceRef,
    build_evidence_bundle_from_excerpts,
)


def test_curated_ingestion_preserves_selection_provenance() -> None:
    excerpt = QuotedTextExcerpt(
        name="quoted_setup",
        quote="The observer measures x and y.",
        source_ref=SourceRef(source="paper", location="Methods", quote="The observer measures x and y."),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        selection_provenance="manual excerpt chosen from methods paragraph",
    )
    result = build_evidence_bundle_from_excerpts(name="bundle", latent_dim=2, excerpts=(excerpt,))
    assert result.bundle.text_claims[0].source_ref.note == "selection provenance: manual excerpt chosen from methods paragraph"
    assert result.audit.exact_items == ("quoted_setup",)
    assert result.audit.load_bearing_exact_items == ("quoted_setup",)


def test_open_ambiguity_excerpt_requires_human_decision() -> None:
    excerpt = ProtocolExcerpt(
        name="protocol_a",
        facts=("measure_x",),
        candidate_family=("obs_1", "obs_2"),
        source_ref=SourceRef(source="paper", location="Methods"),
        extraction=ExtractionRecord(extraction_mode="open_ambiguity", epistemic_status="ambiguous"),
        selection_provenance="manual protocol excerpt",
    )
    result = build_evidence_bundle_from_excerpts(name="bundle", latent_dim=2, excerpts=(excerpt,))
    assert result.unresolved_ambiguities == ("protocol_a",)
    assert result.required_human_decisions
    assert result.audit.ambiguous_items == ("protocol_a",)
    assert result.audit.load_bearing_ambiguous_items == ("protocol_a",)


def test_curated_ingestion_can_emit_mixed_bundle() -> None:
    table = ManualTableExcerpt(
        name="score_table",
        columns=("x", "y"),
        rows=("r1", "r2"),
        values=((1.0, 2.0), (3.0, 4.0)),
        source_ref=SourceRef(source="paper", location="Table 1"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        selection_provenance="manual table transcription",
    )
    matrix = MatrixExcerpt(
        name="cov_obs",
        matrix_role="visible_object",
        matrix_kind="covariance",
        matrix=((1.0, 0.2), (0.2, 1.0)),
        observer_name="obs",
        source_ref=SourceRef(source="paper", location="Appendix"),
        extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
        selection_provenance="derived covariance from extracted table",
    )
    result = build_evidence_bundle_from_excerpts(name="bundle", latent_dim=2, excerpts=(table, matrix))
    assert result.bundle.table_observations[0].name == "score_table"
    assert result.bundle.matrix_observations[0].name == "cov_obs"
    assert result.audit.exact_items == ("score_table",)
    assert result.audit.inferred_items == ("cov_obs",)
    assert result.audit.load_bearing_items == ("cov_obs", "score_table")

