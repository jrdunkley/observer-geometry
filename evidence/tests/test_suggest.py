from __future__ import annotations

from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    encode_matrix_observation,
    infer_observer_candidates,
)


def _ambiguous_bundle() -> EvidenceBundle:
    protocol = ProtocolObservation(
        name="p",
        facts=("a", "b"),
        candidate_family=("h1", "h2"),
        source_ref=SourceRef(source="paper"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
    )
    h1 = ObserverHypothesis(
        name="h1",
        matrix=[[1.0, 0.0]],
        protocol_name="p",
        features=("a",),
        source_ref=SourceRef(source="model"),
        extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
        assumption_statements=("choose h1",),
    )
    h2 = ObserverHypothesis(
        name="h2",
        matrix=[[0.0, 1.0]],
        protocol_name="p",
        features=("b",),
        source_ref=SourceRef(source="model"),
        extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
        assumption_statements=("choose h2",),
    )
    matrix = encode_matrix_observation(
        name="cov_p",
        matrix_role="visible_object",
        matrix_kind="covariance",
        matrix=((1.0,),),
        observer_name="p",
        source="paper",
    )
    note = ExtractionNote(
        name="model_choice",
        statement="observer matrix is inferred from protocol wording",
        load_bearing=True,
        source_ref=SourceRef(source="model"),
        extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred"),
    )
    return EvidenceBundle(
        name="ambiguous",
        latent_dim=2,
        matrix_observations=(matrix,),
        protocol_observations=(protocol,),
        observer_hypotheses=(h1, h2),
        notes=(note,),
    )


def test_encoded_inference_records_assumptions_and_ranks_deterministically() -> None:
    result = infer_observer_candidates(_ambiguous_bundle(), "p")
    assert result.ranked_candidates == ("h1", "h2")
    assert result.audit.theorem_layer.startswith("finite-family evidence suggestion")
    assert result.assumptions
    assert result.audit.inferred_items == ("h1", "h2")
    assert any("heuristic" in note for note in result.audit.notes)
    assert any("top-ranked candidate" in step for step in result.audit.falsification_route)

