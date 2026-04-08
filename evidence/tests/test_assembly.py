from __future__ import annotations

import pytest

from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    assemble_problem_spec,
    encode_matrix_observation,
)
from evidence.exceptions import EvidenceAssemblyError, EvidenceInputError


def _ambiguous_bundle() -> EvidenceBundle:
    protocol = ProtocolObservation(
        name="p",
        facts=("a", "b"),
        candidate_family=("h1", "h2"),
        source_ref=SourceRef(source="paper", location="methods"),
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


def test_ambiguous_bundle_returns_underdetermined() -> None:
    result = assemble_problem_spec(_ambiguous_bundle())
    assert result.classification == "underdetermined_evidence"
    assert result.problem_spec is None
    assert result.required_human_decisions


def test_assembly_succeeds_after_explicit_observer_selection() -> None:
    result = assemble_problem_spec(_ambiguous_bundle(), observer_selection={"p": "h1"})
    assert result.classification == "assembled_problem_spec"
    assert result.problem_spec is not None
    assert result.problem_spec.evidence[0].observer == "h1"
    assert "choose h1" in result.audit.load_bearing_assumptions
    assert "h1" in result.audit.load_bearing_items
    assert "h2" not in result.audit.load_bearing_items
    assert "h1" in result.audit.load_bearing_inferred_items
    assert result.exact is False


def test_insufficient_evidence_returns_structured_refusal() -> None:
    bundle = _ambiguous_bundle()
    bundle = EvidenceBundle(
        name=bundle.name,
        latent_dim=bundle.latent_dim,
        protocol_observations=bundle.protocol_observations,
        observer_hypotheses=bundle.observer_hypotheses,
    )
    result = assemble_problem_spec(bundle, observer_selection={"p": "h1"})
    assert result.classification == "insufficient_evidence"
    assert result.problem_spec is None


def test_conflicting_authoritative_evidence_is_rejected() -> None:
    matrix_a = encode_matrix_observation(
        name="cov_a",
        matrix_role="visible_object",
        matrix_kind="covariance",
        matrix=((1.0,),),
        observer_name="obs",
        source="paper",
    )
    matrix_b = encode_matrix_observation(
        name="cov_b",
        matrix_role="visible_object",
        matrix_kind="covariance",
        matrix=((2.0,),),
        observer_name="obs",
        source="paper",
    )
    protocol = ProtocolObservation(
        name="obs",
        facts=("a",),
        candidate_family=("obs",),
        source_ref=SourceRef(source="paper"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
    )
    hypothesis = ObserverHypothesis(
        name="obs",
        matrix=[[1.0]],
        protocol_name="obs",
        features=("a",),
        source_ref=SourceRef(source="paper"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
    )
    with pytest.raises(EvidenceInputError, match="conflicting authoritative visible evidence"):
        EvidenceBundle(
            name="conflict",
            latent_dim=1,
            matrix_observations=(matrix_a, matrix_b),
            protocol_observations=(protocol,),
            observer_hypotheses=(hypothesis,),
        )


def test_invalid_observer_selection_is_rejected() -> None:
    with pytest.raises(EvidenceAssemblyError, match="not a valid candidate"):
        assemble_problem_spec(_ambiguous_bundle(), observer_selection={"p": "missing"})

