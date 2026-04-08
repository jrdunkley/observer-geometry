from __future__ import annotations

from .specs import (
    ExtractionRecord,
    MatrixObservation,
    ProtocolObservation,
    SourceRef,
    TableObservation,
)


def encode_table_observation(
    *,
    name: str,
    columns: tuple[str, ...],
    rows: tuple[str, ...],
    values: tuple[tuple[float, ...], ...],
    source: str,
    location: str = "",
    quote: str = "",
    extraction_mode: str = "exact_extraction",
    epistemic_status: str = "exact",
    authoritative: bool = True,
    notes: str = "",
) -> TableObservation:
    return TableObservation(
        name=name,
        columns=columns,
        rows=rows,
        values=values,
        source_ref=SourceRef(source=source, location=location, quote=quote),
        extraction=ExtractionRecord(
            extraction_mode=extraction_mode,
            epistemic_status=epistemic_status,
            authoritative=authoritative,
            notes=notes,
        ),
        notes=notes,
    )


def encode_matrix_observation(
    *,
    name: str,
    matrix_role: str,
    matrix_kind: str,
    matrix: list[list[float]] | tuple[tuple[float, ...], ...],
    source: str,
    observer_name: str | None = None,
    location: str = "",
    quote: str = "",
    extraction_mode: str = "exact_extraction",
    epistemic_status: str = "exact",
    authoritative: bool = True,
    notes: str = "",
) -> MatrixObservation:
    return MatrixObservation(
        name=name,
        matrix_role=matrix_role,
        matrix_kind=matrix_kind,
        matrix=matrix,
        observer_name=observer_name,
        source_ref=SourceRef(source=source, location=location, quote=quote),
        extraction=ExtractionRecord(
            extraction_mode=extraction_mode,
            epistemic_status=epistemic_status,
            authoritative=authoritative,
            notes=notes,
        ),
        notes=notes,
    )


def encode_protocol_observation(
    *,
    name: str,
    facts: tuple[str, ...],
    candidate_family: tuple[str, ...],
    source: str,
    location: str = "",
    quote: str = "",
    extraction_mode: str = "exact_extraction",
    epistemic_status: str = "exact",
    authoritative: bool = True,
    notes: str = "",
) -> ProtocolObservation:
    return ProtocolObservation(
        name=name,
        facts=facts,
        candidate_family=candidate_family,
        source_ref=SourceRef(source=source, location=location, quote=quote),
        extraction=ExtractionRecord(
            extraction_mode=extraction_mode,
            epistemic_status=epistemic_status,
            authoritative=authoritative,
            notes=notes,
        ),
        notes=notes,
    )
