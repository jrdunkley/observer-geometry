from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .audit import EvidenceAudit
from .results import AssemblyResult
from .specs import (
    ConstraintCandidate,
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    MatrixObservation,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    TableObservation,
    TextClaim,
)


@dataclass(frozen=True)
class QuotedTextExcerpt:
    name: str
    quote: str
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    notes: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    def to_text_claim(self) -> TextClaim:
        return TextClaim(
            name=self.name,
            statement=self.quote,
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
            tags=self.tags,
        )


@dataclass(frozen=True)
class ManualTableExcerpt:
    name: str
    columns: tuple[str, ...]
    rows: tuple[str, ...]
    values: tuple[tuple[float | str, ...], ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    units: str = ""
    notes: str = ""

    def to_table_observation(self) -> TableObservation:
        numeric_rows = tuple(tuple(float(entry) for entry in row) for row in self.values)
        return TableObservation(
            name=self.name,
            columns=self.columns,
            rows=self.rows,
            values=numeric_rows,
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
            units=self.units,
            notes=self.notes,
        )


@dataclass(frozen=True)
class MatrixExcerpt:
    name: str
    matrix_role: str
    matrix_kind: str
    matrix: tuple[tuple[float, ...], ...]
    observer_name: str | None
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    notes: str = ""

    def to_matrix_observation(self) -> MatrixObservation:
        return MatrixObservation(
            name=self.name,
            matrix_role=self.matrix_role,
            matrix_kind=self.matrix_kind,
            matrix=self.matrix,
            observer_name=self.observer_name,
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
            notes=self.notes,
        )


@dataclass(frozen=True)
class ProtocolExcerpt:
    name: str
    facts: tuple[str, ...]
    candidate_family: tuple[str, ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    notes: str = ""

    def to_protocol_observation(self) -> ProtocolObservation:
        return ProtocolObservation(
            name=self.name,
            facts=self.facts,
            candidate_family=self.candidate_family,
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
            notes=self.notes,
        )


@dataclass(frozen=True)
class BenchmarkExcerpt:
    name: str
    tasks: tuple[str, ...]
    systems: tuple[str, ...]
    scores: tuple[tuple[float, ...], ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    units: str = ""
    notes: str = ""

    def to_table_observation(self) -> TableObservation:
        return TableObservation(
            name=self.name,
            columns=self.tasks,
            rows=self.systems,
            values=self.scores,
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
            units=self.units,
            notes=self.notes,
        )


@dataclass(frozen=True)
class NumericClaimExcerpt:
    name: str
    value: float
    units: str
    source_ref: SourceRef
    extraction: ExtractionRecord
    selection_provenance: str = ""
    notes: str = ""

    def to_text_claim(self) -> TextClaim:
        suffix = f" {self.units}".rstrip()
        return TextClaim(
            name=self.name,
            statement=f"{self.name} = {self.value}{suffix}",
            source_ref=_with_selection_note(self.source_ref, self.selection_provenance),
            extraction=self.extraction,
        )


@dataclass(frozen=True)
class CuratedIngestionResult:
    bundle: EvidenceBundle
    audit: EvidenceAudit
    required_human_decisions: tuple[str, ...]
    unresolved_ambiguities: tuple[str, ...]


def build_evidence_bundle_from_excerpts(
    *,
    name: str,
    latent_dim: int | None,
    excerpts: Iterable[
        QuotedTextExcerpt | ManualTableExcerpt | MatrixExcerpt | ProtocolExcerpt | BenchmarkExcerpt | NumericClaimExcerpt
    ],
    observer_hypotheses: Iterable[ObserverHypothesis] = (),
    matrix_observations: Iterable[MatrixObservation] = (),
    constraint_candidates: Iterable[ConstraintCandidate] = (),
    notes: Iterable[ExtractionNote] = (),
    description: str = "",
    provenance: str = "",
    tags: tuple[str, ...] = (),
) -> CuratedIngestionResult:
    text_claims: list[TextClaim] = []
    table_observations: list[TableObservation] = []
    direct_matrices: list[MatrixObservation] = []
    protocol_observations: list[ProtocolObservation] = []

    exact_items: list[str] = []
    inferred_items: list[str] = []
    ambiguous_items: list[str] = []
    authoritative_items: list[str] = []
    unresolved: list[str] = []
    required_human_decisions: list[str] = []

    for excerpt in excerpts:
        _record_excerpt_status(
            excerpt.name,
            excerpt.extraction,
            exact_items,
            inferred_items,
            ambiguous_items,
            authoritative_items,
            unresolved,
            required_human_decisions,
        )
        if isinstance(excerpt, QuotedTextExcerpt):
            text_claims.append(excerpt.to_text_claim())
        elif isinstance(excerpt, NumericClaimExcerpt):
            text_claims.append(excerpt.to_text_claim())
        elif isinstance(excerpt, ManualTableExcerpt):
            table_observations.append(excerpt.to_table_observation())
        elif isinstance(excerpt, BenchmarkExcerpt):
            table_observations.append(excerpt.to_table_observation())
        elif isinstance(excerpt, MatrixExcerpt):
            direct_matrices.append(excerpt.to_matrix_observation())
        elif isinstance(excerpt, ProtocolExcerpt):
            protocol_observations.append(excerpt.to_protocol_observation())
        else:
            raise TypeError(f"unsupported excerpt type: {type(excerpt)!r}")

    bundle = EvidenceBundle(
        name=name,
        latent_dim=latent_dim,
        text_claims=tuple(text_claims),
        table_observations=tuple(table_observations),
        matrix_observations=tuple([*direct_matrices, *matrix_observations]),
        protocol_observations=tuple(protocol_observations),
        observer_hypotheses=tuple(observer_hypotheses),
        constraint_candidates=tuple(constraint_candidates),
        notes=tuple(notes),
        description=description,
        provenance=provenance,
        tags=tags,
    )
    audit = EvidenceAudit(
        exact_items=tuple(sorted(exact_items)),
        inferred_items=tuple(sorted(inferred_items)),
        ambiguous_items=tuple(sorted(ambiguous_items)),
        authoritative_items=tuple(sorted(authoritative_items)),
        load_bearing_items=tuple(sorted([*exact_items, *inferred_items, *ambiguous_items])),
        load_bearing_exact_items=tuple(sorted(exact_items)),
        load_bearing_inferred_items=tuple(sorted(inferred_items)),
        load_bearing_ambiguous_items=tuple(sorted(ambiguous_items)),
        load_bearing_assumptions=tuple(note.statement for note in notes if note.load_bearing),
        unresolved_items=tuple(sorted(unresolved)),
        residuals={"ambiguous_excerpt_count": float(len(unresolved))},
        theorem_layer="curated ingestion preserves epistemic status but does not itself make a descent conclusion",
        falsification_route=(
            "replace an encoded inference excerpt with an exact excerpt from the source if the encoding is too strong",
            "downgrade any item to open_ambiguity if the source does not determine it cleanly",
        ),
        ambiguity_collapse_route=(
            "add one more curated excerpt that distinguishes between the remaining candidate families",
        ),
        notes=("curated ingestion keeps exact, inferred, and ambiguous material separate before ProblemSpec assembly",),
    )
    return CuratedIngestionResult(
        bundle=bundle,
        audit=audit,
        required_human_decisions=tuple(required_human_decisions),
        unresolved_ambiguities=tuple(sorted(unresolved)),
    )


def _record_excerpt_status(
    name: str,
    extraction: ExtractionRecord,
    exact_items: list[str],
    inferred_items: list[str],
    ambiguous_items: list[str],
    authoritative_items: list[str],
    unresolved: list[str],
    required_human_decisions: list[str],
) -> None:
    if extraction.authoritative:
        authoritative_items.append(name)
    if extraction.extraction_mode == "exact_extraction":
        exact_items.append(name)
    elif extraction.extraction_mode == "encoded_inference":
        inferred_items.append(name)
    else:
        ambiguous_items.append(name)
        unresolved.append(name)
        required_human_decisions.append(f"resolve open ambiguity in excerpt '{name}' before treating it as authoritative")


def _with_selection_note(source_ref: SourceRef, selection_provenance: str) -> SourceRef:
    if not selection_provenance:
        return source_ref
    note_parts = [part for part in (source_ref.note, f"selection provenance: {selection_provenance}") if part]
    return SourceRef(
        source=source_ref.source,
        location=source_ref.location,
        quote=source_ref.quote,
        note=" | ".join(note_parts),
    )
