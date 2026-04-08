from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._utils import as_matrix, matrix_to_list, validate_epistemic_pair
from .exceptions import EvidenceInputError


@dataclass(frozen=True)
class SourceRef:
    source: str
    location: str = ""
    quote: str = ""
    note: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "source": self.source,
            "location": self.location,
            "quote": self.quote,
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "SourceRef":
        return cls(
            source=str(data["source"]),
            location=str(data.get("location", "")),
            quote=str(data.get("quote", "")),
            note=str(data.get("note", "")),
        )


@dataclass(frozen=True)
class ExtractionRecord:
    extraction_mode: str
    epistemic_status: str
    confidence: float = 1.0
    authoritative: bool = False
    parser: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        validate_epistemic_pair(self.extraction_mode, self.epistemic_status, self.authoritative)
        if not 0.0 <= self.confidence <= 1.0:
            raise EvidenceInputError("confidence must lie in [0, 1]")

    def to_dict(self) -> dict[str, object]:
        return {
            "extraction_mode": self.extraction_mode,
            "epistemic_status": self.epistemic_status,
            "confidence": self.confidence,
            "authoritative": self.authoritative,
            "parser": self.parser,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ExtractionRecord":
        return cls(
            extraction_mode=str(data["extraction_mode"]),
            epistemic_status=str(data["epistemic_status"]),
            confidence=float(data.get("confidence", 1.0)),
            authoritative=bool(data.get("authoritative", False)),
            parser=str(data.get("parser", "")),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class TextClaim:
    name: str
    statement: str
    source_ref: SourceRef
    extraction: ExtractionRecord
    tags: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "statement": self.statement,
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "tags": list(self.tags),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "TextClaim":
        return cls(
            name=str(data["name"]),
            statement=str(data["statement"]),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            tags=tuple(str(item) for item in data.get("tags", [])),
        )


@dataclass(frozen=True)
class TableObservation:
    name: str
    columns: tuple[str, ...]
    rows: tuple[str, ...]
    values: tuple[tuple[float, ...], ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    units: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.columns or not self.rows:
            raise EvidenceInputError("table observations must have non-empty rows and columns")
        if len(self.values) != len(self.rows):
            raise EvidenceInputError("table value row count must match rows")
        for row in self.values:
            if len(row) != len(self.columns):
                raise EvidenceInputError("table value column count must match columns")

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "columns": list(self.columns),
            "rows": list(self.rows),
            "values": [list(row) for row in self.values],
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "units": self.units,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "TableObservation":
        return cls(
            name=str(data["name"]),
            columns=tuple(str(item) for item in data["columns"]),
            rows=tuple(str(item) for item in data["rows"]),
            values=tuple(tuple(float(entry) for entry in row) for row in data["values"]),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            units=str(data.get("units", "")),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class MatrixObservation:
    name: str
    matrix_role: str
    matrix_kind: str
    matrix: np.ndarray | list[list[float]]
    observer_name: str | None
    source_ref: SourceRef
    extraction: ExtractionRecord
    notes: str = ""

    def __post_init__(self) -> None:
        if self.matrix_role not in {"visible_object", "ceiling", "latent_matrix"}:
            raise EvidenceInputError("matrix_role must be 'visible_object', 'ceiling', or 'latent_matrix'")
        if self.matrix_kind not in {"covariance", "precision", "generic"}:
            raise EvidenceInputError("matrix_kind must be 'covariance', 'precision', or 'generic'")
        matrix = as_matrix("matrix observation", self.matrix)
        if matrix.shape[0] != matrix.shape[1]:
            raise EvidenceInputError("matrix observations must be square")
        if self.matrix_role != "visible_object" and self.observer_name is not None:
            raise EvidenceInputError("only visible_object observations may reference an observer_name")
        if self.matrix_role == "visible_object" and not self.observer_name:
            raise EvidenceInputError("visible_object observations must reference an observer_name")
        object.__setattr__(self, "matrix", matrix)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "matrix_role": self.matrix_role,
            "matrix_kind": self.matrix_kind,
            "matrix": matrix_to_list(self.matrix),
            "observer_name": self.observer_name,
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "MatrixObservation":
        return cls(
            name=str(data["name"]),
            matrix_role=str(data["matrix_role"]),
            matrix_kind=str(data["matrix_kind"]),
            matrix=data["matrix"],
            observer_name=None if data.get("observer_name") is None else str(data["observer_name"]),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class ProtocolObservation:
    name: str
    facts: tuple[str, ...]
    candidate_family: tuple[str, ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "facts": list(self.facts),
            "candidate_family": list(self.candidate_family),
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ProtocolObservation":
        return cls(
            name=str(data["name"]),
            facts=tuple(str(item) for item in data.get("facts", [])),
            candidate_family=tuple(str(item) for item in data.get("candidate_family", [])),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class ObserverHypothesis:
    name: str
    matrix: np.ndarray | list[list[float]]
    protocol_name: str
    features: tuple[str, ...]
    source_ref: SourceRef
    extraction: ExtractionRecord
    assumption_statements: tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""

    def __post_init__(self) -> None:
        matrix = as_matrix("observer hypothesis matrix", self.matrix)
        if matrix.shape[0] == 0 or matrix.shape[1] == 0:
            raise EvidenceInputError("observer hypothesis matrices must have positive dimensions")
        if np.linalg.matrix_rank(matrix) != matrix.shape[0]:
            raise EvidenceInputError(f"observer hypothesis '{self.name}' must be surjective")
        object.__setattr__(self, "matrix", matrix)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "matrix": matrix_to_list(self.matrix),
            "protocol_name": self.protocol_name,
            "features": list(self.features),
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "assumption_statements": list(self.assumption_statements),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ObserverHypothesis":
        return cls(
            name=str(data["name"]),
            matrix=data["matrix"],
            protocol_name=str(data["protocol_name"]),
            features=tuple(str(item) for item in data.get("features", [])),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            assumption_statements=tuple(str(item) for item in data.get("assumption_statements", [])),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class ConstraintCandidate:
    name: str
    kind: str
    payload: dict[str, object]
    source_ref: SourceRef
    extraction: ExtractionRecord
    exact: bool
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "kind": self.kind,
            "payload": self.payload,
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
            "exact": self.exact,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ConstraintCandidate":
        return cls(
            name=str(data["name"]),
            kind=str(data["kind"]),
            payload=dict(data.get("payload", {})),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
            exact=bool(data["exact"]),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class ExtractionNote:
    name: str
    statement: str
    load_bearing: bool
    source_ref: SourceRef
    extraction: ExtractionRecord

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "statement": self.statement,
            "load_bearing": self.load_bearing,
            "source_ref": self.source_ref.to_dict(),
            "extraction": self.extraction.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ExtractionNote":
        return cls(
            name=str(data["name"]),
            statement=str(data["statement"]),
            load_bearing=bool(data["load_bearing"]),
            source_ref=SourceRef.from_dict(dict(data["source_ref"])),
            extraction=ExtractionRecord.from_dict(dict(data["extraction"])),
        )


@dataclass(frozen=True)
class EvidenceBundle:
    name: str
    latent_dim: int | None
    text_claims: tuple[TextClaim, ...] = field(default_factory=tuple)
    table_observations: tuple[TableObservation, ...] = field(default_factory=tuple)
    matrix_observations: tuple[MatrixObservation, ...] = field(default_factory=tuple)
    protocol_observations: tuple[ProtocolObservation, ...] = field(default_factory=tuple)
    observer_hypotheses: tuple[ObserverHypothesis, ...] = field(default_factory=tuple)
    constraint_candidates: tuple[ConstraintCandidate, ...] = field(default_factory=tuple)
    notes: tuple[ExtractionNote, ...] = field(default_factory=tuple)
    description: str = ""
    provenance: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.latent_dim is not None and self.latent_dim <= 0:
            raise EvidenceInputError("latent_dim must be positive when specified")

        protocol_names = [item.name for item in self.protocol_observations]
        if len(protocol_names) != len(set(protocol_names)):
            raise EvidenceInputError("protocol observation names must be unique")

        hypothesis_names = [item.name for item in self.observer_hypotheses]
        if len(hypothesis_names) != len(set(hypothesis_names)):
            raise EvidenceInputError("observer hypothesis names must be unique")

        protocols = {item.name for item in self.protocol_observations}
        for hypothesis in self.observer_hypotheses:
            if protocols and hypothesis.protocol_name not in protocols:
                raise EvidenceInputError(
                    f"observer hypothesis '{hypothesis.name}' references unknown protocol '{hypothesis.protocol_name}'"
                )

        visible_keys: dict[tuple[str, str, str], np.ndarray] = {}
        for matrix in self.matrix_observations:
            if matrix.matrix_role != "visible_object" or not matrix.extraction.authoritative:
                continue
            key = (matrix.observer_name or "", matrix.matrix_role, matrix.matrix_kind)
            if key in visible_keys and not np.allclose(visible_keys[key], matrix.matrix, atol=1e-10, rtol=1e-10):
                raise EvidenceInputError(
                    f"conflicting authoritative visible evidence for observer '{matrix.observer_name}' and kind '{matrix.matrix_kind}'"
                )
            visible_keys[key] = matrix.matrix

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "latent_dim": self.latent_dim,
            "text_claims": [item.to_dict() for item in self.text_claims],
            "table_observations": [item.to_dict() for item in self.table_observations],
            "matrix_observations": [item.to_dict() for item in self.matrix_observations],
            "protocol_observations": [item.to_dict() for item in self.protocol_observations],
            "observer_hypotheses": [item.to_dict() for item in self.observer_hypotheses],
            "constraint_candidates": [item.to_dict() for item in self.constraint_candidates],
            "notes": [item.to_dict() for item in self.notes],
            "description": self.description,
            "provenance": self.provenance,
            "tags": list(self.tags),
        }

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "EvidenceBundle":
        return cls(
            name=str(data["name"]),
            latent_dim=None if data.get("latent_dim") is None else int(data["latent_dim"]),
            text_claims=tuple(TextClaim.from_dict(item) for item in data.get("text_claims", [])),
            table_observations=tuple(TableObservation.from_dict(item) for item in data.get("table_observations", [])),
            matrix_observations=tuple(MatrixObservation.from_dict(item) for item in data.get("matrix_observations", [])),
            protocol_observations=tuple(ProtocolObservation.from_dict(item) for item in data.get("protocol_observations", [])),
            observer_hypotheses=tuple(ObserverHypothesis.from_dict(item) for item in data.get("observer_hypotheses", [])),
            constraint_candidates=tuple(ConstraintCandidate.from_dict(item) for item in data.get("constraint_candidates", [])),
            notes=tuple(ExtractionNote.from_dict(item) for item in data.get("notes", [])),
            description=str(data.get("description", "")),
            provenance=str(data.get("provenance", "")),
            tags=tuple(str(item) for item in data.get("tags", [])),
        )

    @classmethod
    def from_json(cls, path: str | Path) -> "EvidenceBundle":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
