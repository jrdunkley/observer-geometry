from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .exceptions import DescentInputError


def _as_matrix(name: str, value: Any) -> np.ndarray:
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2:
        raise DescentInputError(f"{name} must be a 2D array")
    return matrix


def _matrix_to_list(value: np.ndarray) -> list[list[float]]:
    return [[float(entry) for entry in row] for row in np.asarray(value, dtype=float)]


@dataclass(frozen=True)
class ObserverSpec:
    name: str
    matrix: np.ndarray | list[list[float]]
    provenance: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        matrix = _as_matrix("observer matrix", self.matrix)
        if matrix.shape[0] == 0 or matrix.shape[1] == 0:
            raise DescentInputError("observer matrix must have positive dimensions")
        if np.linalg.matrix_rank(matrix) != matrix.shape[0]:
            raise DescentInputError(f"observer '{self.name}' must be surjective (full row rank)")
        object.__setattr__(self, "matrix", matrix)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "matrix": _matrix_to_list(self.matrix),
            "provenance": self.provenance,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ObserverSpec":
        return cls(
            name=str(data["name"]),
            matrix=data["matrix"],
            provenance=str(data.get("provenance", "")),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class VisibleEvidenceSpec:
    name: str
    observer: str
    kind: str
    matrix: np.ndarray | list[list[float]]
    exact: bool = True
    tolerance: float | None = None
    authoritative: bool = True
    provenance: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        matrix = _as_matrix("visible evidence matrix", self.matrix)
        if self.kind not in {"covariance", "precision"}:
            raise DescentInputError("visible evidence kind must be 'covariance' or 'precision'")
        if matrix.shape[0] != matrix.shape[1]:
            raise DescentInputError("visible evidence matrix must be square")
        if self.tolerance is not None and self.tolerance < 0.0:
            raise DescentInputError("visible evidence tolerance must be non-negative")
        object.__setattr__(self, "matrix", matrix)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "observer": self.observer,
            "kind": self.kind,
            "matrix": _matrix_to_list(self.matrix),
            "exact": self.exact,
            "tolerance": self.tolerance,
            "authoritative": self.authoritative,
            "provenance": self.provenance,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "VisibleEvidenceSpec":
        return cls(
            name=str(data["name"]),
            observer=str(data["observer"]),
            kind=str(data["kind"]),
            matrix=data["matrix"],
            exact=bool(data.get("exact", True)),
            tolerance=None if data.get("tolerance") is None else float(data["tolerance"]),
            authoritative=bool(data.get("authoritative", True)),
            provenance=str(data.get("provenance", "")),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class CeilingSpec:
    name: str
    observer: str
    matrix: np.ndarray | list[list[float]]
    exact: bool = True
    provenance: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", _as_matrix("ceiling matrix", self.matrix))

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "observer": self.observer,
            "matrix": _matrix_to_list(self.matrix),
            "exact": self.exact,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "CeilingSpec":
        return cls(
            name=str(data["name"]),
            observer=str(data["observer"]),
            matrix=data["matrix"],
            exact=bool(data.get("exact", True)),
            provenance=str(data.get("provenance", "")),
        )


@dataclass(frozen=True)
class ConstraintSpec:
    name: str
    kind: str
    exact: bool
    payload: dict[str, object]
    provenance: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "kind": self.kind,
            "exact": self.exact,
            "payload": self.payload,
            "provenance": self.provenance,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ConstraintSpec":
        return cls(
            name=str(data["name"]),
            kind=str(data["kind"]),
            exact=bool(data["exact"]),
            payload=dict(data.get("payload", {})),
            provenance=str(data.get("provenance", "")),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class AssumptionEntry:
    label: str
    statement: str
    exact: bool
    provenance: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "statement": self.statement,
            "exact": self.exact,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "AssumptionEntry":
        return cls(
            label=str(data["label"]),
            statement=str(data["statement"]),
            exact=bool(data["exact"]),
            provenance=str(data.get("provenance", "")),
        )


@dataclass(frozen=True)
class AssumptionLedger:
    entries: tuple[AssumptionEntry, ...] = field(default_factory=tuple)

    def exact_entries(self) -> tuple[str, ...]:
        return tuple(entry.statement for entry in self.entries if entry.exact)

    def approximate_entries(self) -> tuple[str, ...]:
        return tuple(entry.statement for entry in self.entries if not entry.exact)

    def to_dict(self) -> dict[str, object]:
        return {"entries": [entry.to_dict() for entry in self.entries]}

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "AssumptionLedger":
        return cls(entries=tuple(AssumptionEntry.from_dict(item) for item in data.get("entries", [])))


@dataclass(frozen=True)
class GoalSpec:
    kind: str
    target: str | None = None
    candidates: tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""

    def __post_init__(self) -> None:
        allowed = {
            "factorisation",
            "relation_classification",
            "common_completion",
            "tower_check",
            "minimal_refinement",
            "staged_descent",
        }
        if self.kind not in allowed:
            raise DescentInputError(f"goal kind must be one of {sorted(allowed)}")

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "target": self.target,
            "candidates": list(self.candidates),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "GoalSpec":
        return cls(
            kind=str(data["kind"]),
            target=None if data.get("target") is None else str(data["target"]),
            candidates=tuple(str(item) for item in data.get("candidates", [])),
            notes=str(data.get("notes", "")),
        )


@dataclass(frozen=True)
class ProblemSpec:
    name: str
    latent_dim: int | None
    observers: tuple[ObserverSpec, ...]
    evidence: tuple[VisibleEvidenceSpec, ...] = field(default_factory=tuple)
    ceilings: tuple[CeilingSpec, ...] = field(default_factory=tuple)
    constraints: tuple[ConstraintSpec, ...] = field(default_factory=tuple)
    assumptions: AssumptionLedger = field(default_factory=AssumptionLedger)
    goals: tuple[GoalSpec, ...] = field(default_factory=tuple)
    description: str = ""
    provenance: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.latent_dim is not None and self.latent_dim <= 0:
            raise DescentInputError("latent_dim must be positive when specified")
        observer_names = [observer.name for observer in self.observers]
        if len(observer_names) != len(set(observer_names)):
            raise DescentInputError("observer names must be unique")
        observer_map = {observer.name: observer for observer in self.observers}

        for evidence in self.evidence:
            if evidence.observer not in observer_map:
                raise DescentInputError(f"evidence '{evidence.name}' references unknown observer '{evidence.observer}'")
            if observer_map[evidence.observer].matrix.shape[0] != evidence.matrix.shape[0]:
                raise DescentInputError(
                    f"evidence '{evidence.name}' dimension does not match observer '{evidence.observer}'"
                )
        for ceiling in self.ceilings:
            if ceiling.observer not in observer_map:
                raise DescentInputError(f"ceiling '{ceiling.name}' references unknown observer '{ceiling.observer}'")

    def observer_map(self) -> dict[str, ObserverSpec]:
        return {observer.name: observer for observer in self.observers}

    def evidence_map(self) -> dict[str, VisibleEvidenceSpec]:
        return {evidence.name: evidence for evidence in self.evidence}

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "latent_dim": self.latent_dim,
            "observers": [observer.to_dict() for observer in self.observers],
            "evidence": [item.to_dict() for item in self.evidence],
            "ceilings": [item.to_dict() for item in self.ceilings],
            "constraints": [item.to_dict() for item in self.constraints],
            "assumptions": self.assumptions.to_dict(),
            "goals": [goal.to_dict() for goal in self.goals],
            "description": self.description,
            "provenance": self.provenance,
            "tags": list(self.tags),
        }

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "ProblemSpec":
        return cls(
            name=str(data["name"]),
            latent_dim=None if data.get("latent_dim") is None else int(data["latent_dim"]),
            observers=tuple(ObserverSpec.from_dict(item) for item in data.get("observers", [])),
            evidence=tuple(VisibleEvidenceSpec.from_dict(item) for item in data.get("evidence", [])),
            ceilings=tuple(CeilingSpec.from_dict(item) for item in data.get("ceilings", [])),
            constraints=tuple(ConstraintSpec.from_dict(item) for item in data.get("constraints", [])),
            assumptions=AssumptionLedger.from_dict(dict(data.get("assumptions", {}))),
            goals=tuple(GoalSpec.from_dict(item) for item in data.get("goals", [])),
            description=str(data.get("description", "")),
            provenance=str(data.get("provenance", "")),
            tags=tuple(str(item) for item in data.get("tags", [])),
        )

    @classmethod
    def from_json(cls, path: str | Path) -> "ProblemSpec":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
