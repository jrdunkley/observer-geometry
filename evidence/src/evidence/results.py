from __future__ import annotations

from dataclasses import dataclass, field

from nomodescent import ProblemSpec

from .audit import EvidenceAudit


@dataclass(frozen=True)
class AssemblyResult:
    classification: str
    exact: bool
    problem_spec: ProblemSpec | None
    selected_observers: tuple[str, ...]
    unresolved_observer_groups: tuple[str, ...]
    required_human_decisions: tuple[str, ...]
    audit: EvidenceAudit
    details: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class SuggestionResult:
    target: str
    ranked_candidates: tuple[str, ...]
    scores: dict[str, float]
    assumptions: tuple[str, ...]
    audit: EvidenceAudit
    details: dict[str, object] = field(default_factory=dict)
