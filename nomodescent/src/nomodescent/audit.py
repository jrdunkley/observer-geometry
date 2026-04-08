from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AuditReport:
    exact_assumptions: tuple[str, ...]
    approximate_assumptions: tuple[str, ...]
    authoritative_inputs: tuple[str, ...]
    residuals: dict[str, float]
    theorem_layer: str
    falsification_route: tuple[str, ...]
    notes: tuple[str, ...] = field(default_factory=tuple)
