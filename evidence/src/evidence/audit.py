from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class EvidenceAudit:
    exact_items: tuple[str, ...]
    inferred_items: tuple[str, ...]
    ambiguous_items: tuple[str, ...]
    authoritative_items: tuple[str, ...]
    load_bearing_items: tuple[str, ...]
    load_bearing_exact_items: tuple[str, ...]
    load_bearing_inferred_items: tuple[str, ...]
    load_bearing_ambiguous_items: tuple[str, ...]
    load_bearing_assumptions: tuple[str, ...]
    unresolved_items: tuple[str, ...]
    residuals: dict[str, float]
    theorem_layer: str
    falsification_route: tuple[str, ...]
    ambiguity_collapse_route: tuple[str, ...]
    notes: tuple[str, ...] = field(default_factory=tuple)
