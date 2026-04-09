from __future__ import annotations

from .audit import EvidenceAudit
from .exceptions import EvidenceInputError
from .results import SuggestionResult
from .specs import EvidenceBundle, ProtocolObservation


def infer_observer_candidates(
    bundle: EvidenceBundle,
    protocol_name: str,
    *,
    top_k: int | None = None,
) -> SuggestionResult:
    protocol = _protocol_map(bundle)[protocol_name]
    candidates = [item for item in bundle.observer_hypotheses if item.protocol_name == protocol_name]
    if not candidates:
        raise EvidenceInputError(f"no observer hypotheses available for protocol '{protocol_name}'")

    scores: dict[str, float] = {}
    details: dict[str, object] = {}
    protocol_facts = set(protocol.facts)
    for candidate in candidates:
        feature_set = set(candidate.features)
        matched = len(protocol_facts & feature_set)
        missing = len(protocol_facts - feature_set)
        extra = len(feature_set - protocol_facts)
        score = float(matched - 0.5 * missing - 0.1 * extra)
        scores[candidate.name] = score
        details[candidate.name] = {
            "matched_features": sorted(protocol_facts & feature_set),
            "missing_facts": sorted(protocol_facts - feature_set),
            "extra_features": sorted(feature_set - protocol_facts),
        }

    ordered = tuple(name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0])))
    if top_k is not None:
        ordered = ordered[:top_k]

    assumptions = tuple(
        f"observer hypothesis '{candidate.name}' encodes protocol '{protocol_name}' via features {sorted(candidate.features)}"
        for candidate in candidates
    )
    audit = EvidenceAudit(
        exact_items=(protocol.name,),
        inferred_items=tuple(candidate.name for candidate in candidates),
        ambiguous_items=tuple(name for name in ordered[1:]),
        authoritative_items=(protocol.name,),
        load_bearing_items=(protocol.name,),
        load_bearing_exact_items=(protocol.name,),
        load_bearing_inferred_items=tuple(),
        load_bearing_ambiguous_items=tuple(name for name in ordered[1:]),
        load_bearing_assumptions=assumptions,
        unresolved_items=tuple(name for name in ordered[1:]),
        residuals={"score_gap_top2": _top_gap(scores)},
        theorem_layer="finite-family evidence suggestion; not a theorem-grade descent conclusion",
        falsification_route=(
            "do not treat the top-ranked candidate as authoritative without an explicit observer-level justification",
            "supply an explicit observer matrix from the source and bypass finite-family suggestion",
            "add protocol facts that distinguish between the tied candidate hypotheses",
        ),
        ambiguity_collapse_route=(
            "record one more exact protocol fact that only one candidate family satisfies",
        ),
        notes=(
            "candidate ranking is deterministic and finite-family only",
            "the score is a lightweight feature-overlap heuristic, not a scientific selector",
        ),
    )
    return SuggestionResult(
        target=protocol_name,
        ranked_candidates=ordered,
        scores=scores,
        assumptions=assumptions,
        audit=audit,
        details=details,
    )


def _protocol_map(bundle: EvidenceBundle) -> dict[str, ProtocolObservation]:
    return {item.name: item for item in bundle.protocol_observations}


def _top_gap(scores: dict[str, float]) -> float:
    ordered = sorted(scores.values(), reverse=True)
    if len(ordered) < 2:
        return float("inf")
    return float(ordered[0] - ordered[1])
