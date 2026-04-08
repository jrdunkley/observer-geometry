from __future__ import annotations

from collections import defaultdict
from typing import Iterable

import numpy as np
from nomodescent import (
    AssumptionEntry,
    AssumptionLedger,
    ConstraintSpec,
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
)

from .audit import EvidenceAudit
from .exceptions import EvidenceAssemblyError
from .results import AssemblyResult
from .specs import EvidenceBundle, MatrixObservation, ObserverHypothesis


def assemble_problem_spec(
    bundle: EvidenceBundle,
    *,
    observer_selection: dict[str, str] | None = None,
    include_non_authoritative: bool = False,
) -> AssemblyResult:
    observer_selection = {} if observer_selection is None else dict(observer_selection)
    unresolved_groups: list[str] = []
    required_decisions: list[str] = []

    selected_hypotheses = _select_observers(bundle, observer_selection, unresolved_groups, required_decisions)
    if unresolved_groups:
        audit = _build_audit(
            bundle,
            exact=False,
            theorem_layer="evidence assembly stopped before ProblemSpec construction",
            unresolved_items=tuple(unresolved_groups),
            load_bearing_assumptions=_load_bearing_assumptions(selected_hypotheses.values()),
            load_bearing_names=_load_bearing_names(bundle, selected_hypotheses.values(), ()),
            notes=("assembly refused because at least one protocol observer family remained unresolved",),
        )
        return AssemblyResult(
            classification="underdetermined_evidence",
            exact=False,
            problem_spec=None,
            selected_observers=tuple(),
            unresolved_observer_groups=tuple(unresolved_groups),
            required_human_decisions=tuple(required_decisions),
            audit=audit,
            details={"observer_selection": observer_selection},
        )

    visible = _select_visible_evidence(bundle, include_non_authoritative=include_non_authoritative)
    if not visible:
        audit = _build_audit(
            bundle,
            exact=False,
            theorem_layer="evidence assembly stopped because no visible object survived authoritative filtering",
            unresolved_items=("visible_object",),
            load_bearing_assumptions=_load_bearing_assumptions(selected_hypotheses.values()),
            load_bearing_names=_load_bearing_names(bundle, selected_hypotheses.values(), ()),
            notes=("assembly refused because the bundle did not contain enough visible matrix evidence",),
        )
        return AssemblyResult(
            classification="insufficient_evidence",
            exact=False,
            problem_spec=None,
            selected_observers=tuple(selected_hypotheses.keys()),
            unresolved_observer_groups=tuple(),
            required_human_decisions=("add at least one authoritative visible covariance or precision matrix",),
            audit=audit,
            details={"observer_selection": observer_selection},
        )

    constraints = tuple(
        ConstraintSpec(
            name=item.name,
            kind=item.kind,
            exact=item.exact,
            payload=item.payload,
            provenance=f"{item.source_ref.source}:{item.source_ref.location}".strip(":"),
            notes=item.notes,
        )
        for item in bundle.constraint_candidates
        if include_non_authoritative or item.extraction.authoritative
    )
    assumptions = tuple(_assumption_entries(bundle, selected_hypotheses.values()))
    resolved_observer_names = {
        protocol_name: hypothesis.name for protocol_name, hypothesis in _protocol_to_hypothesis(bundle, selected_hypotheses).items()
    }
    problem = ProblemSpec(
        name=bundle.name,
        latent_dim=bundle.latent_dim,
        observers=tuple(
            ObserverSpec(
                name=item.name,
                matrix=item.matrix,
                provenance=f"{item.source_ref.source}:{item.source_ref.location}".strip(":"),
                notes=item.notes,
            )
            for item in selected_hypotheses.values()
        ),
        evidence=tuple(
            VisibleEvidenceSpec(
                name=item.name,
                observer=_resolve_observer_name(item.observer_name or "", resolved_observer_names, selected_hypotheses),
                kind=item.matrix_kind,
                matrix=item.matrix,
                exact=item.extraction.extraction_mode == "exact_extraction",
                authoritative=item.extraction.authoritative,
                provenance=f"{item.source_ref.source}:{item.source_ref.location}".strip(":"),
                notes=item.notes,
            )
            for item in visible
        ),
        constraints=constraints,
        assumptions=AssumptionLedger(entries=assumptions),
        goals=(GoalSpec(kind="common_completion"),),
        description=bundle.description,
        provenance=bundle.provenance,
        tags=bundle.tags,
    )
    audit = _build_audit(
        bundle,
        exact=all(item.extraction.extraction_mode == "exact_extraction" for item in visible)
        and all(item.extraction.extraction_mode == "exact_extraction" for item in selected_hypotheses.values()),
        theorem_layer="assembled evidence-to-problem encoding; exact only where extraction_mode is exact_extraction",
        unresolved_items=tuple(),
        load_bearing_assumptions=tuple(entry.statement for entry in assumptions if not entry.exact),
        load_bearing_names=_load_bearing_names(bundle, selected_hypotheses.values(), visible),
        notes=("observer hypotheses were promoted into ProblemSpec observers only after explicit selection",),
    )
    return AssemblyResult(
        classification="assembled_problem_spec",
        exact=len(audit.load_bearing_inferred_items) == 0 and len(audit.load_bearing_ambiguous_items) == 0,
        problem_spec=problem,
        selected_observers=tuple(selected_hypotheses.keys()),
        unresolved_observer_groups=tuple(),
        required_human_decisions=tuple(),
        audit=audit,
        details={"observer_selection": observer_selection},
    )


def _select_observers(
    bundle: EvidenceBundle,
    observer_selection: dict[str, str],
    unresolved_groups: list[str],
    required_decisions: list[str],
) -> dict[str, ObserverHypothesis]:
    grouped: dict[str, list[ObserverHypothesis]] = defaultdict(list)
    for hypothesis in bundle.observer_hypotheses:
        grouped[hypothesis.protocol_name].append(hypothesis)

    selected: dict[str, ObserverHypothesis] = {}
    for protocol in bundle.protocol_observations:
        candidates = grouped.get(protocol.name, [])
        authoritative = [item for item in candidates if item.extraction.authoritative]
        pool = authoritative or candidates
        if not pool:
            unresolved_groups.append(protocol.name)
            required_decisions.append(f"add at least one observer hypothesis for protocol '{protocol.name}'")
            continue
        if protocol.name in observer_selection:
            chosen_name = observer_selection[protocol.name]
            match = next((item for item in pool if item.name == chosen_name), None)
            if match is None:
                raise EvidenceAssemblyError(
                    f"observer selection '{chosen_name}' is not a valid candidate for protocol '{protocol.name}'"
                )
            selected[match.name] = match
            continue
        if len(pool) > 1:
            unresolved_groups.append(protocol.name)
            required_decisions.append(
                f"choose one observer hypothesis for protocol '{protocol.name}' from {[item.name for item in pool]}"
            )
            continue
        selected[pool[0].name] = pool[0]
    return selected


def _select_visible_evidence(
    bundle: EvidenceBundle,
    *,
    include_non_authoritative: bool,
) -> tuple[MatrixObservation, ...]:
    visible: list[MatrixObservation] = []
    seen: dict[tuple[str, str], MatrixObservation] = {}
    for matrix in bundle.matrix_observations:
        if matrix.matrix_role != "visible_object":
            continue
        if not include_non_authoritative and not matrix.extraction.authoritative:
            continue
        key = (matrix.observer_name or "", matrix.matrix_kind)
        if key in seen and not np.allclose(seen[key].matrix, matrix.matrix, atol=1e-10, rtol=1e-10):
            raise EvidenceAssemblyError(
                f"conflicting visible observations remain for observer '{matrix.observer_name}' and kind '{matrix.matrix_kind}'"
            )
        seen[key] = matrix
    visible.extend(seen.values())
    return tuple(visible)


def _protocol_to_hypothesis(
    bundle: EvidenceBundle,
    selected_hypotheses: dict[str, ObserverHypothesis],
) -> dict[str, ObserverHypothesis]:
    mapping: dict[str, ObserverHypothesis] = {}
    hypothesis_by_name = dict(selected_hypotheses)
    for protocol in bundle.protocol_observations:
        match = next((item for item in hypothesis_by_name.values() if item.protocol_name == protocol.name), None)
        if match is not None:
            mapping[protocol.name] = match
    return mapping


def _resolve_observer_name(
    raw_name: str,
    protocol_to_hypothesis: dict[str, str],
    selected_hypotheses: dict[str, ObserverHypothesis],
) -> str:
    if raw_name in selected_hypotheses:
        return raw_name
    if raw_name in protocol_to_hypothesis:
        return protocol_to_hypothesis[raw_name]
    raise EvidenceAssemblyError(f"visible evidence references unresolved observer or protocol '{raw_name}'")


def _assumption_entries(bundle: EvidenceBundle, selected_hypotheses: Iterable[ObserverHypothesis]) -> tuple[AssumptionEntry, ...]:
    entries: list[AssumptionEntry] = []
    for note in bundle.notes:
        entries.append(
            AssumptionEntry(
                label=note.name,
                statement=note.statement,
                exact=note.extraction.extraction_mode == "exact_extraction",
                provenance=f"{note.source_ref.source}:{note.source_ref.location}".strip(":"),
            )
        )
    for hypothesis in selected_hypotheses:
        for index, statement in enumerate(hypothesis.assumption_statements):
            entries.append(
                AssumptionEntry(
                    label=f"{hypothesis.name}_assumption_{index}",
                    statement=statement,
                    exact=hypothesis.extraction.extraction_mode == "exact_extraction",
                    provenance=f"{hypothesis.source_ref.source}:{hypothesis.source_ref.location}".strip(":"),
                )
            )
    return tuple(entries)


def _build_audit(
    bundle: EvidenceBundle,
    *,
    exact: bool,
    theorem_layer: str,
    unresolved_items: tuple[str, ...],
    load_bearing_assumptions: tuple[str, ...],
    load_bearing_names: tuple[str, ...],
    notes: tuple[str, ...],
) -> EvidenceAudit:
    exact_items: list[str] = []
    inferred_items: list[str] = []
    ambiguous_items: list[str] = []
    authoritative_items: list[str] = []

    load_bearing_set = set(load_bearing_names)
    load_bearing_exact_items: list[str] = []
    load_bearing_inferred_items: list[str] = []
    load_bearing_ambiguous_items: list[str] = []

    for item_name, mode, authoritative in _flatten_bundle_items(bundle):
        if authoritative:
            authoritative_items.append(item_name)
        if mode == "exact_extraction":
            exact_items.append(item_name)
            if item_name in load_bearing_set:
                load_bearing_exact_items.append(item_name)
        elif mode == "encoded_inference":
            inferred_items.append(item_name)
            if item_name in load_bearing_set:
                load_bearing_inferred_items.append(item_name)
        else:
            ambiguous_items.append(item_name)
            if item_name in load_bearing_set:
                load_bearing_ambiguous_items.append(item_name)

    return EvidenceAudit(
        exact_items=tuple(sorted(exact_items)),
        inferred_items=tuple(sorted(inferred_items)),
        ambiguous_items=tuple(sorted(ambiguous_items)),
        authoritative_items=tuple(sorted(authoritative_items)),
        load_bearing_items=tuple(sorted(load_bearing_set)),
        load_bearing_exact_items=tuple(sorted(load_bearing_exact_items)),
        load_bearing_inferred_items=tuple(sorted(load_bearing_inferred_items)),
        load_bearing_ambiguous_items=tuple(sorted(load_bearing_ambiguous_items)),
        load_bearing_assumptions=load_bearing_assumptions,
        unresolved_items=unresolved_items,
        residuals={"unresolved_count": float(len(unresolved_items)), "exact_flag": 1.0 if exact else 0.0},
        theorem_layer=theorem_layer,
        falsification_route=(
            "supply an explicit observer matrix or remove a conflicting authoritative item",
            "mark any modelling choice as encoded_inference instead of exact_extraction if it is not directly stated",
        ),
        ambiguity_collapse_route=(
            "add an exact source item that distinguishes between the remaining observer hypotheses",
            "promote only one observer hypothesis per protocol into authoritative downstream use",
        ),
        notes=notes,
    )


def _flatten_bundle_items(bundle: EvidenceBundle) -> tuple[tuple[str, str, bool], ...]:
    triples: list[tuple[str, str, bool]] = []
    for collection in (
        bundle.text_claims,
        bundle.table_observations,
        bundle.matrix_observations,
        bundle.protocol_observations,
        bundle.observer_hypotheses,
        bundle.constraint_candidates,
        bundle.notes,
    ):
        for item in collection:
            triples.append((item.name, item.extraction.extraction_mode, item.extraction.authoritative))
    return tuple(triples)


def _load_bearing_assumptions(selected: Iterable[ObserverHypothesis]) -> tuple[str, ...]:
    statements: list[str] = []
    for hypothesis in selected:
        statements.extend(hypothesis.assumption_statements)
    return tuple(statements)


def _load_bearing_names(
    bundle: EvidenceBundle,
    selected_hypotheses: Iterable[ObserverHypothesis],
    visible: Iterable[MatrixObservation],
) -> tuple[str, ...]:
    names = {item.name for item in selected_hypotheses}
    names.update(item.name for item in visible)
    names.update(item.name for item in bundle.constraint_candidates if item.extraction.authoritative)
    names.update(item.name for item in bundle.notes if item.load_bearing)
    return tuple(sorted(names))
