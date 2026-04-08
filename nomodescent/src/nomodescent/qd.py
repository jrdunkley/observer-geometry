from __future__ import annotations

from typing import Iterable

import numpy as np

from .audit import AuditReport
from .engine import common_descent_test, factorisation_test
from .exceptions import DescentInputError
from .results import CommonRefinementResult, FalseCollapseResult, ObstructionCertificate
from .specs import ObserverSpec, ProblemSpec


def classify_qd_relation(
    observer_a: ObserverSpec | np.ndarray,
    observer_b: ObserverSpec | np.ndarray,
    tol: float = 1e-10,
) -> CommonRefinementResult:
    named = _normalise_observers((("observer_a", observer_a), ("observer_b", observer_b)))
    first = factorisation_test(named[0], named[1], tol=tol)
    second = factorisation_test(named[1], named[0], tol=tol)

    if first.classification == "exact_factorisation" and second.classification == "exact_factorisation":
        classification = "equivalent_observers"
        factor_maps = {named[0].name: np.eye(named[0].matrix.shape[0]), named[1].name: np.eye(named[1].matrix.shape[0])}
        refinement = named[0]
        certificates: tuple[ObstructionCertificate, ...] = tuple()
    elif first.classification == "exact_factorisation":
        classification = "exact_factorisation"
        factor_maps = {named[0].name: first.factor_map if first.factor_map is not None else np.empty((0, 0))}
        refinement = named[1]
        certificates = tuple()
    elif second.classification == "exact_factorisation":
        classification = "exact_factorisation"
        factor_maps = {named[1].name: second.factor_map if second.factor_map is not None else np.empty((0, 0))}
        refinement = named[0]
        certificates = tuple()
    else:
        classification = "non_nested_observers"
        refinement_result = common_refinement_test(named, tol=tol)
        factor_maps = refinement_result.factor_maps
        refinement = refinement_result.common_refinement
        certificates = tuple([*first.certificates, *second.certificates])

    audit = AuditReport(
        exact_assumptions=("observers are linear and act on the same latent dimension",),
        approximate_assumptions=(),
        authoritative_inputs=(named[0].name, named[1].name),
        residuals={
            "a_through_b_residual": first.residuals["factorisation_residual"],
            "b_through_a_residual": second.residuals["factorisation_residual"],
        },
        theorem_layer="exact Quotient Descent observer relation classification",
        falsification_route=(
            "exhibit an exact factor map in the claimed factorisation direction if the relation is wrong",
            "provide a vector-space obstruction showing the claimed refinement does not span both observer row spaces",
        ),
    )
    return CommonRefinementResult(
        classification=classification,
        exact=True,
        common_refinement=refinement,
        factor_maps=factor_maps,
        residuals=dict(audit.residuals),
        certificates=certificates,
        audit=audit,
        details={"observer_names": (named[0].name, named[1].name)},
    )


def common_refinement_test(
    observers: Iterable[ObserverSpec | np.ndarray],
    *,
    candidate: ObserverSpec | np.ndarray | None = None,
    tol: float = 1e-10,
    name: str = "common_refinement",
) -> CommonRefinementResult:
    named = _normalise_observers(tuple((f"observer_{index}", observer) for index, observer in enumerate(observers)))
    if len(named) < 2:
        raise DescentInputError("common_refinement_test requires at least two observers")

    if candidate is None:
        refinement = _canonical_common_refinement(named, name=name)
    else:
        refinement = _normalise_observer(name, candidate)
        _validate_same_latent_dim((*named, refinement))

    factor_maps: dict[str, np.ndarray] = {}
    residuals: dict[str, float] = {}
    certificates: list[ObstructionCertificate] = []
    max_residual = 0.0
    for observer in named:
        result = factorisation_test(observer, refinement, tol=tol)
        residual = result.residuals["factorisation_residual"]
        residuals[f"{observer.name}_through_refinement_residual"] = residual
        max_residual = max(max_residual, residual)
        if result.classification == "exact_factorisation" and result.factor_map is not None:
            factor_maps[observer.name] = result.factor_map
            continue
        certificates.append(
            ObstructionCertificate(
                kind="candidate_not_common_refinement",
                exact=True,
                summary=f"candidate '{refinement.name}' does not refine observer '{observer.name}'",
                details={"observer": observer.name, "residual": residual},
            )
        )

    classification = "exact_common_refinement" if not certificates else "candidate_not_common_refinement"
    audit = AuditReport(
        exact_assumptions=("observers are linear and finite-dimensional",),
        approximate_assumptions=(),
        authoritative_inputs=tuple(observer.name for observer in named),
        residuals={"max_refinement_residual": max_residual, **residuals},
        theorem_layer="exact common refinement in linear observer space",
        falsification_route=(
            "show that one observer row is not contained in the claimed refinement row space",
            "provide a smaller row-space basis if the canonical refinement is claimed to be non-minimal",
        ),
        notes=("when no candidate is supplied, the function constructs the minimal stable row-space refinement by greedy rank growth",),
    )
    return CommonRefinementResult(
        classification=classification,
        exact=True,
        common_refinement=refinement if classification == "exact_common_refinement" else None,
        factor_maps=factor_maps,
        residuals=dict(audit.residuals),
        certificates=tuple(certificates),
        audit=audit,
        details={"observer_count": len(named), "refinement_rank": refinement.matrix.shape[0]},
    )


def false_collapse_diagnostic(
    *,
    fine_problem: ProblemSpec,
    coarse_summaries: dict[str, float | np.ndarray | list[float] | list[list[float]]],
    coarse_summary_label: str,
    coarse_certified_compatible: bool = False,
    summary_tol: float = 1e-10,
    completion_kwargs: dict[str, object] | None = None,
) -> FalseCollapseResult:
    if len(coarse_summaries) < 2 and not coarse_certified_compatible:
        raise DescentInputError(
            "false_collapse_diagnostic requires at least two coarse summaries unless coarse_certified_compatible is supplied"
        )
    completion_kwargs = {} if completion_kwargs is None else dict(completion_kwargs)
    fine_result = common_descent_test(fine_problem, **completion_kwargs)
    max_gap = _max_summary_gap(coarse_summaries) if len(coarse_summaries) >= 2 else 0.0
    coarse_agreement = len(coarse_summaries) >= 2 and max_gap <= summary_tol
    fine_incompatible = fine_result.classification in {
        "incompatible_by_linear_inconsistency",
        "incompatible_by_psd_obstruction",
        "incompatible_by_approximate_psd_search",
    }
    detected = fine_incompatible and (coarse_agreement or coarse_certified_compatible)
    classification = "false_collapse_detected" if detected else "no_false_collapse"

    approximate_assumptions = tuple(fine_result.audit.approximate_assumptions)
    if coarse_certified_compatible:
        approximate_assumptions = approximate_assumptions + ("coarse compatibility certification supplied by caller",)

    audit = AuditReport(
        exact_assumptions=tuple(fine_result.audit.exact_assumptions),
        approximate_assumptions=approximate_assumptions,
        authoritative_inputs=tuple(fine_result.audit.authoritative_inputs) + tuple(sorted(coarse_summaries)),
        residuals={
            "coarse_summary_max_gap": max_gap,
            "summary_tol": summary_tol,
            **fine_result.residuals,
        },
        theorem_layer="false-collapse diagnostic comparing coarse summary agreement with richer Quotient Descent compatibility",
        falsification_route=(
            "show that the coarse summaries are not actually equal or not actually coarse-compatible",
            "show that the fine problem is compatible once the richer visible evidence is respected",
        ),
    )
    return FalseCollapseResult(
        classification=classification,
        exact=fine_result.exact and not approximate_assumptions,
        coarse_summary_label=coarse_summary_label,
        coarse_summary_max_gap=max_gap,
        coarse_summary_agreement=coarse_agreement,
        coarse_certified_compatible=coarse_certified_compatible,
        fine_result=fine_result,
        audit=audit,
        details={"coarse_summary_count": len(coarse_summaries)},
    )


def _normalise_observers(named: Iterable[tuple[str, ObserverSpec | np.ndarray]]) -> tuple[ObserverSpec, ...]:
    observers = tuple(_normalise_observer(name, observer) for name, observer in named)
    _validate_same_latent_dim(observers)
    return observers


def _normalise_observer(name: str, observer: ObserverSpec | np.ndarray) -> ObserverSpec:
    if isinstance(observer, ObserverSpec):
        return observer
    return ObserverSpec(name=name, matrix=np.asarray(observer, dtype=float))


def _validate_same_latent_dim(observers: Iterable[ObserverSpec]) -> None:
    observers = tuple(observers)
    dims = {observer.matrix.shape[1] for observer in observers}
    if len(dims) != 1:
        raise DescentInputError("all observers must act on the same latent dimension")


def _canonical_common_refinement(observers: tuple[ObserverSpec, ...], *, name: str) -> ObserverSpec:
    selected_rows: list[np.ndarray] = []
    current = np.zeros((0, observers[0].matrix.shape[1]), dtype=float)
    for observer in observers:
        for row in observer.matrix:
            candidate = np.vstack([current, row])
            if np.linalg.matrix_rank(candidate, tol=1e-12) > current.shape[0]:
                selected_rows.append(np.asarray(row, dtype=float))
                current = candidate
    if not selected_rows:
        raise DescentInputError("failed to build a common refinement from the supplied observers")
    return ObserverSpec(name=name, matrix=np.vstack(selected_rows))


def _max_summary_gap(coarse_summaries: dict[str, float | np.ndarray | list[float] | list[list[float]]]) -> float:
    values = [(name, _summary_array(name, summary)) for name, summary in sorted(coarse_summaries.items())]
    reference_shape = values[0][1].shape
    max_gap = 0.0
    for _name, array in values[1:]:
        if array.shape != reference_shape:
            raise DescentInputError("all coarse summaries must have the same shape")
    for index, (_name, left) in enumerate(values):
        for _other_name, right in values[index + 1 :]:
            max_gap = max(max_gap, float(np.max(np.abs(left - right))))
    return max_gap


def _summary_array(name: str, summary: float | np.ndarray | list[float] | list[list[float]]) -> np.ndarray:
    array = np.asarray(summary, dtype=float)
    if array.ndim == 0:
        return array.reshape(1)
    if array.ndim in {1, 2}:
        return array
    raise DescentInputError(f"coarse summary '{name}' must be scalar, vector, or matrix valued")
