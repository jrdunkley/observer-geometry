from __future__ import annotations

import numpy as np

from nomogeo import observed_covariance, visible_precision

from .audit import AuditReport
from .engine import factorisation_test
from .exceptions import DescentInputError
from .results import ObstructionCertificate, RefinementSearchResult
from .specs import ObserverSpec


def minimal_refinement_search(
    H: np.ndarray,
    targets: list[ObserverSpec],
    candidates: list[ObserverSpec],
    objective: str = "min_dimension_then_logdet",
    tol: float = 1e-10,
) -> RefinementSearchResult:
    if objective not in {"min_dimension_then_logdet", "visible_logdet", "covariance_trace"}:
        raise DescentInputError(
            "objective must be 'min_dimension_then_logdet', 'visible_logdet', or 'covariance_trace'"
        )

    scores: dict[str, float] = {}
    dimensions: dict[str, int] = {}
    certificates: list[ObstructionCertificate] = []
    H0 = np.asarray(H, dtype=float)

    for candidate in candidates:
        factorizations = [factorisation_test(target.matrix, candidate.matrix, tol=tol) for target in targets]
        if any(result.classification != "exact_factorisation" for result in factorizations):
            certificates.append(
                ObstructionCertificate(
                    kind="candidate_not_common_refinement",
                    exact=True,
                    summary=f"candidate '{candidate.name}' does not refine every target observer",
                    details={"candidate": candidate.name},
                )
            )
            continue

        phi_candidate = visible_precision(H0, candidate.matrix)
        dimensions[candidate.name] = candidate.matrix.shape[0]
        if objective == "visible_logdet":
            scores[candidate.name] = -float(np.linalg.slogdet(phi_candidate)[1])
        elif objective == "covariance_trace":
            scores[candidate.name] = float(np.trace(observed_covariance(H0, candidate.matrix)))
        else:
            scores[candidate.name] = -float(np.linalg.slogdet(phi_candidate)[1])

    winner = None
    if scores:
        if objective == "min_dimension_then_logdet":
            winner = sorted(scores.items(), key=lambda item: (dimensions[item[0]], item[1], item[0]))[0][0]
        else:
            winner = sorted(scores.items(), key=lambda item: (item[1], item[0]))[0][0]

    audit = AuditReport(
        exact_assumptions=("latent precision H is treated as authoritative", "candidate family is finite and explicit"),
        approximate_assumptions=(),
        authoritative_inputs=("latent precision H", "target observer family", "candidate observer family"),
        residuals={
            **{f"score_{name}": score for name, score in scores.items()},
            **{f"dimension_{name}": float(dim) for name, dim in dimensions.items()},
        },
        theorem_layer="exact finite-family refinement ranking",
        falsification_route=(
            "Add a candidate observer with a strictly smaller audited score.",
            "Show that the reported winner does not factor every target observer.",
        ),
    )
    return RefinementSearchResult(
        winner=winner,
        objective=objective,
        scores=scores,
        certificates=tuple(certificates),
        audit=audit,
        details={"target_count": len(targets), "candidate_count": len(candidates)},
    )
