from __future__ import annotations

from dataclasses import asdict

import numpy as np
import scipy.linalg as la
import scipy.optimize as opt

from nomogeo import observed_covariance, visible_precision

from ._utils import as_array, devectorize_symmetric, min_eig, null_space, operator_constraint_rows, vectorize_symmetric
from .audit import AuditReport
from .exceptions import DescentInputError
from .results import DescentResult, ObstructionCertificate
from .specs import ObserverSpec, ProblemSpec, VisibleEvidenceSpec


def _observer_matrix(observer: ObserverSpec | np.ndarray) -> np.ndarray:
    if isinstance(observer, ObserverSpec):
        return observer.matrix
    matrix = as_array(observer)
    if matrix.ndim != 2:
        raise DescentInputError("observer must be a 2D array")
    return matrix


def _factor_map(C1: np.ndarray, C2: np.ndarray) -> tuple[np.ndarray, float]:
    gram = C2 @ C2.T
    D = C1 @ C2.T @ np.linalg.inv(gram)
    residual = float(np.linalg.norm(C1 - D @ C2, ord=np.inf))
    return D, residual


def _falsification_route_for_factorisation(exists: bool, residual: float) -> tuple[str, ...]:
    if exists:
        return (
            "Recompute the factor map D and verify C1 - D C2 is below tolerance.",
            "Provide a vector in ker(C2) on which C1 acts nontrivially to break the factorisation.",
        )
    return (
        "Exhibit an exact factor map D with residual below tolerance.",
        "Show that every vector in ker(C2) is also annihilated by C1.",
    )


def factorisation_test(
    coarse: ObserverSpec | np.ndarray,
    fine: ObserverSpec | np.ndarray,
    tol: float = 1e-10,
) -> DescentResult:
    C1 = _observer_matrix(coarse)
    C2 = _observer_matrix(fine)
    if C1.shape[1] != C2.shape[1]:
        raise DescentInputError("observer matrices must act on the same latent dimension")
    if np.linalg.matrix_rank(C2) != C2.shape[0]:
        raise DescentInputError("fine observer must have full row rank")

    D, residual = _factor_map(C1, C2)
    exists = residual <= tol
    classification = "exact_factorisation" if exists else "factorisation_failure"
    certificates: list[ObstructionCertificate] = []
    if not exists:
        null = null_space(C2)
        null_residual = float(np.linalg.norm(C1 @ null, ord=np.inf)) if null.size else residual
        certificates.append(
            ObstructionCertificate(
                kind="factorisation_failure",
                exact=True,
                summary="coarse observer does not factor through fine observer",
                details={"factorisation_residual": residual, "kernel_obstruction": null_residual},
            )
        )

    audit = AuditReport(
        exact_assumptions=("observers are linear", "full-row-rank observers are treated exactly"),
        approximate_assumptions=(),
        authoritative_inputs=("coarse observer matrix", "fine observer matrix"),
        residuals={"factorisation_residual": residual},
        theorem_layer="exact linear observer factorisation",
        falsification_route=_falsification_route_for_factorisation(exists, residual),
    )
    return DescentResult(
        classification=classification,
        exact=True,
        factor_map=D if exists else None,
        residuals={"factorisation_residual": residual},
        certificates=tuple(certificates),
        common_covariance=None,
        common_precision=None,
        audit=audit,
        details={"coarse_shape": C1.shape, "fine_shape": C2.shape},
    )


def classify_relation(
    observer_a: ObserverSpec | np.ndarray,
    observer_b: ObserverSpec | np.ndarray,
    tol: float = 1e-10,
) -> DescentResult:
    a_factors_b = factorisation_test(observer_a, observer_b, tol=tol)
    b_factors_a = factorisation_test(observer_b, observer_a, tol=tol)

    if a_factors_b.classification == "exact_factorisation" and b_factors_a.classification == "exact_factorisation":
        classification = "equivalent_observers"
    elif a_factors_b.classification == "exact_factorisation":
        classification = "a_is_coarsening_of_b"
    elif b_factors_a.classification == "exact_factorisation":
        classification = "b_is_coarsening_of_a"
    else:
        classification = "non_nested_observers"

    certificates: list[ObstructionCertificate] = []
    if classification == "non_nested_observers":
        certificates.extend(a_factors_b.certificates)
        certificates.extend(b_factors_a.certificates)

    audit = AuditReport(
        exact_assumptions=("observer relation test is purely linear",),
        approximate_assumptions=(),
        authoritative_inputs=("observer A matrix", "observer B matrix"),
        residuals={
            "a_through_b_residual": a_factors_b.residuals["factorisation_residual"],
            "b_through_a_residual": b_factors_a.residuals["factorisation_residual"],
        },
        theorem_layer="exact linear observer relation classification",
        falsification_route=(
            "Show an exact factorisation in one direction if a non-nested classification is wrong.",
            "Show both directions factor if an equivalent classification is claimed.",
        ),
    )
    return DescentResult(
        classification=classification,
        exact=True,
        factor_map=a_factors_b.factor_map if classification == "a_is_coarsening_of_b" else b_factors_a.factor_map if classification == "b_is_coarsening_of_a" else None,
        residuals=asdict(audit)["residuals"],
        certificates=tuple(certificates),
        common_covariance=None,
        common_precision=None,
        audit=audit,
        details={},
    )


def staged_descent_check(H: np.ndarray, observer_chain: list[ObserverSpec | np.ndarray], tol: float = 1e-10) -> DescentResult:
    if len(observer_chain) < 2:
        raise DescentInputError("observer_chain must contain at least two observers")
    H0 = as_array(H)
    direct = H0
    composed = _observer_matrix(observer_chain[0])
    for observer in observer_chain[1:]:
        composed = _observer_matrix(observer) @ composed
    direct_phi = visible_precision(H0, composed)

    staged = H0
    for observer in observer_chain:
        staged = visible_precision(staged, _observer_matrix(observer))

    residual = float(np.linalg.norm(direct_phi - staged, ord=np.inf))
    classification = "exact_tower_agreement" if residual <= tol else "tower_disagreement"
    certificates = ()
    if residual > tol:
        certificates = (
            ObstructionCertificate(
                kind="tower_disagreement",
                exact=True,
                summary="staged descent does not match direct quotient descent",
                details={"tower_residual": residual},
            ),
        )
    audit = AuditReport(
        exact_assumptions=("latent object is SPD", "observers are surjective linear maps"),
        approximate_assumptions=(),
        authoritative_inputs=("latent precision H", "observer chain"),
        residuals={"tower_residual": residual},
        theorem_layer="exact quotient tower law",
        falsification_route=(
            "Recompute direct and staged visible precisions and compare them entrywise.",
            "Provide a staged chain whose residual exceeds tolerance if tower agreement is falsely claimed.",
        ),
    )
    return DescentResult(
        classification=classification,
        exact=True,
        factor_map=None,
        residuals={"tower_residual": residual},
        certificates=certificates,
        common_covariance=None,
        common_precision=None,
        audit=audit,
        details={"direct_dimension": int(direct_phi.shape[0])},
    )


def _evidence_covariance(evidence: VisibleEvidenceSpec) -> np.ndarray:
    if evidence.kind == "covariance":
        return evidence.matrix
    return np.linalg.inv(evidence.matrix)


def _build_completion_system(problem: ProblemSpec) -> tuple[np.ndarray, np.ndarray]:
    if problem.latent_dim is None:
        raise DescentInputError("latent_dim is required for common completion problems")
    observer_map = problem.observer_map()
    rows = []
    rhs = []
    for evidence in problem.evidence:
        C = observer_map[evidence.observer].matrix
        operator_rows, _coords = operator_constraint_rows(C, problem.latent_dim)
        rows.append(operator_rows)
        rhs.append(_evidence_covariance(evidence).reshape(-1))
    if not rows:
        raise DescentInputError("common completion requires at least one visible evidence object")
    return np.vstack(rows), np.concatenate(rhs)


def _repeated_marginal_certificate(problem: ProblemSpec, tol: float) -> ObstructionCertificate | None:
    observer_map = problem.observer_map()
    seen: dict[str, float] = {}
    for evidence in problem.evidence:
        cov = _evidence_covariance(evidence)
        C = observer_map[evidence.observer].matrix
        for i in range(C.shape[0]):
            key = np.array2string(C[i], precision=12, separator=",")
            variance = float(cov[i, i])
            if key in seen and abs(seen[key] - variance) > tol:
                return ObstructionCertificate(
                    kind="repeated_marginal_inconsistency",
                    exact=True,
                    summary="the same observed linear functional carries incompatible marginal variance data",
                    details={"row": key, "variance_1": seen[key], "variance_2": variance},
                )
            seen[key] = variance
    return None


def _psd_search_affine(s0: np.ndarray, null: np.ndarray, dim: int, radius: float, grid_size: int) -> tuple[np.ndarray, float]:
    if null.size == 0:
        sigma = devectorize_symmetric(s0, dim)
        return sigma, min_eig(sigma)
    if null.shape[1] > 2:
        raise DescentInputError("approximate PSD search only supports affine nullity up to 2 in v0.2")

    def objective(z: np.ndarray) -> float:
        sigma = devectorize_symmetric(s0 + null @ z, dim)
        return -min_eig(sigma)

    bounds = [(-radius, radius)] * null.shape[1]
    best_z = np.zeros(null.shape[1], dtype=float)
    best_val = objective(best_z)
    grid = np.linspace(-radius, radius, grid_size)
    if null.shape[1] == 1:
        starts = [np.array([value], dtype=float) for value in grid]
    else:
        starts = [np.array([x, y], dtype=float) for x in grid for y in grid]
    for start in starts:
        result = opt.minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        value = float(result.fun)
        if value < best_val:
            best_val = value
            best_z = np.asarray(result.x, dtype=float)
    sigma = devectorize_symmetric(s0 + null @ best_z, dim)
    return sigma, min_eig(sigma)


def common_gaussian_completion(
    problem: ProblemSpec,
    tol: float = 1e-10,
    allow_approximate_psd_search: bool = False,
    psd_search_radius: float = 2.0,
    psd_search_grid_size: int = 41,
) -> DescentResult:
    A, b = _build_completion_system(problem)
    solution, *_ = np.linalg.lstsq(A, b, rcond=None)
    linear_residual = float(np.linalg.norm(A @ solution - b, ord=np.inf))
    certificates: list[ObstructionCertificate] = []

    repeated = _repeated_marginal_certificate(problem, tol)
    if repeated is not None:
        certificates.append(repeated)

    if linear_residual > tol:
        certificates.append(
            ObstructionCertificate(
                kind="linear_constraint_inconsistency",
                exact=True,
                summary="no symmetric covariance satisfies the visible linear constraints",
                details={"linear_residual": linear_residual},
            )
        )
        audit = AuditReport(
            exact_assumptions=problem.assumptions.exact_entries(),
            approximate_assumptions=problem.assumptions.approximate_entries(),
            authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
            residuals={"linear_residual": linear_residual},
            theorem_layer="exact Gaussian linear compatibility test",
            falsification_route=(
                "Exhibit a symmetric covariance matching all visible constraints exactly.",
                "Re-encode the evidence with distinct observers if repeated-marginal mismatch was a modelling error.",
            ),
            notes=("classification is exact because the linear constraints already fail",),
        )
        return DescentResult(
            classification="incompatible_by_linear_inconsistency",
            exact=True,
            factor_map=None,
            residuals={"linear_residual": linear_residual},
            certificates=tuple(certificates),
            common_covariance=None,
            common_precision=None,
            audit=audit,
            details={"nullity": int(A.shape[1] - np.linalg.matrix_rank(A))},
        )

    sigma = devectorize_symmetric(solution, problem.latent_dim)
    null = null_space(A)
    psd_margin = min_eig(sigma)

    if null.size == 0:
        if psd_margin > tol:
            precision = np.linalg.inv(sigma)
            residuals = {"linear_residual": linear_residual, "psd_margin": psd_margin}
            audit = AuditReport(
                exact_assumptions=problem.assumptions.exact_entries(),
                approximate_assumptions=problem.assumptions.approximate_entries(),
                authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
                residuals=residuals,
                theorem_layer="exact Gaussian common completion in finite dimension",
                falsification_route=(
                    "Show that the reconstructed covariance fails one of the visible observation constraints.",
                    "Show that the reconstructed covariance is not positive definite.",
                ),
            )
            return DescentResult(
                classification="exact_common_descent",
                exact=True,
                factor_map=None,
                residuals=residuals,
                certificates=tuple(certificates),
                common_covariance=sigma,
                common_precision=precision,
                audit=audit,
                details={"nullity": 0},
            )
        certificates.append(
            ObstructionCertificate(
                kind="psd_obstruction",
                exact=True,
                summary="the unique affine completion is not positive definite",
                details={"psd_margin": psd_margin},
            )
        )
        audit = AuditReport(
            exact_assumptions=problem.assumptions.exact_entries(),
            approximate_assumptions=problem.assumptions.approximate_entries(),
            authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
            residuals={"linear_residual": linear_residual, "psd_margin": psd_margin},
            theorem_layer="exact Gaussian common completion in finite dimension",
            falsification_route=(
                "Provide a positive definite covariance satisfying the same linear constraints.",
                "Show that the linear system admits a different affine completion family.",
            ),
            notes=("classification is exact because the affine completion is unique and non-PSD",),
        )
        return DescentResult(
            classification="incompatible_by_psd_obstruction",
            exact=True,
            factor_map=None,
            residuals={"linear_residual": linear_residual, "psd_margin": psd_margin},
            certificates=tuple(certificates),
            common_covariance=None,
            common_precision=None,
            audit=audit,
            details={"nullity": 0},
        )

    if allow_approximate_psd_search:
        sigma_best, best_margin = _psd_search_affine(solution, null, problem.latent_dim, psd_search_radius, psd_search_grid_size)
        details = {"nullity": int(null.shape[1]), "psd_search_radius": psd_search_radius, "psd_search_grid_size": psd_search_grid_size}
        if best_margin > tol:
            precision = np.linalg.inv(sigma_best)
            audit = AuditReport(
                exact_assumptions=problem.assumptions.exact_entries(),
                approximate_assumptions=problem.assumptions.approximate_entries() + ("deterministic affine PSD search",),
                authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
                residuals={"linear_residual": linear_residual, "best_psd_margin": best_margin},
                theorem_layer="audited approximate PSD-feasibility search on an exact affine family",
                falsification_route=(
                    "Produce a better affine parameter with a larger PSD margin if the current search is wrong.",
                    "Show that no positive definite point exists in the affine family.",
                ),
            )
            return DescentResult(
                classification="approximate_common_descent",
                exact=False,
                factor_map=None,
                residuals={"linear_residual": linear_residual, "best_psd_margin": best_margin},
                certificates=tuple(certificates),
                common_covariance=sigma_best,
                common_precision=precision,
                audit=audit,
                details=details,
            )
        certificates.append(
            ObstructionCertificate(
                kind="approximate_psd_search_failure",
                exact=False,
                summary="deterministic affine PSD search found no positive definite completion",
                details={"best_psd_margin": best_margin, **details},
            )
        )
        audit = AuditReport(
            exact_assumptions=problem.assumptions.exact_entries(),
            approximate_assumptions=problem.assumptions.approximate_entries() + ("deterministic affine PSD search",),
            authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
            residuals={"linear_residual": linear_residual, "best_psd_margin": best_margin},
            theorem_layer="audited approximate PSD-feasibility search on an exact affine family",
            falsification_route=(
                "Provide a positive definite point in the affine completion family.",
                "Expand the search radius or resolution and check whether the PSD margin changes sign.",
            ),
        )
        return DescentResult(
            classification="incompatible_by_approximate_psd_search",
            exact=False,
            factor_map=None,
            residuals={"linear_residual": linear_residual, "best_psd_margin": best_margin},
            certificates=tuple(certificates),
            common_covariance=None,
            common_precision=None,
            audit=audit,
            details=details,
        )

    certificates.append(
        ObstructionCertificate(
            kind="affine_family_underdetermined",
            exact=True,
            summary="linear constraints define an affine completion family that is not resolved exactly in v0.2",
            details={"nullity": int(null.shape[1]), "psd_margin_at_lstsq_point": psd_margin},
        )
    )
    audit = AuditReport(
        exact_assumptions=problem.assumptions.exact_entries(),
        approximate_assumptions=problem.assumptions.approximate_entries(),
        authoritative_inputs=tuple(item.name for item in problem.evidence if item.authoritative),
        residuals={"linear_residual": linear_residual, "psd_margin_at_lstsq_point": psd_margin},
        theorem_layer="exact affine Gaussian completion family detection",
        falsification_route=(
            "Either produce a positive definite member of the affine family or prove none exists.",
            "Reduce the nullity by supplying additional exact observer evidence.",
        ),
        notes=("current engine deliberately returns 'underdetermined' instead of guessing across a higher-dimensional affine PSD family",),
    )
    return DescentResult(
        classification="underdetermined_affine_family",
        exact=True,
        factor_map=None,
        residuals={"linear_residual": linear_residual, "psd_margin_at_lstsq_point": psd_margin},
        certificates=tuple(certificates),
        common_covariance=None,
        common_precision=None,
        audit=audit,
        details={"nullity": int(null.shape[1])},
    )


def common_descent_test(problem: ProblemSpec, **kwargs: object) -> DescentResult:
    return common_gaussian_completion(problem, **kwargs)
