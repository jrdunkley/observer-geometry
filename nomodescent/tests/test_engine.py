from __future__ import annotations

import numpy as np

from nomodescent import (
    AssumptionEntry,
    AssumptionLedger,
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    classify_relation,
    common_descent_test,
    factorisation_test,
    staged_descent_check,
)


def _bell_problem(delta: float, rho: float) -> ProblemSpec:
    observers = (
        ObserverSpec(name="00", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="01", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
        ObserverSpec(name="10", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="11", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )
    sigma_00 = np.array([[1.0, rho], [rho, 1.0]])
    sigma_01 = np.array([[1.0 + delta, rho * np.sqrt(1.0 + delta)], [rho * np.sqrt(1.0 + delta), 1.0]])
    sigma_10 = np.array([[1.0, rho], [rho, 1.0]])
    sigma_11 = np.array([[1.0, -rho], [-rho, 1.0]])
    evidence = (
        VisibleEvidenceSpec(name="s00", observer="00", kind="covariance", matrix=sigma_00),
        VisibleEvidenceSpec(name="s01", observer="01", kind="covariance", matrix=sigma_01),
        VisibleEvidenceSpec(name="s10", observer="10", kind="covariance", matrix=sigma_10),
        VisibleEvidenceSpec(name="s11", observer="11", kind="covariance", matrix=sigma_11),
    )
    assumptions = AssumptionLedger(
        entries=(AssumptionEntry(label="gaussian", statement="Gaussian common descent model", exact=True),)
    )
    return ProblemSpec(
        name="bell",
        latent_dim=4,
        observers=observers,
        evidence=evidence,
        assumptions=assumptions,
        goals=(GoalSpec(kind="common_completion"),),
    )


def test_exact_factorisation_detection() -> None:
    coarse = ObserverSpec(name="coarse", matrix=[[1.0, 0.0]])
    fine = ObserverSpec(name="fine", matrix=[[1.0, 0.0], [0.0, 1.0]])
    result = factorisation_test(coarse, fine)
    assert result.classification == "exact_factorisation"
    assert result.factor_map is not None
    assert np.allclose(result.factor_map, np.array([[1.0, 0.0]]))


def test_exact_nonfactorisation_detection() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 1.0]])
    observer_b = ObserverSpec(name="b", matrix=[[1.0, 0.0]])
    result = factorisation_test(observer_a, observer_b)
    assert result.classification == "factorisation_failure"
    assert result.certificates[0].kind == "factorisation_failure"


def test_exact_relation_classification_non_nested() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 0.0, 0.0]])
    observer_b = ObserverSpec(name="b", matrix=[[0.0, 1.0, 0.0]])
    result = classify_relation(observer_a, observer_b)
    assert result.classification == "non_nested_observers"


def test_exact_tower_agreement() -> None:
    H = np.array([[3.0, 0.4, 0.1], [0.4, 2.7, 0.3], [0.1, 0.3, 1.9]])
    chain = [
        ObserverSpec(name="mid", matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        ObserverSpec(name="coarse", matrix=[[1.0, 0.0]]),
    ]
    result = staged_descent_check(H, chain)
    assert result.classification == "exact_tower_agreement"
    assert result.residuals["tower_residual"] < 1e-12


def test_exact_common_completion_on_known_compatible_example() -> None:
    problem = ProblemSpec(
        name="compatible",
        latent_dim=2,
        observers=(ObserverSpec(name="full", matrix=[[1.0, 0.0], [0.0, 1.0]]),),
        evidence=(VisibleEvidenceSpec(name="cov", observer="full", kind="covariance", matrix=[[2.0, 0.3], [0.3, 1.5]]),),
        goals=(GoalSpec(kind="common_completion"),),
    )
    result = common_descent_test(problem)
    assert result.classification == "exact_common_descent"
    assert result.common_covariance is not None
    assert np.allclose(result.common_covariance, np.array([[2.0, 0.3], [0.3, 1.5]]))


def test_exact_obstruction_detection_by_linear_inconsistency() -> None:
    result = common_descent_test(_bell_problem(delta=0.22, rho=0.6))
    assert result.classification == "incompatible_by_linear_inconsistency"
    assert any(cert.kind == "repeated_marginal_inconsistency" for cert in result.certificates)


def test_exact_obstruction_detection_by_unique_psd_failure() -> None:
    problem = ProblemSpec(
        name="non_psd",
        latent_dim=2,
        observers=(ObserverSpec(name="full", matrix=[[1.0, 0.0], [0.0, 1.0]]),),
        evidence=(VisibleEvidenceSpec(name="cov", observer="full", kind="covariance", matrix=[[1.0, 2.0], [2.0, 1.0]]),),
        goals=(GoalSpec(kind="common_completion"),),
    )
    result = common_descent_test(problem)
    assert result.classification == "incompatible_by_psd_obstruction"
    assert result.exact


def test_exact_vs_approximate_boundary_handling() -> None:
    problem = _bell_problem(delta=0.0, rho=0.82)
    exact_boundary = common_descent_test(problem, allow_approximate_psd_search=False)
    approx_boundary = common_descent_test(problem, allow_approximate_psd_search=True, psd_search_grid_size=21)
    assert exact_boundary.classification == "underdetermined_affine_family"
    assert exact_boundary.exact
    assert approx_boundary.classification == "incompatible_by_approximate_psd_search"
    assert not approx_boundary.exact


def test_audit_object_completeness() -> None:
    result = common_descent_test(_bell_problem(delta=0.22, rho=0.6))
    audit = result.audit
    assert audit.theorem_layer
    assert audit.authoritative_inputs
    assert audit.falsification_route
    assert "linear_residual" in audit.residuals
