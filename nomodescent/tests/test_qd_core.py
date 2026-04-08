from __future__ import annotations

import math

import numpy as np

from nomodescent import (
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    classify_qd_relation,
    common_refinement_test,
    false_collapse_diagnostic,
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
    return ProblemSpec(
        name="bell",
        latent_dim=4,
        observers=observers,
        evidence=evidence,
        goals=(GoalSpec(kind="common_completion"),),
    )


def test_common_refinement_constructs_minimal_row_space_refinement() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 0.0, 0.0]])
    observer_b = ObserverSpec(name="b", matrix=[[0.0, 1.0, 0.0]])
    result = common_refinement_test((observer_a, observer_b))
    assert result.classification == "exact_common_refinement"
    assert result.common_refinement is not None
    assert result.common_refinement.matrix.shape == (2, 3)
    assert np.allclose(result.common_refinement.matrix, np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))
    assert set(result.factor_maps) == {"a", "b"}


def test_common_refinement_rejects_invalid_candidate() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 0.0, 0.0]])
    observer_b = ObserverSpec(name="b", matrix=[[0.0, 1.0, 0.0]])
    candidate = ObserverSpec(name="bad", matrix=[[1.0, 0.0, 0.0]])
    result = common_refinement_test((observer_a, observer_b), candidate=candidate)
    assert result.classification == "candidate_not_common_refinement"
    assert any(cert.kind == "candidate_not_common_refinement" for cert in result.certificates)


def test_classify_qd_relation_recovers_factorisation_direction() -> None:
    coarse = ObserverSpec(name="coarse", matrix=[[1.0, 0.0]])
    fine = ObserverSpec(name="fine", matrix=[[1.0, 0.0], [0.0, 1.0]])
    result = classify_qd_relation(coarse, fine)
    assert result.classification == "exact_factorisation"
    assert result.common_refinement is not None
    assert result.common_refinement.name == "fine"


def test_classify_qd_relation_non_nested_includes_common_refinement() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 0.0, 0.0]])
    observer_b = ObserverSpec(name="b", matrix=[[0.0, 1.0, 0.0]])
    result = classify_qd_relation(observer_a, observer_b)
    assert result.classification == "non_nested_observers"
    assert result.common_refinement is not None
    assert result.common_refinement.matrix.shape == (2, 3)


def test_false_collapse_detected_when_coarse_method_is_non_obstructive_but_full_problem_is_incompatible() -> None:
    problem = _bell_problem(delta=0.22, rho=0.6)
    rho = 0.6
    arcsine_margin = math.pi - abs(math.asin(rho) + math.asin(rho) + math.asin(rho) - math.asin(-rho))
    result = false_collapse_diagnostic(
        fine_problem=problem,
        coarse_summaries={
            "00": np.array([[rho, rho], [rho, -rho]], dtype=float),
            "01": np.array([[rho, rho], [rho, -rho]], dtype=float),
        },
        coarse_summary_label="normalized_correlator_matrix",
        coarse_certified_compatible=arcsine_margin > 0.0,
    )
    assert result.classification == "false_collapse_detected"
    assert result.fine_result.classification == "incompatible_by_linear_inconsistency"


def test_false_collapse_not_detected_when_fine_problem_is_compatible() -> None:
    problem = _bell_problem(delta=0.0, rho=0.6)
    rho = 0.6
    result = false_collapse_diagnostic(
        fine_problem=problem,
        coarse_summaries={
            "left": np.array([[rho, rho], [rho, -rho]], dtype=float),
            "right": np.array([[rho, rho], [rho, -rho]], dtype=float),
        },
        coarse_summary_label="normalized_correlator_matrix",
        coarse_certified_compatible=True,
        completion_kwargs={"allow_approximate_psd_search": True, "psd_search_grid_size": 21},
    )
    assert result.classification == "no_false_collapse"
    assert result.fine_result.classification == "approximate_common_descent"
