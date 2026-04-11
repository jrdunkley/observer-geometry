import numpy as np
import pytest

from nomogeo import (
    closure_scores,
    declared_frontier_local_certificate,
    declared_ladder_dimension_cost_intervals,
    exact_branch_hessian,
    general_graph_frontier_hessian,
    weighted_family_frontier_scores,
)
from nomogeo.exceptions import InputValidationError


def _graph_score(family: list[np.ndarray], weights: np.ndarray, x: np.ndarray, mu: float) -> float:
    graph = np.vstack([np.eye(x.shape[1]), x])
    q, _ = np.linalg.qr(graph)
    p = q @ q.T
    visible = 0.0
    leakage = 0.0
    for w, a in zip(weights, family):
        visible += float(w * np.linalg.norm(p @ a @ p, ord="fro") ** 2)
        comm = a @ p - p @ a
        leakage += float(w * 0.5 * np.linalg.norm(comm, ord="fro") ** 2)
    return visible - mu * leakage


def test_weighted_family_frontier_scores_energy_split_and_closure_score_agreement() -> None:
    rng = np.random.default_rng(4101)
    family = []
    for _ in range(4):
        raw = rng.normal(size=(5, 5))
        family.append(0.5 * (raw + raw.T))
    q, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    weights = np.array([0.4, 1.2, 0.7, 0.3])

    weighted = weighted_family_frontier_scores(family, q, weights=weights, mu=0.5)

    assert abs(weighted.energy_split_residual) <= 1e-12
    assert weighted.penalized_score == pytest.approx(weighted.visible_score - 0.5 * weighted.leakage)
    assert "not a noncommuting optimiser" in weighted.metadata.notes[-1]

    unweighted = weighted_family_frontier_scores(family, q)
    existing = closure_scores(np.eye(5), family, q)
    assert unweighted.leakage == pytest.approx(existing.leakage)
    assert unweighted.visible_score == pytest.approx(existing.visible_score)
    assert unweighted.captured_curvature == pytest.approx(existing.total_curvature)


def test_weighted_family_frontier_scores_rejects_bad_inputs() -> None:
    with pytest.raises(InputValidationError):
        weighted_family_frontier_scores([np.eye(3), np.eye(4)], np.eye(3, 1))
    with pytest.raises(InputValidationError):
        weighted_family_frontier_scores([np.eye(3)], np.ones((3, 1)))
    with pytest.raises(InputValidationError):
        weighted_family_frontier_scores([np.eye(3)], np.eye(3, 1), weights=np.array([-1.0]))
    with pytest.raises(InputValidationError):
        weighted_family_frontier_scores([np.eye(3)], np.eye(3, 1), mu=-0.1)


def test_exact_branch_hessian_status_and_graph_chart_sign_convention() -> None:
    family = [np.diag([1.0, 0.0, 3.0])]
    weights = np.array([1.0])
    B = np.array([[1.0], [0.0], [0.0]])

    result = exact_branch_hessian(family, B, weights=weights, mu=0.0)

    assert result.status == "saddle"
    assert np.allclose(result.eigenvalues, [-4.0, 2.0], atol=1e-12)
    assert np.allclose(result.second_variation_operator, -2.0 * result.hessian_contract)

    h = 2e-5
    numeric = []
    for row in range(2):
        x = np.zeros((2, 1))
        x[row, 0] = 1.0
        f0 = _graph_score(family, weights, np.zeros((2, 1)), 0.0)
        second = (_graph_score(family, weights, h * x, 0.0) - 2.0 * f0 + _graph_score(family, weights, -h * x, 0.0)) / (h * h)
        numeric.append(second)
    assert np.allclose(np.sort(numeric), np.sort(np.diag(result.second_variation_operator)), atol=1e-5)


def test_exact_branch_hessian_semidefinite_kernel_is_degenerate() -> None:
    result = exact_branch_hessian([np.diag([1.0, 1.0, 2.0])], np.array([[1.0], [0.0], [0.0]]), mu=0.0)

    assert result.status == "degenerate"
    assert result.nullity == 1


def test_exact_branch_hessian_handles_noncommuting_internal_blocks() -> None:
    a1 = np.array([[1.0, 0.7], [0.7, -0.2]])
    a2 = np.array([[0.3, 1.1], [1.1, 0.4]])
    b1 = np.array([[2.0, -0.4], [-0.4, -1.0]])
    b2 = np.array([[-0.8, 0.9], [0.9, 1.7]])
    family = [
        np.block([[a1, np.zeros((2, 2))], [np.zeros((2, 2)), b1]]),
        np.block([[a2, np.zeros((2, 2))], [np.zeros((2, 2)), b2]]),
    ]
    B = np.eye(4, 2)

    result = exact_branch_hessian(family, B, mu=0.4)

    assert result.off_block_norm <= 1e-12
    assert result.hessian_contract.shape == (4, 4)
    assert result.min_eigenvalue > 0.0
    assert result.status == "strict_max"
    assert np.linalg.norm(a1 @ a2 - a2 @ a1, ord="fro") > 1e-2
    assert np.linalg.norm(b1 @ b2 - b2 @ b1, ord="fro") > 1e-2


def test_exact_branch_hessian_rejects_non_exact_branch_and_bad_rank() -> None:
    family = [np.array([[1.0, 0.2], [0.2, 2.0]])]
    with pytest.raises(InputValidationError):
        exact_branch_hessian(family, np.array([[1.0], [0.0]]))
    with pytest.raises(InputValidationError):
        exact_branch_hessian([np.eye(2)], np.eye(2))


def test_declared_ladder_dimension_cost_intervals_reports_phase_diagram() -> None:
    result = declared_ladder_dimension_cost_intervals(
        np.array([4.0, 5.8, 6.4]),
        np.array([1.0, 2.0, 4.0]),
    )

    assert result.winner_at_zero == (2,)
    assert result.interval_nonempty.tolist() == [True, True, True]
    assert result.interval_lower[2] == pytest.approx(0.0)
    assert result.interval_upper[2] == pytest.approx(0.3)
    assert result.interval_lower[1] == pytest.approx(0.3)
    assert result.interval_upper[1] == pytest.approx(1.8)
    assert result.interval_lower[0] == pytest.approx(1.8)
    assert np.isinf(result.interval_upper[0])
    assert "not a global" in result.metadata.notes[-1]


def test_general_graph_frontier_hessian_recovers_exact_branch_hessian() -> None:
    a1 = np.diag([3.0, 2.0])
    d1 = np.diag([0.5, 1.0, 1.5])
    a2 = np.array([[2.5, 0.2], [0.2, 1.7]])
    d2 = np.diag([0.3, 1.3, 2.0])
    family = [
        np.block([[a1, np.zeros((2, 3))], [np.zeros((3, 2)), d1]]),
        np.block([[a2, np.zeros((2, 3))], [np.zeros((3, 2)), d2]]),
    ]
    B = np.eye(5, 2)
    weights = np.array([0.7, 1.2])

    exact = exact_branch_hessian(family, B, weights=weights, mu=0.4)
    graph = general_graph_frontier_hessian(family, B, weights=weights, mu=0.4)

    assert graph.status == "strict_local_max"
    assert graph.stationarity_residual <= 1e-12
    assert np.allclose(graph.hessian_operator, exact.second_variation_operator, atol=1e-12)
    assert np.allclose(np.sort(graph.eigenvalues), np.sort(np.linalg.eigvalsh(exact.second_variation_operator)), atol=1e-12)


def test_general_graph_frontier_hessian_detects_non_exact_stationary_breakpoint() -> None:
    B = np.array([[1.0], [0.0]])
    for e, expected, status in [
        (1.0, -16.0, "strict_local_max"),
        (np.sqrt(3.0), 0.0, "stationary_degenerate"),
        (2.0, 8.0, "strict_local_min"),
    ]:
        family = [
            np.array([[3.0, e], [e, 1.0]]),
            np.array([[3.0, -e], [-e, 1.0]]),
        ]
        result = general_graph_frontier_hessian(family, B, weights=np.array([0.5, 0.5]), mu=0.0)

        assert result.stationarity_residual <= 1e-12
        assert result.hessian_operator.shape == (1, 1)
        assert result.hessian_operator[0, 0] == pytest.approx(expected, abs=1e-12)
        assert result.status == status


def test_declared_frontier_local_certificate_tracks_non_exact_stationary_sign() -> None:
    B = np.array([[1.0], [0.0]])

    max_family = [
        np.array([[3.0, 1.0], [1.0, 1.0]]),
        np.array([[3.0, -1.0], [-1.0, 1.0]]),
    ]
    max_cert = declared_frontier_local_certificate(max_family, B, weights=np.array([0.5, 0.5]), mode="max")
    assert max_cert.certificate_passes
    assert max_cert.certificate_kind == "sufficient_local_max"

    deg_family = [
        np.array([[3.0, np.sqrt(3.0)], [np.sqrt(3.0), 1.0]]),
        np.array([[3.0, -np.sqrt(3.0)], [-np.sqrt(3.0), 1.0]]),
    ]
    deg_cert = declared_frontier_local_certificate(deg_family, B, weights=np.array([0.5, 0.5]), mode="max")
    assert not deg_cert.certificate_passes
    assert deg_cert.certificate_kind == "vacuous"

    min_family = [
        np.array([[3.0, 2.0], [2.0, 1.0]]),
        np.array([[3.0, -2.0], [-2.0, 1.0]]),
    ]
    min_cert = declared_frontier_local_certificate(min_family, B, weights=np.array([0.5, 0.5]), mode="min")
    assert min_cert.certificate_passes
    assert min_cert.certificate_kind == "sufficient_local_min"


def test_general_graph_frontier_hessian_rejects_bad_complement_and_certificate_inputs() -> None:
    family = [np.eye(3)]
    B = np.eye(3, 1)
    with pytest.raises(InputValidationError):
        general_graph_frontier_hessian(family, B, complement_basis=np.eye(3, 1))
    with pytest.raises(InputValidationError):
        declared_frontier_local_certificate(family, B, mode="largest")
    with pytest.raises(InputValidationError):
        declared_frontier_local_certificate(family, B, rho=2.0)
