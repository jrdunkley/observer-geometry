from __future__ import annotations

import itertools

import numpy as np
import pytest

from nomogeo import (
    closure_adapted_observer,
    closure_scores,
    compare_observers,
    leakage_channels,
    observer_from_subspace,
    visible_geometry,
    whitened_perturbation,
)
from nomogeo.exceptions import InputValidationError

from .helpers import orthogonal_matrix, random_spd


def _sqrt_spd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T


def _inv_sqrt_spd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (eigenvectors * (1.0 / np.sqrt(eigenvalues))) @ eigenvectors.T


def _commuting_family(
    rng: np.random.Generator,
    n: int,
    task_count: int,
    *,
    active_dim: int,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray], np.ndarray]:
    H = random_spd(rng, n)
    H_half = _sqrt_spd(H)
    U = orthogonal_matrix(rng, n)
    lambdas = rng.normal(size=(task_count, n))
    lambdas[:, active_dim:] = 0.0
    family = [H_half @ U @ np.diag(row) @ U.T @ H_half for row in lambdas]
    return H, U, family, lambdas


def test_whitened_perturbation_and_observer_normal_form() -> None:
    rng = np.random.default_rng(1101)
    H = random_spd(rng, 6)
    B = orthogonal_matrix(rng, 6)[:, :3]
    raw = rng.normal(size=(6, 6))
    Delta = 0.5 * (raw + raw.T)

    whitened = whitened_perturbation(H, Delta)
    inv_half = _inv_sqrt_spd(H)
    expected_whitened = 0.5 * (inv_half @ Delta @ inv_half + inv_half @ Delta @ inv_half.T)
    assert np.allclose(whitened, whitened.T, atol=1e-10, rtol=1e-10)
    assert np.allclose(whitened, expected_whitened, atol=1e-9, rtol=1e-9)

    C = observer_from_subspace(H, B)
    geometry = visible_geometry(H, C)
    H_half = _sqrt_spd(H)
    Pi = B @ B.T
    expected_lift = _inv_sqrt_spd(H) @ B
    expected_projector = np.eye(H.shape[0]) - _inv_sqrt_spd(H) @ Pi @ H_half

    assert np.allclose(geometry.phi, np.eye(B.shape[1]), atol=1e-9, rtol=1e-9)
    assert np.allclose(geometry.lift, expected_lift, atol=1e-9, rtol=1e-9)
    assert np.allclose(geometry.projector, expected_projector, atol=1e-9, rtol=1e-9)


def test_closure_scores_exact_invariant_family_and_zero_curvature_edge_case() -> None:
    rng = np.random.default_rng(1102)
    H, U, family, _lambdas = _commuting_family(rng, 7, 4, active_dim=3)
    B_star = U[:, :3]
    scores_star = closure_scores(H, family, B_star)
    scores_zero = closure_scores(H, [np.zeros_like(H)], B_star)
    B_random = orthogonal_matrix(rng, 7)[:, :3]
    scores_random = closure_scores(H, family, B_random)

    assert scores_star.leakage < 1e-10
    assert scores_star.eta < 1e-10
    assert scores_star.visible_score > 0.0
    assert scores_random.leakage > 1e-3
    assert scores_random.eta > 1e-3

    assert scores_zero.leakage == 0.0
    assert scores_zero.visible_score == 0.0
    assert scores_zero.total_curvature == 0.0
    assert scores_zero.eta == 0.0
    assert "numerically zero" in scores_zero.metadata.notes[0]


def test_commuting_exact_solver_matches_exhaustive_search() -> None:
    rng = np.random.default_rng(1103)
    H, U, family, lambdas = _commuting_family(rng, 6, 3, active_dim=5)
    result = closure_adapted_observer(H, family, 2)

    rows = []
    for subset in itertools.combinations(range(6), 2):
        B = U[:, list(subset)]
        scores = closure_scores(H, family, B)
        rows.append((subset, scores.visible_score, scores.leakage, scores.eta))
    exact_rows = [row for row in rows if abs(row[2]) < 1e-10]
    best_subset, best_score, _, _ = max(exact_rows, key=lambda row: row[1])

    assert result.scores.leakage < 1e-10
    assert result.scores.eta < 1e-10
    assert np.isclose(result.scores.visible_score, best_score, atol=1e-10, rtol=1e-9)

    mu = np.sum(lambdas**2, axis=0)
    expected_score = float(np.sum(np.sort(mu)[::-1][:2]))
    assert np.isclose(result.scores.visible_score, expected_score, atol=1e-10, rtol=1e-9)


def test_commuting_solver_handles_repeated_common_eigenspace() -> None:
    H = np.diag([2.0, 3.0, 5.0, 7.0])
    U = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, np.sqrt(0.5), -np.sqrt(0.5)],
            [0.0, 0.0, np.sqrt(0.5), np.sqrt(0.5)],
        ]
    )
    H_half = _sqrt_spd(H)
    D1 = U @ np.diag([4.0, 4.0, 1.5, -0.5]) @ U.T
    D2 = U @ np.diag([3.0, 3.0, 0.0, 2.0]) @ U.T
    family = [H_half @ D1 @ H_half, H_half @ D2 @ H_half]

    result = closure_adapted_observer(H, family, 2)
    projector = result.projector
    whitened_family = [whitened_perturbation(H, Delta) for Delta in family]

    for whitened in whitened_family:
        commutator = whitened @ projector - projector @ whitened
        assert np.linalg.norm(commutator, ord="fro") < 1e-10
    assert result.scores.eta < 1e-10
    assert result.scores.visible_score > 0.0


def test_leakage_channels_recover_q_and_hidden_coupling_rank() -> None:
    H = np.eye(4)
    Delta = np.diag([3.0, 0.0, 0.0, -1.0])
    theta = 0.35
    B = np.array(
        [
            [np.cos(theta), 0.0],
            [0.0, 1.0],
            [np.sin(theta), 0.0],
            [0.0, 0.0],
        ]
    )

    channels = leakage_channels(H, Delta, B)
    scores = closure_scores(H, [Delta], B)

    assert np.allclose(channels.leakage_gram, channels.coupling.T @ channels.coupling, atol=1e-10, rtol=1e-10)
    assert np.isclose(np.sum(channels.singular_values**2), scores.leakage, atol=1e-10, rtol=1e-10)
    assert channels.visible_channel_basis.shape[1] == channels.hidden_channel_basis.shape[1] == 1
    assert channels.metadata.support_rank == 1


def test_compare_observers_detects_same_rank_dominance_and_condition_warning() -> None:
    rng = np.random.default_rng(1104)
    H, U, family, _lambdas = _commuting_family(rng, 6, 3, active_dim=3)
    B_star = U[:, :2]
    B_bad = orthogonal_matrix(rng, 6)[:, :2]
    comparison = compare_observers(H, family, B_star, B_bad)

    assert comparison.left_dominates
    assert not comparison.right_dominates
    assert comparison.leakage_delta < 0.0
    assert comparison.visible_score_delta > 0.0

    H_bad = np.diag([1.0, 1e7])
    B = np.eye(2)[:, :1]
    warning_scores = closure_scores(H_bad, [np.diag([1.0, 0.0])], B)
    assert warning_scores.metadata.condition_number is not None
    assert warning_scores.metadata.condition_number >= 1e7
    assert any("ill-conditioned" in note or "conditioned" in note for note in warning_scores.metadata.notes)


def test_closure_adapted_input_rejection() -> None:
    H = np.eye(3)
    family = [np.diag([1.0, 0.0, -1.0])]
    B_bad = np.array([[1.0, 1.0], [0.0, 0.0], [0.0, 0.0]])

    with pytest.raises(InputValidationError, match="orthonormal"):
        observer_from_subspace(H, B_bad)

    with pytest.raises(InputValidationError, match="at least one perturbation"):
        closure_scores(H, [], np.eye(3)[:, :1])

    with pytest.raises(InputValidationError, match="same shape as H"):
        closure_scores(H, [np.eye(2)], np.eye(3)[:, :1])

    with pytest.raises(InputValidationError, match="pairwise commuting"):
        closure_adapted_observer(
            H,
            [
                np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
                np.diag([0.0, 1.0, -1.0]),
            ],
            1,
        )

    with pytest.raises(InputValidationError, match="currently implemented"):
        closure_adapted_observer(H, family, 1, mode="frontier")

    with pytest.raises(InputValidationError, match="rank must satisfy"):
        closure_adapted_observer(H, family, 0)

    with pytest.raises(InputValidationError, match="same visible rank"):
        compare_observers(H, family, np.eye(3)[:, :1], np.eye(3)[:, :2])
