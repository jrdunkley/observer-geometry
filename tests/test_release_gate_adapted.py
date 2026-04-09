from __future__ import annotations

import itertools

import numpy as np
import pytest

from nomogeo import (
    closure_adapted_observer,
    closure_scores,
    compare_observers,
    dv_bridge,
    leakage_channels,
    local_visible_calculus,
    observer_from_subspace,
    visible_precision,
)
from nomogeo.exceptions import InputValidationError

from .helpers import orthogonal_matrix, random_spd


def _sqrt_spd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T


def test_commuting_solver_exhaustive_agreement_over_seeded_panel() -> None:
    rng = np.random.default_rng(1201)
    for _ in range(12):
        n = int(rng.integers(5, 8))
        m = int(rng.integers(1, min(4, n)))
        task_count = int(rng.integers(2, 5))
        H = random_spd(rng, n)
        H_half = _sqrt_spd(H)
        U = orthogonal_matrix(rng, n)
        lambdas = rng.normal(size=(task_count, n))
        lambdas[:, m + 1 :] = 0.0
        family = [H_half @ U @ np.diag(row) @ U.T @ H_half for row in lambdas]

        result = closure_adapted_observer(H, family, m)

        exhaustive = []
        for subset in itertools.combinations(range(n), m):
            B = U[:, list(subset)]
            scores = closure_scores(H, family, B)
            exhaustive.append((subset, scores.visible_score, scores.leakage))
        exact_rows = [row for row in exhaustive if abs(row[2]) < 1e-10]
        _best_subset, best_score, _ = max(exact_rows, key=lambda row: row[1])

        assert result.scores.leakage < 1e-10
        assert result.scores.eta < 1e-10
        assert np.isclose(result.scores.visible_score, best_score, atol=1e-10, rtol=1e-9)


def test_bridge_adapted_observer_is_fisher_tight_and_beats_seeded_random_observers() -> None:
    rng = np.random.default_rng(1202)
    leakage_random = []
    score_random = []

    for _ in range(8):
        H0 = random_spd(rng, 7)
        J = rng.normal(size=(7, 7))
        J = J - J.T
        bridge = dv_bridge(H0, J)

        adapted = closure_adapted_observer(H0, [bridge.delta_dv], 2)
        local = local_visible_calculus(H0, adapted.C, bridge.delta_dv)
        assert adapted.scores.leakage < 1e-10
        assert adapted.scores.eta < 1e-10
        assert np.linalg.norm(local.Q, ord="fro") < 1e-10

        eps = 0.2
        phi_eps = visible_precision(H0 + (eps * eps) * bridge.delta_dv, adapted.C)
        residual = phi_eps - np.eye(2) - (eps * eps) * local.V
        assert np.linalg.norm(residual, ord="fro") / (eps**4) < 1e-8

        for _ in range(80):
            B = orthogonal_matrix(rng, 7)[:, :2]
            C = observer_from_subspace(H0, B)
            scores = closure_scores(H0, [bridge.delta_dv], B)
            leakage_random.append(scores.eta)
            score_random.append(scores.visible_score)
            assert C.shape == (2, 7)

        assert adapted.scores.visible_score > np.max(score_random[-80:])
        assert adapted.scores.eta < np.min(leakage_random[-80:])

        random_basis = orthogonal_matrix(rng, 7)[:, :2]
        comparison = compare_observers(H0, [bridge.delta_dv], adapted.B, random_basis)
        assert comparison.left_dominates

        channels = leakage_channels(H0, bridge.delta_dv, random_basis)
        assert np.isclose(
            np.sum(channels.singular_values**2),
            closure_scores(H0, [bridge.delta_dv], random_basis).leakage,
            atol=1e-10,
            rtol=1e-9,
        )


def test_total_captured_curvature_identity_and_tower_obstruction() -> None:
    rng = np.random.default_rng(1203)
    for _ in range(40):
        n = 6
        m = 2
        family = []
        for _ in range(3):
            raw = rng.normal(size=(n, n))
            family.append(0.5 * (raw + raw.T))
        B = orthogonal_matrix(rng, n)[:, :m]
        scores = closure_scores(np.eye(n), family, B)
        Pi = B @ B.T
        M = sum(delta @ delta for delta in family)
        assert np.isclose(np.trace(Pi @ M), scores.total_curvature, atol=1e-10, rtol=1e-9)

    D1 = np.array([[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    D2 = np.array([[3.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    rank_one = closure_scores(np.eye(3), [D1, D2], np.array([[1.0], [0.0], [0.0]]))
    assert rank_one.eta == 0.0
    quotient_commutator = D1[1:, 1:] @ D2[1:, 1:] - D2[1:, 1:] @ D1[1:, 1:]
    assert np.linalg.norm(quotient_commutator, ord="fro") > 1.0

    best_eta = float("inf")
    for theta in np.linspace(0.0, 2.0 * np.pi, 2001):
        v = np.array([0.0, np.cos(theta), np.sin(theta)])
        B = np.column_stack([np.array([1.0, 0.0, 0.0]), v])
        scores = closure_scores(np.eye(3), [D1, D2], B)
        best_eta = min(best_eta, scores.eta)
    assert best_eta > 0.05


def test_commuting_solver_rejects_noncommuting_input_loudly() -> None:
    with pytest.raises(InputValidationError, match="pairwise commuting"):
        closure_adapted_observer(
            np.eye(4),
            [
                np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, -1.0]]),
                np.diag([2.0, 1.0, -1.0, 0.5]),
            ],
            2,
        )
