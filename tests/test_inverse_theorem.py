from __future__ import annotations

import numpy as np
import pytest

from nomogeo import hidden_load, inverse_visible_class, visible_from_hidden_load
from nomogeo.exceptions import InputValidationError

from .helpers import orthogonal_matrix, random_psd, random_psd_on_subspace


def _relative_error(a: np.ndarray, b: np.ndarray) -> float:
    denom = max(1.0, float(np.linalg.norm(b, ord="fro")))
    return float(np.linalg.norm(a - b, ord="fro") / denom)


def _psd_violation(matrix: np.ndarray) -> float:
    if matrix.size == 0:
        return 0.0
    return float(max(0.0, -np.min(np.linalg.eigvalsh(matrix))))


def test_inverse_visible_class_is_exact_alias_of_visible_from_hidden_load() -> None:
    rng = np.random.default_rng(701)
    T, _basis = random_psd_on_subspace(rng, 6, 4)
    lambda_reduced = random_psd(rng, 4)

    x_alias = inverse_visible_class(T, lambda_reduced, lambda_representation="reduced")
    x_direct = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    assert np.allclose(x_alias, x_direct, atol=1e-10, rtol=1e-10)


def test_inverse_theorem_roundtrips_hold_across_support_ranks_and_rotated_supports() -> None:
    rng = np.random.default_rng(702)

    for n in range(1, 7):
        for rank in range(0, n + 1):
            T, _basis = random_psd_on_subspace(rng, n, rank)
            if rank == 0:
                X = inverse_visible_class(T, np.zeros((0, 0), dtype=float), lambda_representation="reduced")
                recovered = hidden_load(T, X)
                assert np.allclose(X, np.zeros_like(T))
                assert recovered.reduced_lambda.shape == (0, 0)
                continue

            lambda_reduced = 0.2 * random_psd(rng, rank)
            X = inverse_visible_class(T, lambda_reduced, lambda_representation="reduced")
            recovered = hidden_load(T, X)

            assert _relative_error(recovered.reduced_lambda, lambda_reduced) < 1e-10
            assert _relative_error(inverse_visible_class(T, recovered.reduced_lambda, lambda_representation="reduced"), X) < 1e-10

            U = orthogonal_matrix(rng, rank)
            support_basis = recovered.support_basis
            T_s = support_basis.T @ T @ support_basis
            X_s = support_basis.T @ X @ support_basis
            T_rot = U.T @ T_s @ U
            X_rot = U.T @ X_s @ U
            recovered_rot = hidden_load(T_rot, X_rot, support_mode="ambient")
            assert _relative_error(recovered_rot.reduced_lambda, U.T @ lambda_reduced @ U) < 1e-10


def test_inverse_theorem_order_reversal_holds_beneath_fixed_ceiling() -> None:
    rng = np.random.default_rng(703)

    for n in range(2, 7):
        for rank in range(1, n + 1):
            T, _basis = random_psd_on_subspace(rng, n, rank)
            lambda_1 = 0.1 * random_psd(rng, rank)
            lambda_2 = lambda_1 + 0.15 * random_psd(rng, rank)
            x_1 = inverse_visible_class(T, lambda_1, lambda_representation="reduced")
            x_2 = inverse_visible_class(T, lambda_2, lambda_representation="reduced")
            assert _psd_violation(x_1 - x_2) < 1e-10


def test_inverse_theorem_rank_and_clock_identities_hold_on_active_support() -> None:
    rng = np.random.default_rng(704)

    for n in range(2, 7):
        for rank in range(1, n + 1):
            T, _basis = random_psd_on_subspace(rng, n, rank)
            lambda_reduced = np.diag(np.linspace(0.0, 1.2, rank))
            X = inverse_visible_class(T, lambda_reduced, lambda_representation="reduced")
            load = hidden_load(T, X)

            gap_rank = np.linalg.matrix_rank(T - X, tol=max(load.metadata.rank_tol, 1e-12))
            assert load.rank == gap_rank

            support_basis = load.support_basis
            t_s = support_basis.T @ T @ support_basis
            x_s = support_basis.T @ X @ support_basis
            lhs = np.linalg.slogdet(t_s)[1] - np.linalg.slogdet(x_s)[1]
            rhs = np.linalg.slogdet(np.eye(rank) + lambda_reduced)[1]
            assert np.isclose(lhs, rhs, atol=1e-10, rtol=1e-10)


def test_inverse_theorem_full_rank_ambiguity_requires_explicit_representation() -> None:
    T = np.array([[2.0, 0.4], [0.4, 1.7]], dtype=float)
    lambda_reduced = np.diag([0.1, 0.6])
    evals, evecs = np.linalg.eigh(T)
    lambda_ambient = evecs @ lambda_reduced @ evecs.T

    with pytest.raises(InputValidationError):
        inverse_visible_class(T, lambda_reduced)

    x_reduced = inverse_visible_class(T, lambda_reduced, lambda_representation="reduced")
    x_ambient = inverse_visible_class(T, lambda_ambient, lambda_representation="ambient")
    assert np.allclose(x_reduced, x_ambient, atol=1e-10, rtol=1e-10)


def test_inverse_theorem_handles_near_ceiling_and_large_load_regimes() -> None:
    rng = np.random.default_rng(705)
    T, _basis = random_psd_on_subspace(rng, 6, 4)

    near_ceiling = 1e-12 * np.diag([1.0, 2.0, 3.0, 4.0])
    x_near = inverse_visible_class(T, near_ceiling, lambda_representation="reduced")
    recovered_near = hidden_load(T, x_near).reduced_lambda
    assert np.linalg.norm(x_near - T, ord="fro") < 1e-10
    assert np.linalg.norm(recovered_near - near_ceiling, ord="fro") < 1e-10

    large_load = np.diag([1e1, 1e2, 1e3, 1e4])
    x_large = inverse_visible_class(T, large_load, lambda_representation="reduced")
    recovered_large = hidden_load(T, x_large).reduced_lambda
    assert _relative_error(recovered_large, large_load) < 1e-9
