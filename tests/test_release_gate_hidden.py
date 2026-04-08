from __future__ import annotations

import numpy as np

from nomogeo import (
    clock,
    hidden_load,
    hidden_contraction,
    load_from_hidden_contraction,
    transport_hidden_load,
    visible_from_hidden_load,
)

from .helpers import orthogonal_matrix, random_psd, random_psd_on_subspace


def _psd_violation(matrix: np.ndarray) -> float:
    if matrix.size == 0:
        return 0.0
    return float(max(0.0, -np.min(np.linalg.eigvalsh(matrix))))


def test_hidden_representation_consistency_and_support_gauge_covariance() -> None:
    rng = np.random.default_rng(601)
    T, _basis = random_psd_on_subspace(rng, 7, 4)
    lambda_reduced = random_psd(rng, 4)
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    load = hidden_load(T, X)

    lambda_ambient = load.support_basis @ load.reduced_lambda @ load.support_basis.T
    x_from_ambient = visible_from_hidden_load(T, lambda_ambient, lambda_representation="ambient")
    x_from_reduced = visible_from_hidden_load(T, load.reduced_lambda, lambda_representation="reduced")
    assert np.allclose(x_from_ambient, x_from_reduced, atol=1e-9, rtol=1e-9)

    U = orthogonal_matrix(rng, load.support_basis.shape[1])
    B_prime = load.support_basis @ U
    T_support = load.support_basis.T @ T @ load.support_basis
    X_support = load.support_basis.T @ X @ load.support_basis
    T_rot = U.T @ T_support @ U
    X_rot = U.T @ X_support @ U
    load_rot = hidden_load(T_rot, X_rot, support_mode="ambient")

    assert np.allclose(load_rot.reduced_lambda, U.T @ load.reduced_lambda @ U, atol=1e-9, rtol=1e-9)
    x_rot_ambient = B_prime @ X_rot @ B_prime.T
    assert np.allclose(x_rot_ambient, X, atol=1e-9, rtol=1e-9)


def test_hidden_rank_clock_loewner_and_near_ceiling_regimes_across_support_ranks() -> None:
    rng = np.random.default_rng(602)
    eps_values = np.array([1e-1, 5e-2, 2e-2, 1e-2], dtype=float)

    for n in range(2, 7):
        for rank in range(0, n + 1):
            T, _basis = random_psd_on_subspace(rng, n, rank)
            lambda_reduced = random_psd(rng, rank) if rank > 0 else np.zeros((0, 0), dtype=float)
            X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
            load = hidden_load(T, X)

            gap_rank = np.linalg.matrix_rank(T - X, tol=max(load.metadata.rank_tol, 1e-12))
            assert load.rank == gap_rank

            if rank == 0:
                assert load.clock == 0.0
                continue

            support_basis = load.support_basis
            t_s = support_basis.T @ T @ support_basis
            x_s = support_basis.T @ X @ support_basis
            lhs = np.linalg.slogdet(t_s)[1] - np.linalg.slogdet(x_s)[1]
            rhs = np.linalg.slogdet(np.eye(rank) + load.reduced_lambda)[1]
            assert np.isclose(lhs, rhs, atol=1e-9, rtol=1e-9)

            extra = random_psd(rng, rank)
            x_1 = visible_from_hidden_load(T, load.reduced_lambda + extra, lambda_representation="reduced")
            x_2 = visible_from_hidden_load(T, load.reduced_lambda, lambda_representation="reduced")
            lambda_1 = hidden_load(T, x_1).reduced_lambda
            lambda_2 = hidden_load(T, x_2).reduced_lambda
            assert _psd_violation(x_2 - x_1) < 1e-10
            assert _psd_violation(lambda_1 - lambda_2) < 1e-10

            A = random_psd(rng, rank)
            residuals = []
            for eps in eps_values:
                x_eps = visible_from_hidden_load(T, eps * A, lambda_representation="reduced")
                lambda_eps = hidden_load(T, x_eps).reduced_lambda
                residuals.append(np.linalg.norm(lambda_eps - eps * A, ord="fro"))
            assert max(residuals) < 1e-8


def test_hidden_transport_identity_downstairs_semigroup_long_chain_and_scalar_commuting_law() -> None:
    rng = np.random.default_rng(603)
    zero = np.zeros((3, 3), dtype=float)
    lambda_a = random_psd(rng, 3)
    lambda_b = random_psd(rng, 3)
    lambda_c = random_psd(rng, 3)

    assert np.allclose(transport_hidden_load(lambda_a, zero), lambda_a, atol=1e-9, rtol=1e-9)
    assert np.allclose(transport_hidden_load(zero, lambda_a), lambda_a, atol=1e-9, rtol=1e-9)

    two_step = transport_hidden_load(lambda_a, lambda_b)
    from_factors = load_from_hidden_contraction(hidden_contraction(lambda_b) @ hidden_contraction(lambda_a))
    assert np.allclose(two_step, from_factors, atol=1e-9, rtol=1e-9)

    k_a = hidden_contraction(lambda_a)
    k_b = hidden_contraction(lambda_b)
    k_c = hidden_contraction(lambda_c)
    left_factor = k_c @ (k_b @ k_a)
    right_factor = (k_c @ k_b) @ k_a
    left_load = load_from_hidden_contraction(left_factor)
    right_load = load_from_hidden_contraction(right_factor)
    assert np.allclose(left_factor, right_factor, atol=1e-12, rtol=1e-12)
    assert np.allclose(left_load, right_load, atol=1e-9, rtol=1e-9)

    chain = [0.05 * random_psd(rng, 3) for _ in range(20)]
    factor_total = np.eye(3, dtype=float)
    for load in chain:
        factor_total = hidden_contraction(load) @ factor_total
    total_from_chain = load_from_hidden_contraction(factor_total)
    pi_total = factor_total.T @ factor_total
    total_from_pi = hidden_load(np.eye(3), pi_total, support_mode="ambient").reduced_lambda

    assert np.allclose(total_from_chain, total_from_pi, atol=1e-8, rtol=1e-8)
    assert np.isclose(clock(total_from_chain), sum(clock(load) for load in chain), atol=1e-8, rtol=1e-8)

    diag_a = np.diag([0.1, 0.2, 0.3])
    diag_b = np.diag([0.4, 0.5, 0.6])
    diag_total = transport_hidden_load(diag_a, diag_b)
    expected = np.diag(np.diag(diag_a) + np.diag(diag_b) + np.diag(diag_a) * np.diag(diag_b))
    assert np.allclose(diag_total, expected, atol=1e-9, rtol=1e-9)


def test_local_hidden_birth_first_order_fit_and_rate() -> None:
    rng = np.random.default_rng(604)
    V = np.array([[2.0, 0.3], [0.3, 1.6]])
    Q = random_psd(rng, 2)
    evals, evecs = np.linalg.eigh(V)
    inv_sqrt_v = (evecs * (1.0 / np.sqrt(evals))) @ evecs.T
    A = inv_sqrt_v @ Q @ inv_sqrt_v

    ts = np.array([1e-1, 5e-2, 2e-2, 1e-2, 5e-3], dtype=float)
    residuals = []
    for t in ts:
        X_t = t * V - (t * t) * Q
        T_t = t * V
        lambda_t = hidden_load(T_t, X_t, support_mode="ambient").reduced_lambda
        residuals.append(np.linalg.norm(lambda_t - t * A, ord="fro"))
    slope, _ = np.polyfit(np.log(ts), np.log(residuals), 1)

    assert slope > 1.7
    assert residuals[-1] < residuals[0]
