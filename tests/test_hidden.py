from __future__ import annotations

import numpy as np
import pytest

from nomogeo.exceptions import InputValidationError
from nomogeo import (
    canonical_hidden_realisation,
    clock,
    hidden_load,
    minimal_hidden_realisation,
    transport_hidden_load,
    visible_from_hidden_load,
)

from .helpers import orthogonal_matrix, random_psd_on_subspace, schur_complement


def test_hidden_load_forward_inverse_are_mutual_inverses() -> None:
    rng = np.random.default_rng(201)
    T, _basis = random_psd_on_subspace(rng, 5, 3)
    lambda_reduced = np.diag([0.2, 0.5, 1.0])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    result = hidden_load(T, X)

    rebuilt = visible_from_hidden_load(T, result.lambda_, lambda_representation="ambient")
    assert np.allclose(rebuilt, X, atol=1e-9, rtol=1e-9)
    assert np.allclose(result.reduced_lambda, lambda_reduced, atol=1e-9, rtol=1e-9)


def test_canonical_hidden_realisation_has_correct_schur_complement() -> None:
    rng = np.random.default_rng(202)
    T, _basis = random_psd_on_subspace(rng, 6, 4)
    lambda_reduced = np.diag([0.1, 0.3, 0.6, 0.9])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    canonical = canonical_hidden_realisation(T, X)
    load = hidden_load(T, X)

    schur = schur_complement(canonical.matrix, visible_dim=4)
    x_support = load.support_basis.T @ X @ load.support_basis
    assert np.allclose(schur, x_support, atol=1e-9, rtol=1e-9)
    assert np.min(np.linalg.eigvalsh(canonical.matrix)) > 0.0


def test_minimal_hidden_realisation_has_correct_rank_and_gauge_invariance() -> None:
    rng = np.random.default_rng(203)
    T, _basis = random_psd_on_subspace(rng, 5, 3)
    lambda_reduced = np.diag([0.0, 0.4, 1.2])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    minimal = minimal_hidden_realisation(T, X)
    load = hidden_load(T, X)

    q = orthogonal_matrix(rng, minimal.factor.shape[1])
    rotated = minimal.factor @ q

    assert minimal.rank == 2
    assert np.allclose(minimal.factor @ minimal.factor.T, load.reduced_lambda, atol=1e-9)
    assert np.allclose(rotated @ rotated.T, load.reduced_lambda, atol=1e-9)


def test_transport_law_and_clock_additivity_hold() -> None:
    lambda_a = np.diag([0.1, 0.5, 0.9])
    lambda_b = np.diag([0.3, 0.2, 0.4])
    total = transport_hidden_load(lambda_a, lambda_b)

    pi = np.linalg.inv(np.eye(3) + lambda_a)
    xi = np.linalg.inv(np.eye(3) + lambda_b)
    sequential = np.linalg.multi_dot([np.linalg.cholesky(pi), xi, np.linalg.cholesky(pi).T])
    expected_total = np.linalg.inv(sequential) - np.eye(3)

    assert np.allclose(total, expected_total, atol=1e-9, rtol=1e-9)
    assert np.isclose(clock(total), clock(lambda_a) + clock(lambda_b), atol=1e-9, rtol=1e-9)


def test_local_hidden_birth_matches_first_order_generator() -> None:
    V = np.array([[2.0, 0.4], [0.4, 1.5]])
    Q = np.array([[0.5, 0.1], [0.1, 0.3]])
    t = 1e-4
    X_t = t * V - t * t * Q
    T_t = t * V
    result = hidden_load(T_t, X_t, support_mode="ambient")

    v_eigs, v_vecs = np.linalg.eigh(V)
    inv_sqrt_v = (v_vecs * (1.0 / np.sqrt(v_eigs))) @ v_vecs.T
    generator = inv_sqrt_v @ Q @ inv_sqrt_v

    assert np.allclose(result.reduced_lambda / t, generator, atol=5e-4, rtol=5e-4)


def test_zero_support_inverse_reconstruction_returns_zero_matrix() -> None:
    T = np.zeros((3, 3), dtype=float)
    X = visible_from_hidden_load(T, np.zeros((0, 0), dtype=float), lambda_representation="reduced")
    assert np.allclose(X, np.zeros((3, 3), dtype=float))


def test_full_rank_nondiagonal_ceiling_requires_explicit_representation() -> None:
    T = np.array([[2.0, 0.7], [0.7, 1.5]])
    evals, evecs = np.linalg.eigh(T)
    lambda_reduced = np.diag([0.2, 0.6])
    lambda_ambient = evecs @ lambda_reduced @ evecs.T

    with pytest.raises(InputValidationError):
        visible_from_hidden_load(T, lambda_reduced)

    x_from_reduced = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    x_from_ambient = visible_from_hidden_load(T, lambda_ambient, lambda_representation="ambient")
    assert np.allclose(x_from_reduced, x_from_ambient, atol=1e-9, rtol=1e-9)


def test_clock_rejects_indefinite_hidden_load() -> None:
    with pytest.raises(InputValidationError):
        clock(np.diag([-0.1, 0.2]))


def test_transport_rejects_indefinite_hidden_loads() -> None:
    with pytest.raises(InputValidationError):
        transport_hidden_load(np.diag([-0.1, 0.2]), np.diag([0.1, 0.3]))
    with pytest.raises(InputValidationError):
        transport_hidden_load(np.diag([0.1, 0.2]), np.diag([-0.05, 0.3]))


def test_hidden_load_rank_matches_rank_gap_identity() -> None:
    rng = np.random.default_rng(204)
    T, _basis = random_psd_on_subspace(rng, 6, 4)
    lambda_reduced = np.diag([0.0, 0.2, 0.7, 1.4])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    load = hidden_load(T, X)
    gap_rank = np.linalg.matrix_rank(T - X, tol=load.metadata.rank_tol)
    assert load.rank == gap_rank


def test_hidden_load_clock_matches_pdet_identity() -> None:
    rng = np.random.default_rng(205)
    T, _basis = random_psd_on_subspace(rng, 5, 3)
    lambda_reduced = np.diag([0.2, 0.4, 0.9])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    load = hidden_load(T, X)

    support_basis = load.support_basis
    t_s = support_basis.T @ T @ support_basis
    x_s = support_basis.T @ X @ support_basis
    lhs = np.linalg.slogdet(t_s)[1] - np.linalg.slogdet(x_s)[1]
    rhs = np.linalg.slogdet(np.eye(load.reduced_lambda.shape[0]) + load.reduced_lambda)[1]
    assert np.isclose(lhs, rhs, atol=1e-9, rtol=1e-9)
