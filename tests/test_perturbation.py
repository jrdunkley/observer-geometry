import numpy as np
import pytest

from nomogeo import (
    rank_k_covariance_perturbation,
    rank_one_covariance_perturbation,
    residual_margin_ordering,
    simple_spectrum_closure_certificate,
)
from nomogeo.exceptions import InputValidationError


def _hidden_gap(sigma: np.ndarray, visible_dim: int) -> np.ndarray:
    h = np.linalg.inv(sigma)
    phi = np.linalg.inv(sigma[:visible_dim, :visible_dim])
    return 0.5 * ((h[:visible_dim, :visible_dim] - phi) + (h[:visible_dim, :visible_dim] - phi).T)


def test_rank_one_covariance_perturbation_white_sector_is_one_channel() -> None:
    rng = np.random.default_rng(1001)
    n, m = 8, 3
    eps = 0.29
    signal = rng.normal(size=n)

    result = rank_one_covariance_perturbation(np.eye(n), signal, m, eps)
    direct = _hidden_gap(np.eye(n) + eps * np.outer(signal, signal), m) - _hidden_gap(np.eye(n), m)

    assert result.update_rank == 1
    assert result.one_channel
    assert np.allclose(result.hidden_gap_increment, direct, atol=1e-12, rtol=1e-10)
    assert result.formula_residual <= 1e-12

    a = signal[:m]
    expected_coeff = eps**2 * float(signal[m:] @ signal[m:]) / (
        (1.0 + eps * float(a @ a)) * (1.0 + eps * float(signal @ signal))
    )
    assert np.allclose(result.hidden_gap_increment, expected_coeff * np.outer(a, a), atol=1e-12, rtol=1e-10)


def test_rank_one_covariance_perturbation_coloured_sector_can_be_two_rank() -> None:
    rng = np.random.default_rng(1002)
    n, m = 9, 4
    eps = 0.37

    found = None
    for _ in range(200):
        A = rng.normal(size=(n, n))
        sigma0 = A.T @ A / n + 0.8 * np.eye(n)
        signal = rng.normal(size=n)
        result = rank_one_covariance_perturbation(sigma0, signal, m, eps)
        if result.update_rank == 2:
            found = result
            break

    assert found is not None
    assert not found.one_channel
    assert found.formula_residual <= 1e-11
    assert found.direction_alignment < 1.0 - 1e-4


def test_rank_one_covariance_perturbation_block_aligned_coloured_sector_is_one_channel() -> None:
    rng = np.random.default_rng(1003)
    n, m = 7, 3
    eps = 0.41
    visible_block = np.diag([1.3, 1.7, 2.1])
    hidden_block = np.diag([0.8, 1.1, 1.5, 1.9])
    sigma0 = np.block(
        [
            [visible_block, np.zeros((m, n - m))],
            [np.zeros((n - m, m)), hidden_block],
        ]
    )
    signal = np.zeros(n)
    signal[1] = 0.7
    signal[m:] = rng.normal(size=n - m)

    result = rank_one_covariance_perturbation(sigma0, signal, m, eps)

    assert result.update_rank == 1
    assert result.one_channel
    assert result.direction_alignment == pytest.approx(1.0)
    assert result.formula_residual <= 1e-12


def test_rank_k_covariance_perturbation_matches_woodbury_formula_and_rank_bound() -> None:
    rng = np.random.default_rng(1005)
    n, m, k = 10, 4, 3
    raw = rng.normal(size=(n, n))
    sigma0 = raw.T @ raw / n + np.eye(n)
    factor = rng.normal(size=(n, k)) / 3.0

    result = rank_k_covariance_perturbation(sigma0, factor, m)
    direct = _hidden_gap(sigma0 + factor @ factor.T, m) - _hidden_gap(sigma0, m)

    assert np.allclose(result.hidden_gap_increment, direct, atol=1e-12, rtol=1e-10)
    assert np.allclose(result.hidden_gap_increment, result.formula_increment, atol=1e-12, rtol=1e-10)
    assert result.formula_residual <= 1e-12
    assert result.update_rank <= result.rank_bound == 2 * k
    assert "not a generic non-Gaussian" in result.metadata.notes[-1]


def test_rank_one_covariance_perturbation_is_rank_k_special_case() -> None:
    rng = np.random.default_rng(1006)
    n, m = 7, 3
    eps = 0.23
    raw = rng.normal(size=(n, n))
    sigma0 = raw.T @ raw / n + 1.2 * np.eye(n)
    signal = rng.normal(size=n)

    rank_one = rank_one_covariance_perturbation(sigma0, signal, m, eps)
    rank_k = rank_k_covariance_perturbation(sigma0, np.sqrt(eps) * signal[:, None], m)

    assert np.allclose(rank_one.hidden_gap_increment, rank_k.hidden_gap_increment)
    assert np.allclose(rank_one.formula_increment, rank_k.formula_increment)
    assert rank_one.update_rank == rank_k.update_rank


def test_rank_one_covariance_perturbation_rejects_bad_inputs() -> None:
    with pytest.raises(InputValidationError):
        rank_one_covariance_perturbation(np.eye(3), np.ones(2), 1, 0.1)
    with pytest.raises(InputValidationError):
        rank_one_covariance_perturbation(np.eye(3), np.ones(3), 0, 0.1)
    with pytest.raises(InputValidationError):
        rank_one_covariance_perturbation(np.eye(3), np.ones(3), 3, 0.1)
    with pytest.raises(InputValidationError):
        rank_one_covariance_perturbation(np.eye(3), np.ones(3), 1, -0.1)
    with pytest.raises(InputValidationError):
        rank_k_covariance_perturbation(np.eye(3), np.ones((2, 1)), 1)
    with pytest.raises(InputValidationError):
        rank_k_covariance_perturbation(np.eye(3), np.ones((3, 0)), 1)


def test_residual_margin_ordering_certifies_and_refuses_thresholds() -> None:
    safe = residual_margin_ordering(0.08, 0.03)
    assert safe.robust
    assert not safe.adversarial_reversal_possible
    assert safe.margin == pytest.approx(0.02)

    exact_threshold = residual_margin_ordering(0.08, 0.04)
    assert not exact_threshold.robust
    assert exact_threshold.adversarial_reversal_possible
    assert exact_threshold.margin == pytest.approx(0.0)

    unsafe = residual_margin_ordering(0.08, 0.05)
    assert not unsafe.robust
    assert unsafe.adversarial_reversal_possible
    assert unsafe.worst_case_gap == pytest.approx(-0.02)

    with pytest.raises(InputValidationError):
        residual_margin_ordering(-0.1, 0.01)
    with pytest.raises(InputValidationError):
        residual_margin_ordering(0.1, -0.01)


def test_simple_spectrum_closure_certificate_finds_exact_coordinate_subspace() -> None:
    D1 = np.diag([1.0, 2.0, 3.0, 4.0])
    D2 = np.array(
        [
            [1.0, 0.2, 0.0, 0.0],
            [0.2, 1.4, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.3],
            [0.0, 0.0, 0.3, 2.5],
        ]
    )

    result = simple_spectrum_closure_certificate([D1, D2], rank=2)

    assert result.exact_common_subspace_exists
    assert result.obstruction_certified is False
    assert result.min_cross_block_norm <= 1e-12
    assert result.best_indices in ((0, 1), (2, 3))


def test_simple_spectrum_closure_certificate_certifies_generic_obstruction() -> None:
    rng = np.random.default_rng(1004)
    D1 = np.diag([1.0, 2.0, 3.0, 4.0, 5.0])
    A = rng.normal(size=(5, 5))
    D2 = 0.5 * (A + A.T)

    result = simple_spectrum_closure_certificate([D1, D2], rank=2)

    assert not result.exact_common_subspace_exists
    assert result.obstruction_certified
    assert result.min_cross_block_norm > 1e-6
    assert result.checked_subset_count == 10


def test_simple_spectrum_closure_certificate_rejects_degenerate_anchor() -> None:
    D1 = np.diag([1.0, 1.0, 2.0])
    D2 = np.diag([3.0, 4.0, 5.0])

    with pytest.raises(InputValidationError):
        simple_spectrum_closure_certificate([D1, D2], rank=1)
