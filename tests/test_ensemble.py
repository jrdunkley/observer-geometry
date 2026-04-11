import numpy as np
import scipy.linalg as la
import pytest

from nomogeo import (
    ceiling_mediated_local_quadratic_ensemble,
    coordinate_local_quadratic_ensemble,
    hidden_load,
    intrinsic_local_quadratic_ensemble,
    visible_precision,
)
from nomogeo.exceptions import InputValidationError


def _sample_hessians() -> np.ndarray:
    rng = np.random.default_rng(2001)
    n = 6
    samples = []
    for strength in [0.0, 0.08, 0.18, 0.31]:
        coupling = np.zeros((n, n))
        raw = rng.normal(size=(3, 3))
        coupling[:3, 3:] = raw
        coupling[3:, :3] = raw.T
        h = np.eye(n) + strength * coupling + (0.6 + strength**2) * np.eye(n)
        samples.append(0.5 * (h + h.T))
    return np.stack(samples)


def test_coordinate_local_quadratic_ensemble_matches_samplewise_primitives() -> None:
    hessians = _sample_hessians()
    visible_dim = 3

    result = coordinate_local_quadratic_ensemble(hessians, visible_dim)

    assert result.sample_count == hessians.shape[0]
    assert result.phis.shape == (4, 3, 3)
    assert result.lambdas.shape == (4, 3, 3)
    assert result.hidden_ranks.shape == (4,)

    C = np.eye(visible_dim, hessians.shape[1])
    for idx, H in enumerate(hessians):
        phi = visible_precision(H, C)
        load = hidden_load(H[:visible_dim, :visible_dim], phi, support_mode="ambient")
        assert np.allclose(result.phis[idx], phi, atol=1e-12, rtol=1e-10)
        assert np.allclose(result.lambdas[idx], load.lambda_, atol=1e-12, rtol=1e-10)
        assert result.clocks[idx] == pytest.approx(load.clock)

    assert result.mean_clock == pytest.approx(float(np.mean(result.clocks)))
    assert result.std_clock == pytest.approx(float(np.std(result.clocks)))
    assert result.min_clock == pytest.approx(float(np.min(result.clocks)))
    assert result.max_clock == pytest.approx(float(np.max(result.clocks)))
    assert "not full-law cumulants" in result.metadata.notes[-1]


def test_intrinsic_local_quadratic_ensemble_is_latent_basis_covariant() -> None:
    hessians = _sample_hessians()
    C = np.eye(3, 6)
    rng = np.random.default_rng(2002)
    S = rng.normal(size=(6, 6))
    while abs(np.linalg.det(S)) < 0.1:
        S = rng.normal(size=(6, 6))
    transformed = np.stack([S.T @ H @ S for H in hessians])

    original = intrinsic_local_quadratic_ensemble(hessians, C)
    moved = intrinsic_local_quadratic_ensemble(transformed, C @ S)

    assert np.allclose(original.phis, moved.phis, atol=1e-10, rtol=1e-9)
    assert np.allclose(original.logdet_phis, moved.logdet_phis, atol=1e-10, rtol=1e-9)


def test_ceiling_mediated_ensemble_requires_explicit_ceiling_and_matches_coordinate_helper() -> None:
    hessians = _sample_hessians()
    C = np.eye(3, 6)
    ceilings = hessians[:, :3, :3]

    mediated = ceiling_mediated_local_quadratic_ensemble(hessians, C, ceilings)
    coordinate = coordinate_local_quadratic_ensemble(hessians, 3)

    assert np.allclose(mediated.phis, coordinate.phis, atol=1e-12, rtol=1e-10)
    assert np.allclose(mediated.lambdas, coordinate.lambdas, atol=1e-12, rtol=1e-10)
    assert np.allclose(mediated.clocks, coordinate.clocks, atol=1e-12, rtol=1e-10)
    assert "explicitly supplied ceilings" in mediated.metadata.notes[0]


def test_ceiling_mediated_hidden_load_clock_is_visible_basis_invariant() -> None:
    hessians = _sample_hessians()
    C = np.eye(3, 6)
    ceilings = hessians[:, :3, :3]
    U = np.array([[1.0, 0.2, 0.0], [0.1, 1.2, -0.1], [0.0, 0.3, 0.9]])
    U_inv = np.linalg.inv(U)
    transformed_ceilings = np.stack([U_inv.T @ T @ U_inv for T in ceilings])

    original = ceiling_mediated_local_quadratic_ensemble(hessians, C, ceilings)
    moved = ceiling_mediated_local_quadratic_ensemble(hessians, U @ C, transformed_ceilings)

    assert np.allclose(original.clocks, moved.clocks, atol=1e-10, rtol=1e-9)
    for idx in range(original.sample_count):
        eig0 = la.eigvalsh(original.lambdas[idx])
        eig1 = la.eigvalsh(moved.lambdas[idx])
        assert np.allclose(eig0, eig1, atol=1e-10, rtol=1e-9)


def test_coordinate_local_quadratic_ensemble_detects_nontrivial_distribution() -> None:
    hessians = _sample_hessians()
    result = coordinate_local_quadratic_ensemble(hessians, 3)

    assert result.max_clock > result.min_clock
    assert result.std_clock > 0.0

    mean_hessian = np.mean(hessians, axis=0)
    mean_phi = visible_precision(mean_hessian, np.eye(3, 6))
    mean_load = hidden_load(mean_hessian[:3, :3], mean_phi, support_mode="ambient")
    assert abs(mean_load.clock - result.mean_clock) > 1e-5


def test_coordinate_local_quadratic_ensemble_rejects_bad_inputs() -> None:
    with pytest.raises(InputValidationError):
        coordinate_local_quadratic_ensemble(np.eye(3), 1)
    with pytest.raises(InputValidationError):
        coordinate_local_quadratic_ensemble(np.zeros((0, 3, 3)), 1)
    with pytest.raises(InputValidationError):
        coordinate_local_quadratic_ensemble(np.stack([np.eye(3)]), 0)
    with pytest.raises(InputValidationError):
        coordinate_local_quadratic_ensemble(np.stack([np.eye(3)]), 3)
    bad = np.stack([np.diag([1.0, 0.0, 2.0])])
    with pytest.raises(InputValidationError):
        coordinate_local_quadratic_ensemble(bad, 1)


def test_ceiling_mediated_ensemble_rejects_bad_ceilings() -> None:
    hessians = _sample_hessians()
    C = np.eye(3, 6)
    with pytest.raises(InputValidationError):
        ceiling_mediated_local_quadratic_ensemble(hessians, C, np.eye(2))
    too_small = np.stack([0.1 * np.eye(3) for _ in range(hessians.shape[0])])
    with pytest.raises(InputValidationError):
        ceiling_mediated_local_quadratic_ensemble(hessians, C, too_small)
