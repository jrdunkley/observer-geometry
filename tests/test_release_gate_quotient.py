from __future__ import annotations

import numpy as np

from nomogeo import gaussian_data_processing_contraction, observer_collapse_descends, visible_precision

from .helpers import random_invertible, random_spd


def test_gaussian_contraction_over_seeded_ensembles_and_equality_cases() -> None:
    rng = np.random.default_rng(801)

    strict_witness = False
    for _ in range(12):
        H1 = random_spd(rng, 4)
        H2 = random_spd(rng, 4)
        C1 = np.eye(4)
        D = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0]])
        result = gaussian_data_processing_contraction(H1, H2, C1, D)

        assert result.forward_kl_coarse <= result.forward_kl_fine + 1e-10
        assert result.reverse_kl_coarse <= result.reverse_kl_fine + 1e-10
        assert result.bhattacharyya_coarse <= result.bhattacharyya_fine + 1e-10
        assert result.hellinger_sq_coarse <= result.hellinger_sq_fine + 1e-10

        if result.forward_kl_coarse < result.forward_kl_fine - 1e-8:
            strict_witness = True

    assert strict_witness

    H1 = np.diag([1.0, 2.0])
    H2 = np.diag([1.5, 2.5])
    C1 = np.eye(2)
    D = random_invertible(rng, 2)
    equality = gaussian_data_processing_contraction(H1, H2, C1, D)

    assert np.isclose(equality.forward_kl_fine, equality.forward_kl_coarse, atol=1e-9, rtol=1e-9)
    assert np.isclose(equality.reverse_kl_fine, equality.reverse_kl_coarse, atol=1e-9, rtol=1e-9)
    assert np.isclose(equality.bhattacharyya_fine, equality.bhattacharyya_coarse, atol=1e-9, rtol=1e-9)
    assert np.isclose(equality.hellinger_sq_fine, equality.hellinger_sq_coarse, atol=1e-9, rtol=1e-9)


def test_collapse_asymmetry_examples_in_both_directions() -> None:
    H1 = np.diag([2.0, 3.0, 5.0])
    b = np.array([[1.0], [2.0]])
    d = np.array([[5.0]])
    A2 = np.diag([2.0, 3.0]) + b @ np.linalg.solve(d, b.T)
    H2 = np.block([[A2, b], [b.T, d]])
    C_rich = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    D = np.array([[1.0, 0.0]])

    assert np.allclose(visible_precision(H1, C_rich), visible_precision(H2, C_rich), atol=1e-9, rtol=1e-9)
    assert observer_collapse_descends(H1, H2, C_rich, D)

    H_a = np.array([[1.0, 0.0], [0.0, 4.0]])
    H_b = np.array([[1.0, 0.0], [0.0, 1.0]])
    coarse = np.array([[1.0, 0.0]])
    rich = (1.0 / np.sqrt(2.0)) * np.array([[1.0, 1.0]])

    assert np.allclose(visible_precision(H_a, coarse), visible_precision(H_b, coarse), atol=1e-9, rtol=1e-9)
    assert not np.allclose(visible_precision(H_a, rich), visible_precision(H_b, rich), atol=1e-9, rtol=1e-9)
