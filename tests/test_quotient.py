from __future__ import annotations

import numpy as np

from nomogeo import (
    gaussian_data_processing_contraction,
    hidden_load,
    observer_collapse_descends,
    visible_precision,
    visible_from_hidden_load,
)

from .helpers import random_psd_on_subspace, random_spd, random_surjective


def test_hidden_load_parametrisation_relative_to_reference_ceiling() -> None:
    rng = np.random.default_rng(401)
    H_ref = random_spd(rng, 5)
    H = random_spd(rng, 5)
    C = random_surjective(rng, 3, 5)

    T = visible_precision(H_ref, C)
    lambda_reduced = np.diag([0.2, 0.4, 0.8])
    X = visible_from_hidden_load(
        T,
        lambda_reduced,
        support_mode="ambient",
        lambda_representation="reduced",
    )
    load = hidden_load(T, X, support_mode="ambient")

    assert np.allclose(load.reduced_lambda, lambda_reduced, atol=1e-9, rtol=1e-9)


def test_collapse_asymmetry_numerically() -> None:
    H_a = np.array([[1.0, 0.0], [0.0, 4.0]])
    H_b = np.array([[1.0, 0.0], [0.0, 1.0]])
    coarse = np.array([[1.0, 0.0]])
    rich = (1.0 / np.sqrt(2.0)) * np.array([[1.0, 1.0]])

    assert np.allclose(visible_precision(H_a, coarse), visible_precision(H_b, coarse))
    assert not np.allclose(visible_precision(H_a, rich), visible_precision(H_b, rich))


def test_equality_under_richer_observer_descends_to_coarsening() -> None:
    H = np.diag([1.0, 2.0, 3.0])
    C1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    D = np.array([[1.0, 1.0]])
    assert observer_collapse_descends(H, H, C1, D)


def test_gaussian_data_processing_contracts_under_coarsening() -> None:
    rng = np.random.default_rng(402)
    H1 = random_spd(rng, 4)
    H2 = random_spd(rng, 4)
    C1 = np.eye(4)
    D = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 0.0],
        ]
    )
    result = gaussian_data_processing_contraction(H1, H2, C1, D)

    assert result.forward_kl_coarse <= result.forward_kl_fine + 1e-10
    assert result.reverse_kl_coarse <= result.reverse_kl_fine + 1e-10
    assert result.hellinger_sq_coarse <= result.hellinger_sq_fine + 1e-10
    assert result.bhattacharyya_coarse <= result.bhattacharyya_fine + 1e-10
