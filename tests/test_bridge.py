from __future__ import annotations

import numpy as np

from nomogeo import dv_bridge, hidden_load, local_visible_calculus, visible_precision

from .helpers import random_skew, random_spd, random_surjective


def test_dv_bridge_matches_gram_factor_and_is_positive() -> None:
    rng = np.random.default_rng(301)
    H0 = random_spd(rng, 5)
    Jhat = random_skew(rng, 5)
    result = dv_bridge(H0, Jhat)

    assert np.allclose(result.delta_dv, result.gram_factor @ result.gram_factor.T, atol=1e-9)
    assert np.min(np.linalg.eigvalsh(result.delta_dv)) >= -1e-9
    assert np.min(np.linalg.eigvalsh(result.h_dv)) > 0.0


def test_dv_bridge_has_quadratic_visible_onset_and_quartic_hidden_birth() -> None:
    rng = np.random.default_rng(302)
    H0 = random_spd(rng, 4)
    J1 = random_skew(rng, 4)
    C = random_surjective(rng, 2, 4)

    delta2 = dv_bridge(H0, J1).delta_dv
    local = local_visible_calculus(H0, C, delta2)
    eps = 5e-3
    H_eps = H0 + (eps * eps) * delta2
    X_eps = visible_precision(H_eps, C) - visible_precision(H0, C)

    assert np.allclose(X_eps / (eps * eps), local.V, atol=1e-5, rtol=1e-3)
    assert np.allclose((X_eps - (eps * eps) * local.V) / (eps**4), -local.Q, atol=3e-3, rtol=2e-2)

    load = hidden_load((eps * eps) * local.V, X_eps, support_mode="auto")
    v_eigs, v_vecs = np.linalg.eigh(local.V)
    keep = v_eigs > 1e-10
    v_basis = v_vecs[:, keep]
    v_s = v_basis.T @ local.V @ v_basis
    q_s = v_basis.T @ local.Q @ v_basis
    evals, evecs = np.linalg.eigh(v_s)
    inv_sqrt_v = (evecs * (1.0 / np.sqrt(evals))) @ evecs.T
    generator = inv_sqrt_v @ q_s @ inv_sqrt_v

    assert np.allclose(load.reduced_lambda / (eps * eps), generator, atol=2e-2, rtol=2e-2)

