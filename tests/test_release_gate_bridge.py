from __future__ import annotations

import numpy as np

from nomogeo import dv_bridge, hidden_load, local_visible_calculus, visible_precision

from .helpers import orthogonal_matrix, random_skew, random_spd, random_surjective


def test_dv_bridge_trivial_current_sign_symmetry_and_basis_covariance() -> None:
    rng = np.random.default_rng(701)
    H0 = random_spd(rng, 5)
    Jhat = random_skew(rng, 5)

    zero_bridge = dv_bridge(H0, np.zeros_like(H0))
    assert np.allclose(zero_bridge.delta_dv, np.zeros_like(H0), atol=1e-12, rtol=1e-12)
    assert np.allclose(zero_bridge.h_dv, H0, atol=1e-12, rtol=1e-12)

    bridge = dv_bridge(H0, Jhat)
    bridge_neg = dv_bridge(H0, -Jhat)
    assert np.allclose(bridge.delta_dv, bridge.gram_factor @ bridge.gram_factor.T, atol=1e-9, rtol=1e-9)
    assert np.allclose(bridge.delta_dv, bridge_neg.delta_dv, atol=1e-9, rtol=1e-9)

    G = orthogonal_matrix(rng, H0.shape[0])
    bridge_cov = dv_bridge(G.T @ H0 @ G, G.T @ Jhat @ G)
    assert np.allclose(bridge_cov.delta_dv, G.T @ bridge.delta_dv @ G, atol=1e-9, rtol=1e-9)
    assert np.allclose(bridge_cov.h_dv, G.T @ bridge.h_dv @ G, atol=1e-9, rtol=1e-9)


def test_dv_bridge_small_parameter_visible_and_hidden_rates() -> None:
    rng = np.random.default_rng(702)
    H0 = random_spd(rng, 4)
    J1 = random_skew(rng, 4)
    C = random_surjective(rng, 2, 4)
    delta2 = dv_bridge(H0, J1).delta_dv
    local = local_visible_calculus(H0, C, delta2)

    eps_values = np.array([1e-1, 5e-2, 2e-2, 1e-2], dtype=float)
    visible_residuals = []
    hidden_normalised_residuals = []
    v_eigs, v_vecs = np.linalg.eigh(local.V)
    keep = v_eigs > 1e-10
    v_basis = v_vecs[:, keep]
    v_s = v_basis.T @ local.V @ v_basis
    q_s = v_basis.T @ local.Q @ v_basis
    evals, evecs = np.linalg.eigh(v_s)
    inv_sqrt_v = (evecs * (1.0 / np.sqrt(evals))) @ evecs.T
    A = inv_sqrt_v @ q_s @ inv_sqrt_v

    for eps in eps_values:
        bridge_eps = dv_bridge(H0, eps * J1)
        H_eps = bridge_eps.h_dv
        X_eps = visible_precision(H_eps, C) - visible_precision(H0, C)
        visible_residuals.append(np.linalg.norm(X_eps - (eps * eps) * local.V + (eps**4) * local.Q, ord="fro"))
        lambda_eps = hidden_load((eps * eps) * local.V, X_eps, support_mode="auto").reduced_lambda
        hidden_normalised_residuals.append(np.linalg.norm(lambda_eps / (eps * eps) - A, ord="fro"))

    visible_slope, _ = np.polyfit(np.log(eps_values), np.log(visible_residuals), 1)

    assert visible_slope > 3.0
    assert max(hidden_normalised_residuals) < 2e-6
