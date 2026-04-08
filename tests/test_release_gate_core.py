from __future__ import annotations

import numpy as np

from nomogeo import local_visible_calculus, visible_geometry, visible_precision

from .helpers import null_space_basis, random_invertible, random_spd, random_surjective


def _kappa(H: np.ndarray, C: np.ndarray) -> float:
    phi = visible_precision(H, C)
    sign, value = np.linalg.slogdet(phi)
    assert sign > 0
    return float(-value)


def test_visible_geometry_full_projector_and_energy_split() -> None:
    rng = np.random.default_rng(501)
    H = random_spd(rng, 6)
    C = random_surjective(rng, 3, 6)
    geometry = visible_geometry(H, C)
    phi = geometry.phi
    L = geometry.lift
    P = geometry.projector
    N = null_space_basis(C)

    assert np.allclose(C @ L, np.eye(C.shape[0]), atol=1e-9, rtol=1e-9)
    assert np.allclose(C @ P, np.zeros_like(C @ P), atol=1e-9, rtol=1e-9)
    assert np.allclose(P @ P, P, atol=1e-9, rtol=1e-9)
    assert np.allclose(P @ L, np.zeros_like(L), atol=1e-9, rtol=1e-9)
    assert np.allclose(L @ C + P, np.eye(H.shape[0]), atol=1e-9, rtol=1e-9)
    assert np.allclose(P.T @ H, H @ P, atol=1e-9, rtol=1e-9)

    y = rng.normal(size=C.shape[0])
    z = N @ rng.normal(size=N.shape[1])
    x = L @ y + z
    split = float(y.T @ phi @ y + z.T @ H @ z)

    assert np.allclose(x.T @ H @ x, split, atol=1e-9, rtol=1e-9)
    assert np.allclose(L.T @ H @ N, np.zeros((C.shape[0], N.shape[1])), atol=1e-9, rtol=1e-9)


def test_local_visible_calculus_pure_visible_and_hidden_annihilating_families() -> None:
    rng = np.random.default_rng(502)
    H = random_spd(rng, 6)
    C = random_surjective(rng, 3, 6)
    geometry = visible_geometry(H, C)

    A = rng.normal(size=(3, 3))
    A = 0.5 * (A + A.T)
    delta_visible = C.T @ A @ C
    visible_result = local_visible_calculus(H, C, delta_visible)

    null_lift = null_space_basis(geometry.lift.T)
    B = rng.normal(size=(null_lift.shape[1], null_lift.shape[1]))
    B = 0.5 * (B + B.T)
    delta_hidden = null_lift @ B @ null_lift.T
    hidden_result = local_visible_calculus(H, C, delta_hidden)

    assert np.allclose(visible_result.V, A, atol=1e-9, rtol=1e-9)
    assert np.allclose(visible_result.Q, np.zeros_like(A), atol=1e-9, rtol=1e-9)
    assert np.allclose(delta_hidden @ geometry.lift, np.zeros_like(geometry.lift), atol=1e-9, rtol=1e-9)
    assert np.allclose(hidden_result.V, np.zeros_like(hidden_result.V), atol=1e-9, rtol=1e-9)
    assert np.allclose(hidden_result.Q, np.zeros_like(hidden_result.Q), atol=1e-9, rtol=1e-9)


def test_determinant_curvature_theorem_and_rate() -> None:
    rng = np.random.default_rng(503)
    H = random_spd(rng, 5)
    C = random_surjective(rng, 2, 5)
    A = rng.normal(size=(5, 5))
    Delta = 0.08 * (A + A.T) / np.linalg.norm(A + A.T, ord=2)
    local = local_visible_calculus(H, C, Delta)

    step = 1e-5
    curvature_fd = (_kappa(H + step * Delta, C) - 2.0 * _kappa(H, C) + _kappa(H - step * Delta, C)) / (step * step)
    assert np.isclose(curvature_fd, local.det_split, atol=5e-5, rtol=5e-4)

    ts = np.array([1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3], dtype=float)
    residuals = []
    for t in ts:
        phi_t = visible_precision(H + t * Delta, C)
        residual = phi_t - local.phi - t * local.V + (t * t) * local.Q
        residuals.append(np.linalg.norm(residual, ord="fro"))
    slope, _ = np.polyfit(np.log(ts), np.log(residuals), 1)
    assert slope > 2.7


def test_multistage_tower_parenthesisations_and_inverse_relation() -> None:
    rng = np.random.default_rng(504)
    H = random_spd(rng, 8)
    C1 = random_surjective(rng, 6, 8)
    C2 = random_surjective(rng, 5, 6)
    C3 = random_surjective(rng, 3, 5)
    C4 = random_surjective(rng, 2, 3)

    direct = visible_precision(H, C4 @ C3 @ C2 @ C1)
    staged_left = visible_precision(visible_precision(visible_precision(visible_precision(H, C1), C2), C3), C4)
    staged_split = visible_precision(visible_precision(H, C2 @ C1), C4 @ C3)
    staged_right = visible_precision(H, C4 @ (C3 @ (C2 @ C1)))

    gram = (C4 @ C3 @ C2 @ C1) @ np.linalg.solve(H, (C4 @ C3 @ C2 @ C1).T)

    assert np.allclose(direct, staged_left, atol=1e-9, rtol=1e-9)
    assert np.allclose(direct, staged_split, atol=1e-9, rtol=1e-9)
    assert np.allclose(direct, staged_right, atol=1e-9, rtol=1e-9)
    assert np.allclose(gram, np.linalg.inv(direct), atol=1e-9, rtol=1e-9)


def test_latent_and_visible_basis_covariance() -> None:
    rng = np.random.default_rng(505)
    H = random_spd(rng, 6)
    C = random_surjective(rng, 3, 6)
    A = rng.normal(size=(6, 6))
    Delta = 0.5 * (A + A.T)

    geometry = visible_geometry(H, C)
    local = local_visible_calculus(H, C, Delta)

    G = random_invertible(rng, 6)
    q_g, _ = np.linalg.qr(G)
    H_latent = q_g.T @ H @ q_g
    C_latent = C @ q_g
    geometry_latent = visible_geometry(H_latent, C_latent)

    assert np.allclose(geometry_latent.phi, geometry.phi, atol=1e-9, rtol=1e-9)
    assert np.allclose(geometry_latent.lift, q_g.T @ geometry.lift, atol=1e-9, rtol=1e-9)
    assert np.allclose(geometry_latent.projector, q_g.T @ geometry.projector @ q_g, atol=1e-9, rtol=1e-9)

    U = random_invertible(rng, 3)
    C_visible = U @ C
    local_visible = local_visible_calculus(H, C_visible, Delta)
    U_inv = np.linalg.inv(U)
    transform = U_inv.T

    assert np.allclose(local_visible.phi, transform @ local.phi @ U_inv, atol=1e-9, rtol=1e-9)
    assert np.allclose(local_visible.V, transform @ local.V @ U_inv, atol=1e-9, rtol=1e-9)
    assert np.allclose(local_visible.Q, transform @ local.Q @ U_inv, atol=1e-9, rtol=1e-9)
