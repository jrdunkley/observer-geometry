from __future__ import annotations

import numpy as np

from nomogeo import canonical_lift, hidden_projector, local_visible_calculus, visible_precision

from .helpers import random_spd, random_surjective


def test_visible_precision_matches_constrained_variational_solution() -> None:
    rng = np.random.default_rng(101)
    H = random_spd(rng, 5)
    C = random_surjective(rng, 3, 5)
    phi = visible_precision(H, C)
    y = rng.normal(size=3)

    kkt = np.block([[H, C.T], [C, np.zeros((3, 3), dtype=float)]])
    rhs = np.concatenate([np.zeros(5, dtype=float), y])
    solution = np.linalg.solve(kkt, rhs)
    x_star = solution[:5]

    assert np.allclose(C @ x_star, y)
    assert np.allclose(x_star.T @ H @ x_star, y.T @ phi @ y)


def test_schur_complement_matches_visible_precision() -> None:
    rng = np.random.default_rng(102)
    H = random_spd(rng, 5)
    C = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
        ]
    )
    phi = visible_precision(H, C)
    H_vv = H[:2, :2]
    H_vh = H[:2, 2:]
    H_hh = H[2:, 2:]
    schur = H_vv - H_vh @ np.linalg.solve(H_hh, H_vh.T)
    assert np.allclose(phi, schur)


def test_tower_law_holds_for_surjective_maps() -> None:
    rng = np.random.default_rng(103)
    H = random_spd(rng, 6)
    C1 = random_surjective(rng, 4, 6)
    C2 = random_surjective(rng, 2, 4)
    composed = C2 @ C1

    direct = visible_precision(H, composed)
    staged = visible_precision(visible_precision(H, C1), C2)
    assert np.allclose(direct, staged)


def test_lift_and_projector_satisfy_identities() -> None:
    rng = np.random.default_rng(104)
    H = random_spd(rng, 5)
    C = random_surjective(rng, 3, 5)
    L = canonical_lift(H, C)
    P = hidden_projector(H, C)

    assert np.allclose(C @ L, np.eye(3), atol=1e-9)
    assert np.allclose(C @ P, np.zeros((3, 5)), atol=1e-9)
    assert np.allclose(P @ P, P, atol=1e-9)
    assert np.allclose(P @ L, np.zeros((5, 3)), atol=1e-9)


def test_local_visible_calculus_matches_finite_differences_and_q_gram() -> None:
    rng = np.random.default_rng(105)
    H = random_spd(rng, 5)
    C = random_surjective(rng, 3, 5)
    A = rng.normal(size=(5, 5))
    Delta = A + A.T
    result = local_visible_calculus(H, C, Delta)

    step = 1e-5
    phi_plus = visible_precision(H + step * Delta, C)
    phi_minus = visible_precision(H - step * Delta, C)
    v_fd = (phi_plus - phi_minus) / (2.0 * step)
    q_fd = -(phi_plus + phi_minus - 2.0 * result.phi) / (2.0 * step * step)

    gram_temp = result.projector.T @ Delta @ result.lift
    q_gram = gram_temp.T @ np.linalg.solve(H, gram_temp)

    assert np.allclose(result.V, v_fd, atol=5e-7, rtol=5e-5)
    assert np.allclose(result.Q, q_fd, atol=5e-6, rtol=5e-4)
    assert np.allclose(result.Q, q_gram, atol=1e-9, rtol=1e-9)
    assert np.min(np.linalg.eigvalsh(result.Q)) >= -1e-8
    assert result.det_split >= -1e-8

