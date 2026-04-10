from __future__ import annotations

import numpy as np

import nomogeo
from nomogeo import (
    kernel_schur_jet_from_coefficients,
    lambda_rhs,
    local_coupled_birth,
    restart_hidden_load_birth,
    sampled_interval_leakage,
    support_stratum_transport,
)


def test_field_public_surface_and_version_gate() -> None:
    assert nomogeo.__version__ == "0.30.0"

    Lambda = np.diag([0.2, 0.5])
    A = np.diag([0.3, 0.7])
    transport = support_stratum_transport(Lambda, A)
    assert transport.generator_psd
    assert np.allclose(transport.lambda_rhs, lambda_rhs(Lambda, A), atol=1e-10, rtol=1e-10)

    old_basis = np.eye(3)[:, :2]
    new_basis = np.eye(3)
    restart = restart_hidden_load_birth(Lambda, old_basis, new_basis)
    assert restart.lambda_after.shape == (3, 3)
    assert np.allclose(restart.lambda_after[:2, :2], Lambda, atol=1e-10, rtol=1e-10)
    assert np.allclose(restart.lambda_after[2:, :], 0.0, atol=1e-10, rtol=1e-10)

    V0 = np.diag([0.0, 1.0])
    V1 = np.array([[0.0, 1.0], [1.0, 0.0]])
    V2 = np.array([[1.0, 0.0], [0.0, 0.0]])
    V3 = np.array([[2.0, 0.0], [0.0, 0.0]])
    jet = kernel_schur_jet_from_coefficients([V0, V1, V2, V3])
    assert jet.order == 3

    H = np.array([[2.0 + 0.5 * 0.2 * 0.2, 0.2], [0.2, 1.0]])
    H_dot = np.array([[0.5, 1.0], [1.0, 0.0]])
    H_ddot = np.array([[0.0, 0.0], [0.0, 0.0]])
    local = local_coupled_birth(H, H_dot, H_ddot, np.array([[1.0, 0.0]]), np.zeros((1, 2)), np.array([[0.0], [1.0]]))
    assert local.a_cpl.shape == (1, 1)

    family = [np.diag([2.0, 1.0]), np.diag([3.0, 0.5])]
    interval = sampled_interval_leakage(family, np.eye(2)[:, :1])
    assert interval.sampled_exact_closure
