from __future__ import annotations

import numpy as np
import pytest

from nomogeo.exceptions import InputValidationError
from nomogeo.field import (
    classify_support_event_from_jet,
    comparison_envelope_bounds,
    hidden_load_from_pi,
    interval_hessian_at_exact_family,
    kernel_schur_jet_from_coefficients,
    lambda_rhs,
    local_coupled_birth,
    pi_from_hidden_load,
    pi_rhs,
    restart_hidden_load_birth,
    restart_hidden_load_death,
    sampled_interval_closure_check,
    sampled_interval_leakage,
    sampled_interval_stationarity,
    semisimple_event_block,
    support_stratum_transport,
)


def _sym(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def _rot(theta: float) -> np.ndarray:
    return np.array(
        [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ],
        dtype=float,
    )


def test_support_stratum_transport_and_envelopes() -> None:
    Lambda = np.diag([0.3, 0.8])
    A = np.array([[0.4, 0.1], [0.1, 0.2]])
    result = support_stratum_transport(Lambda, A)
    Pi = pi_from_hidden_load(Lambda)
    Lambda_roundtrip = hidden_load_from_pi(Pi)

    assert np.allclose(result.pi, Pi, atol=1e-10, rtol=1e-10)
    assert np.allclose(Lambda_roundtrip, Lambda, atol=1e-10, rtol=1e-10)
    assert np.allclose(result.lambda_rhs, lambda_rhs(Lambda, A), atol=1e-10, rtol=1e-10)
    assert np.allclose(result.pi_rhs, pi_rhs(Pi, A), atol=1e-10, rtol=1e-10)
    assert np.isclose(result.clock_rate, np.trace(A), atol=1e-12)
    assert result.generator_psd

    lower, upper = comparison_envelope_bounds(Pi, 0.2, 0.6)
    assert np.all(np.linalg.eigvalsh(upper - lower) >= -1e-12)
    assert np.all(np.linalg.eigvalsh(Pi - upper) >= -1e-12)

    indefinite = support_stratum_transport(Lambda, np.diag([0.2, -0.1]))
    assert not indefinite.generator_psd
    assert "not positive semidefinite" in indefinite.metadata.notes[0]
    with pytest.raises(InputValidationError, match="positive semidefinite"):
        support_stratum_transport(Lambda, np.diag([0.2, -0.1]), require_psd_generator=True)


def test_restart_maps_handle_rotated_nested_bases() -> None:
    lambda_before = np.array([[0.6, 0.1], [0.1, 0.3]])
    old_basis = np.eye(3)[:, :2]
    theta = 0.37
    new_basis = np.column_stack([old_basis @ _rot(theta), np.eye(3)[:, 2]])

    birth = restart_hidden_load_birth(lambda_before, old_basis, new_basis)
    expected_ambient = old_basis @ lambda_before @ old_basis.T
    assert birth.event_kind == "birth"
    assert birth.lambda_after.shape == (3, 3)
    assert np.allclose(new_basis @ birth.lambda_after @ new_basis.T, expected_ambient, atol=1e-10, rtol=1e-10)

    survivor_basis = old_basis @ _rot(theta)[:, :1]
    death = restart_hidden_load_death(lambda_before, old_basis, survivor_basis)
    expected_death = survivor_basis.T @ old_basis @ lambda_before @ old_basis.T @ survivor_basis
    assert death.event_kind == "death"
    assert death.lambda_after.shape == (1, 1)
    assert np.allclose(death.lambda_after, expected_death, atol=1e-10, rtol=1e-10)

    with pytest.raises(InputValidationError, match="span"):
        restart_hidden_load_birth(lambda_before, old_basis, np.eye(3)[:, [0, 2]])


def test_kernel_jet_first_second_third_and_semisimple_blocks() -> None:
    V0 = np.diag([0.0, 0.0, 3.0])
    V1 = np.diag([2.0, -1.0, 0.0])
    first = kernel_schur_jet_from_coefficients([V0, V1])
    assert first.order == 1
    assert first.birth_count_forward == 1
    assert first.death_count_forward == 1
    assert first.event_kind == "mixed"
    assert classify_support_event_from_jet(first) == "mixed"

    V0_second = np.diag([0.0, 1.0])
    V1_second = np.array([[0.0, 1.0], [1.0, 0.0]])
    V2_second = np.zeros((2, 2))
    second = kernel_schur_jet_from_coefficients([V0_second, V1_second, V2_second])
    assert second.order == 2
    assert np.allclose(second.leading_effective, np.array([[-1.0]]), atol=1e-10, rtol=1e-10)
    assert second.event_kind == "degenerate"

    V2_third = np.array([[1.0, 0.0], [0.0, 0.0]])
    V3_third = np.array([[2.0, 0.0], [0.0, 0.0]])
    third = kernel_schur_jet_from_coefficients([V0_second, V1_second, V2_third, V3_third])
    assert third.order == 3
    assert np.allclose(third.leading_effective, np.array([[2.0]]), atol=1e-10, rtol=1e-10)

    touch = kernel_schur_jet_from_coefficients([V0_second, np.zeros((2, 2)), np.diag([1.5, 0.0])])
    assert touch.order == 2
    assert touch.event_kind == "touch"
    right_block = semisimple_event_block(touch, side_sign=1)
    left_block = semisimple_event_block(touch, side_sign=-1)
    assert right_block.birth_like
    assert np.isclose(right_block.pole_coefficient, -1.0)
    assert left_block.death_like
    assert np.isclose(left_block.clock_log_coefficient, 1.0)
    assert np.isclose(left_block.desingularisation_power, 1.0)

    degenerate = kernel_schur_jet_from_coefficients([V0_second, np.zeros((2, 2)), np.zeros((2, 2))])
    assert degenerate.event_kind == "degenerate"
    assert degenerate.order is None


def test_local_coupled_birth_extractor_is_hidden_basis_invariant() -> None:
    H0 = np.array([[3.0, 0.4, -0.2], [0.4, 2.2, 0.3], [-0.2, 0.3, 1.7]], dtype=float)
    H1 = np.array([[0.2, -0.1, 0.05], [-0.1, 0.15, 0.07], [0.05, 0.07, -0.08]], dtype=float)
    H2 = np.array([[0.06, 0.03, -0.02], [0.03, -0.04, 0.01], [-0.02, 0.01, 0.05]], dtype=float)

    def H(t: float) -> np.ndarray:
        return H0 + t * H1 + 0.5 * t * t * H2

    def H_dot(t: float) -> np.ndarray:
        return H1 + t * H2

    def C(t: float) -> np.ndarray:
        return np.array([[1.0, t, 0.2 * t * t]], dtype=float)

    def C_dot(t: float) -> np.ndarray:
        return np.array([[0.0, 1.0, 0.4 * t]], dtype=float)

    def Z(t: float) -> np.ndarray:
        return np.array([[-t, -0.2 * t * t], [1.0, 0.0], [0.0, 1.0]], dtype=float)

    def T(t: float) -> np.ndarray:
        return np.array([[np.exp(0.3 * t), 0.2 * t], [-0.1 * t, 1.0 + 0.15 * t]], dtype=float)

    t0 = 0.25
    base = local_coupled_birth(H(t0), H_dot(t0), H2, C(t0), C_dot(t0), Z(t0))
    changed = local_coupled_birth(H(t0), H_dot(t0), H2, C(t0), C_dot(t0), Z(t0) @ T(t0))

    assert np.allclose(base.Q, changed.Q, atol=1e-10, rtol=1e-10)
    assert np.allclose(base.observer_tensor, changed.observer_tensor, atol=1e-10, rtol=1e-10)
    assert np.allclose(base.W, changed.W, atol=1e-10, rtol=1e-10)
    assert np.allclose(base.a_cpl, changed.a_cpl, atol=1e-10, rtol=1e-10)

    h = 1e-6

    def V_of(t: float) -> np.ndarray:
        result = local_coupled_birth(H(t), H_dot(t), H2, C(t), C_dot(t), Z(t))
        return result.V

    phi = base.phi
    lift = base.lift
    alpha = -C_dot(t0) @ lift
    direct_w = (V_of(t0 + h) - V_of(t0 - h)) / (2.0 * h) - alpha.T @ (lift.T @ H_dot(t0) @ lift) - (lift.T @ H_dot(t0) @ lift) @ alpha
    assert np.allclose(base.W, direct_w, atol=5e-8, rtol=5e-8)
    assert phi.shape == (1, 1)


def test_local_coupled_birth_scalar_normal_form_recovers_pole() -> None:
    phi_star = 2.0
    a = 1.6
    r = 3.0
    b = 0.9
    t = 0.04
    H = np.array([[phi_star + 0.5 * a * t * t + b * t * b * t * r, b * t * r], [b * t * r, r]])
    H_dot = np.array([[a * t + 2.0 * b * b * t * r, b * r], [b * r, 0.0]])
    H_ddot = np.array([[a + 2.0 * b * b * r, 0.0], [0.0, 0.0]])
    C = np.array([[1.0, 0.0]])
    C_dot = np.zeros_like(C)
    Z = np.array([[0.0], [1.0]])

    result = local_coupled_birth(H, H_dot, H_ddot, C, C_dot, Z)
    assert np.allclose(result.a_cpl, np.array([[-0.5 / t]]), atol=1e-10, rtol=1e-10)


def test_sampled_interval_leakage_stationarity_and_hessian() -> None:
    family = [
        np.diag([4.0, 2.0, 0.4]),
        np.diag([3.5, 2.5, 0.6]),
        np.diag([3.0, 2.8, 0.8]),
    ]
    basis = np.eye(3)[:, :2]
    result = sampled_interval_leakage(family, basis)
    assert result.sampled_exact_closure
    assert result.leakage == 0.0
    assert np.allclose(sampled_interval_stationarity(family, basis), 0.0, atol=1e-12, rtol=1e-12)
    assert sampled_interval_closure_check(family, basis)

    broken = [member.copy() for member in family]
    broken[1][0, 2] = broken[1][2, 0] = 0.2
    broken_result = sampled_interval_leakage(broken, basis)
    assert broken_result.leakage > 0.0
    assert not broken_result.sampled_exact_closure

    hessian = interval_hessian_at_exact_family(family, basis)
    assert hessian.locally_rigid
    assert hessian.spectral_gap is not None
    assert hessian.rigidity_lower_bound is not None
    assert np.min(np.linalg.eigvalsh(hessian.hessian_operator)) >= hessian.rigidity_lower_bound - 1e-10


def test_field_input_rejection() -> None:
    with pytest.raises(InputValidationError, match="Pi <= I"):
        hidden_load_from_pi(2.0 * np.eye(2))
    with pytest.raises(InputValidationError, match="same shape"):
        lambda_rhs(np.eye(2), np.eye(3))
    with pytest.raises(InputValidationError, match="same shape as C"):
        local_coupled_birth(np.eye(2), np.eye(2), np.eye(2), np.array([[1.0, 0.0]]), np.eye(2))
    with pytest.raises(InputValidationError, match="coefficients"):
        kernel_schur_jet_from_coefficients([np.eye(2)])
    with pytest.raises(InputValidationError, match="weights"):
        sampled_interval_leakage([np.eye(2)], np.eye(2)[:, :1], weights=[1.0, 2.0])
