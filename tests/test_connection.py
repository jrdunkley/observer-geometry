from __future__ import annotations

import numpy as np
import pytest

from nomogeo import (
    connection_current,
    fixed_observer_coordinates,
    forcing_from_current,
    observer_transition,
    reconstruct_precision_from_fixed_observer_coordinates,
    visible_precision,
)
from nomogeo.exceptions import InputValidationError

from .helpers import random_spd, random_surjective


def _random_symmetric(rng: np.random.Generator, n: int, scale: float = 1.0) -> np.ndarray:
    raw = rng.normal(size=(n, n))
    symmetric = 0.5 * (raw + raw.T)
    return scale * symmetric / max(np.linalg.norm(symmetric, ord=2), 1e-15)


def _delta_from_chart_variation(
    chart: object,
    dphi: np.ndarray,
    dr: np.ndarray,
    dk: np.ndarray,
) -> np.ndarray:
    phi = chart.phi
    hidden = chart.hidden_block
    coupling = chart.coupling
    basis = chart.adapted_basis
    visible_dim = phi.shape[0]
    hidden_dim = hidden.shape[0]
    triangular = np.block(
        [
            [np.eye(visible_dim, dtype=float), coupling],
            [np.zeros((hidden_dim, visible_dim), dtype=float), np.eye(hidden_dim, dtype=float)],
        ]
    )
    triangular_dot = np.block(
        [
            [np.zeros((visible_dim, visible_dim), dtype=float), dk],
            [np.zeros((hidden_dim, visible_dim), dtype=float), np.zeros((hidden_dim, hidden_dim), dtype=float)],
        ]
    )
    block_diag = np.block(
        [
            [phi, np.zeros((visible_dim, hidden_dim), dtype=float)],
            [np.zeros((hidden_dim, visible_dim), dtype=float), hidden],
        ]
    )
    block_diag_dot = np.block(
        [
            [dphi, np.zeros((visible_dim, hidden_dim), dtype=float)],
            [np.zeros((hidden_dim, visible_dim), dtype=float), dr],
        ]
    )
    delta_tilde = triangular_dot @ block_diag @ triangular.T + triangular @ block_diag_dot @ triangular.T + triangular @ block_diag @ triangular_dot.T
    basis_inv = np.linalg.inv(basis)
    return 0.5 * (basis_inv.T @ delta_tilde @ basis_inv + basis_inv.T @ delta_tilde.T @ basis_inv)


def test_fixed_observer_coordinates_roundtrip_and_phi_recovery() -> None:
    rng = np.random.default_rng(2101)
    H = random_spd(rng, 6)
    C = random_surjective(rng, 2, 6)

    chart = fixed_observer_coordinates(H, C)
    H_recovered = reconstruct_precision_from_fixed_observer_coordinates(
        chart.phi,
        chart.hidden_block,
        chart.coupling,
        chart.adapted_basis,
    )

    assert np.allclose(chart.phi, visible_precision(H, C), atol=1e-10, rtol=1e-10)
    assert np.allclose(H_recovered, H, atol=1e-10, rtol=1e-10)
    assert chart.hidden_block.shape == (4, 4)
    assert chart.coupling.shape == (2, 4)
    assert chart.metadata.method == "fixed-observer-chart"


def test_observer_transition_matches_direct_right_chart() -> None:
    rng = np.random.default_rng(2102)
    H = random_spd(rng, 5)
    C_left = random_surjective(rng, 2, 5)
    C_right = random_surjective(rng, 2, 5)

    transition = observer_transition(H, C_left, C_right)

    assert np.allclose(transition.right_phi_from_left, transition.right.phi, atol=1e-10, rtol=1e-10)
    assert np.allclose(transition.right_hidden_from_left, transition.right.hidden_block, atol=1e-10, rtol=1e-10)
    assert np.allclose(transition.right_coupling_from_left, transition.right.coupling, atol=1e-10, rtol=1e-10)
    assert transition.residual < 1e-10
    assert transition.metadata.method == "observer-transition-law"


def test_connection_current_recovers_forcing_and_fixed_k_zero_edge_case() -> None:
    rng = np.random.default_rng(2103)
    H = random_spd(rng, 6)
    C = random_surjective(rng, 2, 6)
    chart = fixed_observer_coordinates(H, C)

    dphi = _random_symmetric(rng, 2, scale=0.15)
    dr = _random_symmetric(rng, 4, scale=0.08)
    dk = rng.normal(size=(2, 4))
    Delta = _delta_from_chart_variation(chart, dphi, dr, dk)

    current_result = connection_current(H, C, Delta)
    expected_forcing = forcing_from_current(chart.phi, chart.hidden_block, current_result.current)

    assert np.allclose(current_result.forcing, expected_forcing, atol=1e-10, rtol=1e-10)
    assert np.allclose(current_result.forcing, current_result.q, atol=1e-10, rtol=1e-10)
    assert np.linalg.norm(current_result.current, ord="fro") > 1e-6

    Delta_fixed_k = _delta_from_chart_variation(
        chart,
        _random_symmetric(rng, 2, scale=0.12),
        _random_symmetric(rng, 4, scale=0.06),
        np.zeros((2, 4), dtype=float),
    )
    fixed_k_result = connection_current(H, C, Delta_fixed_k)

    assert np.allclose(fixed_k_result.coupling_velocity, 0.0, atol=1e-10, rtol=1e-10)
    assert np.allclose(fixed_k_result.current, 0.0, atol=1e-10, rtol=1e-10)
    assert np.allclose(fixed_k_result.forcing, 0.0, atol=1e-10, rtol=1e-10)
    assert np.allclose(fixed_k_result.q, 0.0, atol=1e-10, rtol=1e-10)


def test_full_visible_observer_zero_hidden_sector() -> None:
    H = np.diag([2.0, 3.0, 5.0])
    C = np.eye(3)
    Delta = np.diag([1.0, -0.5, 0.25])

    chart = fixed_observer_coordinates(H, C)
    reconstructed = reconstruct_precision_from_fixed_observer_coordinates(
        chart.phi,
        chart.hidden_block,
        chart.coupling,
        chart.adapted_basis,
    )
    current_result = connection_current(H, C, Delta)
    transition = observer_transition(H, C, C)

    assert chart.hidden_block.shape == (0, 0)
    assert chart.coupling.shape == (3, 0)
    assert np.allclose(chart.phi, H, atol=1e-10, rtol=1e-10)
    assert np.allclose(reconstructed, H, atol=1e-10, rtol=1e-10)
    assert current_result.current.shape == (3, 0)
    assert np.allclose(current_result.current, 0.0, atol=1e-10, rtol=1e-10)
    assert np.allclose(current_result.forcing, 0.0, atol=1e-10, rtol=1e-10)
    assert np.allclose(current_result.q, 0.0, atol=1e-10, rtol=1e-10)
    assert transition.residual < 1e-12


def test_fixed_observer_coordinates_emit_conditioning_note() -> None:
    H = np.diag([1.0, 2.0, 1e7, 2e7])
    C = np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]])

    chart = fixed_observer_coordinates(H, C)

    assert chart.metadata.condition_number is not None
    assert chart.metadata.condition_number >= 1e7
    assert any("condition" in note for note in chart.metadata.notes)


def test_connection_input_rejection() -> None:
    H = np.eye(4)
    C = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    chart = fixed_observer_coordinates(H, C)

    with pytest.raises(InputValidationError, match="same visible rank"):
        observer_transition(H, C, np.array([[1.0, 0.0, 0.0, 0.0]]))

    with pytest.raises(InputValidationError, match="incompatible shape"):
        reconstruct_precision_from_fixed_observer_coordinates(
            chart.phi,
            chart.hidden_block,
            np.zeros((1, 1), dtype=float),
            chart.adapted_basis,
        )

    with pytest.raises(InputValidationError, match="invertible"):
        reconstruct_precision_from_fixed_observer_coordinates(
            chart.phi,
            chart.hidden_block,
            chart.coupling,
            np.zeros_like(chart.adapted_basis),
        )

    with pytest.raises(InputValidationError, match="same shape as H"):
        connection_current(H, C, np.eye(3))

    with pytest.raises(InputValidationError, match="incompatible shape"):
        forcing_from_current(chart.phi, chart.hidden_block, np.zeros((3, 3), dtype=float))
