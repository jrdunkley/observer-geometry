from __future__ import annotations

import numpy as np
import scipy.linalg as la

from .core import local_visible_calculus, visible_precision
from .exceptions import InputValidationError
from .types import (
    ConnectionCurrentResult,
    FixedObserverCoordinatesResult,
    LinearAlgebraMetadata,
    ObserverTransitionResult,
)
from .validation import (
    Tolerances,
    resolve_tolerances,
    symmetrize,
    to_float_array,
    validate_spd_matrix,
    validate_surjective_map,
    validate_symmetric_matrix,
)


def fixed_observer_coordinates(
    H: np.ndarray,
    C: np.ndarray,
    tolerances: Tolerances | None = None,
) -> FixedObserverCoordinatesResult:
    """Return the exact fixed-observer chart (Phi, R, K) for one observer."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    c = validate_surjective_map(C, h.shape[0], tol)
    adapted_basis = _fixed_adapted_basis(c)
    phi, hidden_block, coupling = _coordinates_from_adapted_basis(h, c.shape[0], adapted_basis)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="fixed-observer-chart",
        ambient_dim=h.shape[0],
        support_rank=hidden_block.shape[0],
        visible_dim=c.shape[0],
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=_condition_notes(h),
    )
    return FixedObserverCoordinatesResult(
        phi=phi,
        hidden_block=hidden_block,
        coupling=coupling,
        adapted_basis=adapted_basis,
        metadata=metadata,
    )


def reconstruct_precision_from_fixed_observer_coordinates(
    phi: np.ndarray,
    hidden_block: np.ndarray,
    coupling: np.ndarray,
    adapted_basis: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Reconstruct H from a fixed-observer chart and its adapted basis."""
    tol = resolve_tolerances(tolerances)
    phi_matrix = validate_spd_matrix("phi", phi, tol)
    hidden_candidate = to_float_array("hidden_block", hidden_block)
    if hidden_candidate.shape == (0, 0):
        hidden = hidden_candidate
    else:
        hidden = validate_spd_matrix("hidden_block", hidden_candidate, tol)
    basis = _validate_adapted_basis(adapted_basis, phi_matrix.shape[0] + hidden.shape[0], tol)
    coupling_matrix = to_float_array("coupling", coupling)
    if coupling_matrix.shape != (phi_matrix.shape[0], hidden.shape[0]):
        raise InputValidationError("coupling has incompatible shape")

    hidden_dim = hidden.shape[0]
    visible_dim = phi_matrix.shape[0]
    triangular = np.block(
        [
            [np.eye(visible_dim, dtype=float), coupling_matrix],
            [np.zeros((hidden_dim, visible_dim), dtype=float), np.eye(hidden_dim, dtype=float)],
        ]
    )
    block_diag = np.block(
        [
            [phi_matrix, np.zeros((visible_dim, hidden_dim), dtype=float)],
            [np.zeros((hidden_dim, visible_dim), dtype=float), hidden],
        ]
    )
    h_tilde = triangular @ block_diag @ triangular.T
    basis_inv = np.linalg.inv(basis)
    h = symmetrize(basis_inv.T @ h_tilde @ basis_inv)
    return validate_spd_matrix("reconstructed_H", h, tol)


def observer_transition(
    H: np.ndarray,
    C_left: np.ndarray,
    C_right: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ObserverTransitionResult:
    """Return the exact observer-to-observer transition data between two fixed observers."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    left_obs = validate_surjective_map(C_left, h.shape[0], tol)
    right_obs = validate_surjective_map(C_right, h.shape[0], tol)
    if left_obs.shape[0] != right_obs.shape[0]:
        raise InputValidationError("C_left and C_right must have the same visible rank")

    left = fixed_observer_coordinates(h, left_obs, tolerances=tol)
    right = fixed_observer_coordinates(h, right_obs, tolerances=tol)
    transform = np.linalg.inv(left.adapted_basis) @ right.adapted_basis
    visible_dim = left_obs.shape[0]
    a = transform[:visible_dim, :visible_dim]
    b = transform[:visible_dim, visible_dim:]
    c = transform[visible_dim:, :visible_dim]
    d = transform[visible_dim:, visible_dim:]
    c_hat = c + left.coupling.T @ a
    d_hat = d + left.coupling.T @ b

    right_hidden_formula = symmetrize(b.T @ left.phi @ b + d_hat.T @ left.hidden_block @ d_hat)
    right_coupling_formula = (a.T @ left.phi @ b + c_hat.T @ left.hidden_block @ d_hat) @ np.linalg.inv(right_hidden_formula)
    right_phi_formula = symmetrize(
        a.T @ left.phi @ a
        + c_hat.T @ left.hidden_block @ c_hat
        - (a.T @ left.phi @ b + c_hat.T @ left.hidden_block @ d_hat) @ np.linalg.inv(right_hidden_formula) @ (a.T @ left.phi @ b + c_hat.T @ left.hidden_block @ d_hat).T
    )
    residual = max(
        float(np.linalg.norm(right.phi - right_phi_formula, ord="fro")),
        float(np.linalg.norm(right.hidden_block - right_hidden_formula, ord="fro")),
        float(np.linalg.norm(right.coupling - right_coupling_formula, ord="fro")),
    )
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="observer-transition-law",
        ambient_dim=h.shape[0],
        support_rank=h.shape[0] - visible_dim,
        visible_dim=visible_dim,
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=_condition_notes(h),
    )
    return ObserverTransitionResult(
        left=left,
        right=right,
        transform=transform,
        a=a,
        b=b,
        c=c,
        d=d,
        c_hat=c_hat,
        d_hat=d_hat,
        right_phi_from_left=right_phi_formula,
        right_hidden_from_left=right_hidden_formula,
        right_coupling_from_left=right_coupling_formula,
        residual=residual,
        metadata=metadata,
    )


def connection_current(
    H: np.ndarray,
    C: np.ndarray,
    Delta: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ConnectionCurrentResult:
    """Return K-velocity, current J, and forcing Q for one tangent direction."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    c = validate_surjective_map(C, h.shape[0], tol)
    delta = validate_symmetric_matrix("Delta", Delta, tol)
    if delta.shape != h.shape:
        raise InputValidationError("Delta must have the same shape as H")

    chart = fixed_observer_coordinates(h, c, tolerances=tol)
    delta_tilde = chart.adapted_basis.T @ delta @ chart.adapted_basis
    visible_dim = c.shape[0]
    delta_vh = delta_tilde[:visible_dim, visible_dim:]
    delta_hh = delta_tilde[visible_dim:, visible_dim:]
    coupling_velocity = (delta_vh - chart.coupling @ delta_hh) @ np.linalg.inv(chart.hidden_block)
    current = np.linalg.inv(chart.phi) @ coupling_velocity @ chart.hidden_block
    forcing = symmetrize(chart.phi @ current @ np.linalg.inv(chart.hidden_block) @ current.T @ chart.phi)
    local = local_visible_calculus(h, c, delta, tolerances=tol)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="connection-current",
        ambient_dim=h.shape[0],
        support_rank=chart.hidden_block.shape[0],
        visible_dim=visible_dim,
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=_condition_notes(h),
    )
    return ConnectionCurrentResult(
        phi=chart.phi,
        hidden_block=chart.hidden_block,
        coupling=chart.coupling,
        coupling_velocity=coupling_velocity,
        current=current,
        forcing=forcing,
        q=local.Q,
        adapted_basis=chart.adapted_basis,
        metadata=metadata,
    )


def forcing_from_current(
    phi: np.ndarray,
    hidden_block: np.ndarray,
    current: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return the exact forcing term Phi J R^{-1} J^T Phi."""
    tol = resolve_tolerances(tolerances)
    phi_matrix = validate_spd_matrix("phi", phi, tol)
    hidden_candidate = to_float_array("hidden_block", hidden_block)
    if hidden_candidate.shape == (0, 0):
        hidden = hidden_candidate
    else:
        hidden = validate_spd_matrix("hidden_block", hidden_candidate, tol)
    current_matrix = to_float_array("current", current)
    if current_matrix.shape != (phi_matrix.shape[0], hidden.shape[0]):
        raise InputValidationError("current has incompatible shape")
    if hidden.shape == (0, 0):
        return np.zeros((phi_matrix.shape[0], phi_matrix.shape[0]), dtype=float)
    return symmetrize(phi_matrix @ current_matrix @ np.linalg.inv(hidden) @ current_matrix.T @ phi_matrix)


def _fixed_adapted_basis(C: np.ndarray) -> np.ndarray:
    right_inverse = C.T @ np.linalg.inv(C @ C.T)
    null = la.null_space(C, rcond=1e-12)
    return np.concatenate([right_inverse, null], axis=1)


def _coordinates_from_adapted_basis(H: np.ndarray, visible_dim: int, adapted_basis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    h_tilde = adapted_basis.T @ H @ adapted_basis
    h_vv = h_tilde[:visible_dim, :visible_dim]
    h_vh = h_tilde[:visible_dim, visible_dim:]
    h_hh = h_tilde[visible_dim:, visible_dim:]
    phi = symmetrize(h_vv - h_vh @ np.linalg.solve(h_hh, h_vh.T))
    coupling = h_vh @ np.linalg.inv(h_hh)
    return phi, symmetrize(h_hh), coupling


def _validate_adapted_basis(adapted_basis: np.ndarray, ambient_dim: int, tolerances: Tolerances) -> np.ndarray:
    basis = to_float_array("adapted_basis", adapted_basis)
    if basis.shape != (ambient_dim, ambient_dim):
        raise InputValidationError("adapted_basis must be square with the ambient dimension")
    determinant = float(np.linalg.det(basis))
    cutoff = 10.0 * max(tolerances.atol, tolerances.rtol)
    if abs(determinant) <= cutoff:
        raise InputValidationError("adapted_basis must be invertible")
    return basis


def _condition_number_spd(H: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(H)
    return float(np.max(eigenvalues) / np.min(eigenvalues))


def _condition_notes(H: np.ndarray) -> tuple[str, ...]:
    condition = _condition_number_spd(H)
    if condition >= 1e9:
        return (f"H is extremely ill-conditioned (cond={condition:.3e}); fixed-observer coordinates may be numerically unreliable",)
    if condition >= 1e8:
        return (f"H is severely ill-conditioned (cond={condition:.3e}); interpret fixed-observer metric numerics with care",)
    if condition >= 1e6:
        return (f"H is highly conditioned (cond={condition:.3e}); fixed-observer diagnostics should be checked",)
    return ()
