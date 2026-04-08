from __future__ import annotations

import numpy as np

from .types import LinearAlgebraMetadata, LocalCalculusResult, VisibleGeometryResult
from .validation import (
    Tolerances,
    range_basis,
    resolve_tolerances,
    solve_spd,
    symmetrize,
    validate_spd_matrix,
    validate_surjective_map,
    validate_symmetric_matrix,
)


def visible_geometry(H: np.ndarray, C: np.ndarray, tolerances: Tolerances | None = None) -> VisibleGeometryResult:
    """Return Phi_C(H), the canonical lift, and the hidden projector."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    c = validate_surjective_map(C, h.shape[0], tol)

    h_inv_ct = solve_spd(h, c.T)
    gram = symmetrize(c @ h_inv_ct)
    phi = solve_spd(gram, np.eye(gram.shape[0], dtype=float))
    phi = symmetrize(phi)
    lift = h_inv_ct @ phi
    projector = symmetrize(np.eye(h.shape[0], dtype=float) - lift @ c)
    projector = np.eye(h.shape[0], dtype=float) - lift @ c

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="cholesky",
        ambient_dim=h.shape[0],
        support_rank=h.shape[0],
        visible_dim=c.shape[0],
        support_restricted=False,
    )
    return VisibleGeometryResult(phi=phi, lift=lift, projector=projector, metadata=metadata)


def visible_precision(H: np.ndarray, C: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Compute Phi_C(H) = (C H^{-1} C^T)^{-1}."""
    return visible_geometry(H, C, tolerances=tolerances).phi


def canonical_lift(H: np.ndarray, C: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Compute L_{C,H} = H^{-1} C^T Phi_C(H)."""
    return visible_geometry(H, C, tolerances=tolerances).lift


def hidden_projector(H: np.ndarray, C: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Compute P_{C,H} = I - L_{C,H} C."""
    return visible_geometry(H, C, tolerances=tolerances).projector


def local_visible_calculus(
    H: np.ndarray,
    C: np.ndarray,
    Delta: np.ndarray,
    tolerances: Tolerances | None = None,
) -> LocalCalculusResult:
    """Compute the first visible response V and the positive quartic defect Q."""
    tol = resolve_tolerances(tolerances)
    geometry = visible_geometry(H, C, tolerances=tol)
    delta = validate_symmetric_matrix("Delta", Delta, tol)
    h = validate_spd_matrix("H", H, tol)

    V = symmetrize(geometry.lift.T @ delta @ geometry.lift)
    temp = geometry.projector.T @ delta @ geometry.lift
    Q = symmetrize(temp.T @ solve_spd(h, temp))

    phi_inv_v = solve_spd(geometry.phi, V)
    phi_inv_q = solve_spd(geometry.phi, Q)
    det_split = float(np.trace(phi_inv_v @ phi_inv_v) + 2.0 * np.trace(phi_inv_q))
    active_support = range_basis(V, tol)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="cholesky",
        ambient_dim=h.shape[0],
        support_rank=active_support.shape[1],
        visible_dim=geometry.phi.shape[0],
        support_restricted=False,
    )
    return LocalCalculusResult(
        phi=geometry.phi,
        lift=geometry.lift,
        projector=geometry.projector,
        V=V,
        Q=Q,
        det_split=det_split,
        active_support=active_support,
        metadata=metadata,
    )

