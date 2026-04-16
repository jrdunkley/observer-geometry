"""
Exact local chart reduction — 0.3.3 Technical Note §6 (Theorem 6.1).

This module implements Layer 2 of the 5-layer stack: given ambient local
data (full Hessian, active face indices, orbit tangent directions), produce
a genuinely reduced ReducedLocalDatum on the transverse slice with correct
Jacobian and orbit-volume bookkeeping.

The key responsibilities are:

    1. Active-face restriction: project H onto the span of the active
       boundary face (if the optimum sits on a boundary stratum).
    2. Orbit tangent removal: given an explicit basis for the orbit
       tangent directions T_{G·v̂}, construct the orthogonal complement
       and project the active Hessian onto the transverse slice.
    3. Jacobian convention: record the slice-chart Jacobian so that the
       evidence dispatcher can include it as an external factor.
    4. Orbit-volume bookkeeping: carry log Vol(G·v̂) as metadata.
    5. Production of a genuinely transverse PSD Hessian on the active
       slice, ready for the quadratic regime detector.

The output is a ReducedLocalDatum with active_dim = transverse dimension
and h_active = transverse Hessian H_⊥.

References
----------
0.3.3 Technical Note:
  - §6  Definition 6.1 (regular quotient point)
  - §6  Theorem 6.1  (exact local slice reduction)
  - §7  Corollary 7.1 (typed local evidence template)
"""
from __future__ import annotations

from typing import Mapping, Any

import numpy as np
from numpy.typing import NDArray

from .exceptions import InputValidationError
from .regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    OrbitSpec,
    ReducedLocalDatum,
)
from .validation import Tolerances, resolve_tolerances

Array = NDArray[np.float64]


# ═══════════════════════════════════════════════════════════════════════
# Primary API: reduce_local_chart
# ═══════════════════════════════════════════════════════════════════════

def reduce_local_chart(
    H_ambient: Array,
    *,
    active_face_indices: Array | None = None,
    orbit_tangent_basis: Array | None = None,
    log_orbit_volume: float | None = None,
    log_slice_jacobian: float | None = None,
    jacobian_convention: JacobianConvention = JacobianConvention.SLICE_LEBESGUE,
    cone: ConeSpec | None = None,
    tolerances: Tolerances | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> ReducedLocalDatum:
    """Construct a genuinely reduced ReducedLocalDatum from ambient local data.

    This is the exact local chart reduction of Theorem 6.1: starting from
    the full ambient Hessian at an optimum, it performs active-face
    restriction and orbit-tangent removal, producing a transverse datum
    ready for the quadratic regime detector.

    Parameters
    ----------
    H_ambient : array, shape (d, d)
        Symmetric PSD Hessian of the visible action in ambient coordinates
        at the optimum.  This is the full Hessian before any reduction.
    active_face_indices : array of int, optional
        Indices of the ambient coordinates that are active (interior to
        the boundary face).  If None, all coordinates are active (no
        boundary).  If provided, H is first restricted to the active
        subspace before any orbit reduction.
    orbit_tangent_basis : array, shape (d_active, r), optional
        Columns are an orthonormal basis for the orbit tangent directions
        T_{G·v̂} at the optimum, expressed in the active coordinate
        system (after face restriction, if any).  None means no quotient.
        The columns must be orthonormal (verified within tolerance).
    log_orbit_volume : float, optional
        log Vol(G·v̂) — the orbit volume in the chosen Haar measure.
    log_slice_jacobian : float, optional
        log of the slice-chart Jacobian J_sl at the identity coset.
    jacobian_convention : JacobianConvention
        Which convention the Jacobian and orbit volume are stated in.
    cone : ConeSpec, optional
        Admissible tangent cone in the *ambient* coordinates.  If
        active_face_indices is provided, the cone is restricted to the
        active subspace.  If None, defaults to FULL_SPACE.
    tolerances : Tolerances, optional
        Numerical tolerances for orthonormality checks and PSD verification.
    metadata : dict, optional
        User-supplied metadata to carry through.

    Returns
    -------
    ReducedLocalDatum
        With active_dim = transverse dimension, h_active = transverse
        Hessian H_⊥, cone restricted to the transverse space, and
        orbit data properly separated.

    Notes
    -----
    The caller is responsible for providing a correct orbit tangent basis.
    This function verifies orthonormality and that the tangent directions
    lie in the kernel of H (within tolerance), but does not derive the
    orbit structure from the Hessian alone.
    """
    tol = resolve_tolerances(tolerances)
    H = np.asarray(H_ambient, dtype=float)

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise InputValidationError("H_ambient must be a square matrix")
    d = H.shape[0]

    # Symmetry check.
    if not np.allclose(H, H.T, atol=tol.atol, rtol=tol.rtol):
        raise InputValidationError("H_ambient must be symmetric")

    # ── Step 1: Active-face restriction ─────────────────────────────
    if active_face_indices is not None:
        idx = np.asarray(active_face_indices, dtype=int)
        if idx.ndim != 1:
            raise InputValidationError("active_face_indices must be 1-D")
        if len(idx) == 0:
            raise InputValidationError("active_face_indices must be non-empty")
        if np.any(idx < 0) or np.any(idx >= d):
            raise InputValidationError(
                f"active_face_indices out of range [0, {d})"
            )
        # Restrict H to the active subspace.
        H_active = H[np.ix_(idx, idx)].copy()
        d_active = len(idx)
    else:
        H_active = H.copy()
        d_active = d

    # ── Step 2: Cone handling ───────────────────────────────────────
    if cone is None:
        if active_face_indices is not None:
            # The face restriction already handled the boundary.
            # The active subspace is the interior of the face → FULL_SPACE.
            active_cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d_active)
        else:
            active_cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d_active)
    else:
        if active_face_indices is not None:
            # Restrict the cone to the active subspace.
            active_cone = _restrict_cone_to_face(cone, active_face_indices, d_active)
        else:
            active_cone = cone

    # ── Step 3: Orbit tangent removal ───────────────────────────────
    if orbit_tangent_basis is None:
        # No quotient — the active datum is already the final datum.
        return ReducedLocalDatum(
            active_dim=d_active,
            h_active=H_active,
            cone=active_cone,
            orbit=None,
            metadata=metadata,
        )

    V = np.asarray(orbit_tangent_basis, dtype=float)
    if V.ndim != 2:
        raise InputValidationError("orbit_tangent_basis must be 2-D")
    if V.shape[0] != d_active:
        raise InputValidationError(
            f"orbit_tangent_basis has {V.shape[0]} rows, expected {d_active}"
        )
    r = V.shape[1]  # orbit dimension
    if r >= d_active:
        raise InputValidationError(
            f"orbit_dim={r} >= active_dim={d_active}: no transverse directions"
        )

    # Verify orthonormality.
    gram = V.T @ V
    if not np.allclose(gram, np.eye(r), atol=tol.atol, rtol=tol.rtol):
        raise InputValidationError(
            "orbit_tangent_basis columns must be orthonormal"
        )

    # Verify orbit tangent directions lie in (or near) the kernel of H_active.
    # This is the mathematical requirement for a removable symmetry:
    # H · v_orbit = 0 for each orbit tangent direction.
    HV = H_active @ V
    hv_norm = float(np.linalg.norm(HV))
    h_norm = float(np.linalg.norm(H_active))
    if h_norm > tol.atol and hv_norm > tol.atol + tol.rtol * h_norm:
        raise InputValidationError(
            f"orbit tangent directions are not in the kernel of H_active: "
            f"||H @ V|| = {hv_norm:.2e}, ||H|| = {h_norm:.2e}.  "
            f"For a removable quotient, the orbit tangent must lie in "
            f"ker(H_active)."
        )

    # Construct the orthogonal complement of the orbit tangent space.
    # Q_perp has shape (d_active, m) where m = d_active - r.
    Q_perp = _orthogonal_complement(V, d_active, tol)
    m = d_active - r  # transverse dimension

    # Project H onto the transverse subspace.
    H_perp = Q_perp.T @ H_active @ Q_perp  # shape (m, m)

    # Symmetrise to kill numerical asymmetry.
    H_perp = 0.5 * (H_perp + H_perp.T)

    # Restrict cone to the transverse subspace.
    transverse_cone = _restrict_cone_to_transverse(active_cone, Q_perp, m)

    # Build OrbitSpec.
    orbit = OrbitSpec(
        orbit_dim=r,
        jacobian_convention=jacobian_convention,
        log_orbit_volume=log_orbit_volume,
        log_slice_jacobian=log_slice_jacobian,
    )

    return ReducedLocalDatum(
        active_dim=m,
        h_active=H_perp,
        cone=transverse_cone,
        orbit=orbit,
        metadata=metadata,
    )


# ═══════════════════════════════════════════════════════════════════════
# Active-face restriction for cones
# ═══════════════════════════════════════════════════════════════════════

def active_face_restriction(
    H_ambient: Array,
    active_face_indices: Array,
    tolerances: Tolerances | None = None,
) -> tuple[Array, int]:
    """Restrict the ambient Hessian to an active boundary face.

    Returns the restricted Hessian and the face dimension.  This is a
    simpler entry point when only face restriction is needed (no orbit).

    Parameters
    ----------
    H_ambient : array, shape (d, d)
    active_face_indices : array of int
        Indices of the active (interior-to-face) coordinates.

    Returns
    -------
    (H_face, face_dim) : tuple
    """
    tol = resolve_tolerances(tolerances)
    H = np.asarray(H_ambient, dtype=float)
    idx = np.asarray(active_face_indices, dtype=int)

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise InputValidationError("H_ambient must be a square matrix")
    d = H.shape[0]
    if idx.ndim != 1 or len(idx) == 0:
        raise InputValidationError("active_face_indices must be a non-empty 1-D array")
    if np.any(idx < 0) or np.any(idx >= d):
        raise InputValidationError(f"active_face_indices out of range [0, {d})")

    H_face = H[np.ix_(idx, idx)].copy()
    return H_face, len(idx)


# ═══════════════════════════════════════════════════════════════════════
# Transverse complement construction
# ═══════════════════════════════════════════════════════════════════════

def transverse_complement(
    orbit_tangent_basis: Array,
    ambient_dim: int,
    tolerances: Tolerances | None = None,
) -> Array:
    """Compute an orthonormal basis for the transverse complement.

    Given an orthonormal basis V for the orbit tangent space, return Q_⊥
    such that [V | Q_⊥] spans the full ambient space and Q_⊥ is
    orthonormal.

    Parameters
    ----------
    orbit_tangent_basis : array, shape (d, r)
        Orthonormal columns spanning the orbit tangent space.
    ambient_dim : int
        Dimension of the ambient space.

    Returns
    -------
    Q_perp : array, shape (d, d-r)
        Orthonormal columns spanning the transverse complement.
    """
    tol = resolve_tolerances(tolerances)
    V = np.asarray(orbit_tangent_basis, dtype=float)
    if V.ndim != 2 or V.shape[0] != ambient_dim:
        raise InputValidationError(
            f"orbit_tangent_basis shape {V.shape} incompatible with "
            f"ambient_dim={ambient_dim}"
        )
    return _orthogonal_complement(V, ambient_dim, tol)


# ═══════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════

def _orthogonal_complement(V: Array, d: int, tol: Tolerances) -> Array:
    """Orthonormal basis for the complement of span(V) in R^d.

    Uses a full QR of [V | I] projected out of V, or equivalently
    the null space of V^T.
    """
    r = V.shape[1]
    m = d - r
    if m == 0:
        return np.empty((d, 0), dtype=float)

    # Project the identity onto the complement of V.
    # P_perp = I - V V^T
    P_perp = np.eye(d) - V @ V.T
    # QR on P_perp to extract an orthonormal basis for the range.
    Q, R = np.linalg.qr(P_perp, mode="complete")
    # The first m columns of Q with nonzero diagonal in R
    # span the complement.  But we need to be careful: the
    # diagonal of R tells us which columns are in the range.
    diag_R = np.abs(np.diag(R[:d, :]))
    threshold = max(tol.atol, tol.rtol * float(np.max(diag_R)))
    mask = diag_R > threshold
    Q_perp = Q[:, mask]

    if Q_perp.shape[1] != m:
        # Fallback: SVD of V^T to get the null space.
        _, S, Vt = np.linalg.svd(V.T, full_matrices=True)
        Q_perp = Vt[r:].T  # shape (d, m)

    return Q_perp


def _restrict_cone_to_face(
    cone: ConeSpec,
    face_indices: Array,
    face_dim: int,
) -> ConeSpec:
    """Restrict a cone to the active face subspace."""
    if cone.kind == ConeKind.FULL_SPACE:
        return ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=face_dim)

    if cone.kind == ConeKind.ABSTRACT:
        return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=face_dim)

    if cone.kind == ConeKind.ORTHANT:
        # Restrict orthant signs to the active indices.
        if cone.orthant_signs is not None:
            restricted_signs = cone.orthant_signs[face_indices]
            return ConeSpec(
                kind=ConeKind.ORTHANT,
                ambient_dim=face_dim,
                orthant_signs=restricted_signs,
            )
        return ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=face_dim)

    if cone.kind == ConeKind.POLYHEDRAL:
        if cone.halfspace_normals is not None:
            # Restrict normals to active coordinates.
            A_face = cone.halfspace_normals[:, face_indices]
            # Drop rows that became all-zero (constraints on inactive coords).
            norms = np.linalg.norm(A_face, axis=1)
            active_rows = norms > 1e-15
            if np.any(active_rows):
                A_face = A_face[active_rows]
                return ConeSpec(
                    kind=ConeKind.POLYHEDRAL,
                    ambient_dim=face_dim,
                    halfspace_normals=A_face,
                )
        return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=face_dim)

    return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=face_dim)


def _restrict_cone_to_transverse(
    cone: ConeSpec,
    Q_perp: Array,
    transverse_dim: int,
) -> ConeSpec:
    """Restrict a cone to the transverse subspace via projection.

    For polyhedral cones {x : Ax >= 0}, the transverse restriction is
    {u : A Q_⊥ u >= 0} where Q_⊥ maps transverse coords to active coords.
    """
    if cone.kind == ConeKind.FULL_SPACE:
        return ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=transverse_dim)

    if cone.kind == ConeKind.ABSTRACT:
        return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=transverse_dim)

    if cone.kind in (ConeKind.ORTHANT, ConeKind.POLYHEDRAL):
        # For orthant or polyhedral cones, the transverse restriction
        # is in general a polyhedral cone.  Build the projected normals.
        if cone.kind == ConeKind.ORTHANT:
            d = cone.ambient_dim
            if cone.orthant_signs is not None:
                A = np.diag(cone.orthant_signs.astype(float))
            else:
                A = np.eye(d)
        else:
            if cone.halfspace_normals is None:
                return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=transverse_dim)
            A = cone.halfspace_normals

        A_perp = A @ Q_perp  # shape (num_constraints, transverse_dim)
        # Drop zero rows.
        norms = np.linalg.norm(A_perp, axis=1)
        active = norms > 1e-15
        if np.any(active):
            return ConeSpec(
                kind=ConeKind.POLYHEDRAL,
                ambient_dim=transverse_dim,
                halfspace_normals=A_perp[active],
                metadata={"derived_from": "transverse_cone_restriction"},
            )
        return ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=transverse_dim)

    return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=transverse_dim)
