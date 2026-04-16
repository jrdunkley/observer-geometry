"""Tests for exact local chart reduction (slice.py, 0.3.3)."""
import numpy as np
import pytest

from nomogeo.slice import (
    reduce_local_chart,
    active_face_restriction,
    transverse_complement,
)
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    ReducedLocalDatum,
)
from nomogeo.exceptions import InputValidationError


def _spd(d: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(d, d))
    return raw.T @ raw + 0.5 * np.eye(d)


def _psd_with_kernel(d: int, kernel_dim: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    eigs = np.zeros(d)
    eigs[:d - kernel_dim] = rng.uniform(0.5, 3.0, size=d - kernel_dim)
    return Q @ np.diag(eigs) @ Q.T


# ═══════════════════════════════════════════════════════════════════════
# No reduction (passthrough)
# ═══════════════════════════════════════════════════════════════════════

def test_no_reduction_returns_full_datum() -> None:
    """Without face restriction or orbit, get back the full Hessian."""
    H = _spd(4, seed=1)
    datum = reduce_local_chart(H)
    assert isinstance(datum, ReducedLocalDatum)
    assert datum.active_dim == 4
    assert datum.h_active.shape == (4, 4)
    assert datum.orbit is None
    assert datum.cone.kind == ConeKind.FULL_SPACE
    assert np.allclose(datum.h_active, H)


# ═══════════════════════════════════════════════════════════════════════
# Active-face restriction
# ═══════════════════════════════════════════════════════════════════════

def test_face_restriction_submatrix() -> None:
    """Face restriction should give the correct submatrix."""
    H = _spd(5, seed=2)
    idx = np.array([1, 3, 4])
    datum = reduce_local_chart(H, active_face_indices=idx)
    assert datum.active_dim == 3
    expected = H[np.ix_(idx, idx)]
    assert np.allclose(datum.h_active, expected)


def test_face_restriction_standalone() -> None:
    """active_face_restriction utility function."""
    H = _spd(4, seed=3)
    idx = np.array([0, 2])
    H_face, dim = active_face_restriction(H, idx)
    assert dim == 2
    assert np.allclose(H_face, H[np.ix_(idx, idx)])


def test_face_restriction_preserves_symmetry() -> None:
    H = _spd(6, seed=4)
    idx = np.array([0, 2, 4])
    datum = reduce_local_chart(H, active_face_indices=idx)
    assert np.allclose(datum.h_active, datum.h_active.T)


# ═══════════════════════════════════════════════════════════════════════
# Orbit tangent removal
# ═══════════════════════════════════════════════════════════════════════

def test_orbit_removal_reduces_dimension() -> None:
    """Removing a 1-dim orbit from a 4-dim space gives 3-dim transverse."""
    # Build a PSD matrix with a known kernel direction (orbit tangent).
    H = _psd_with_kernel(4, kernel_dim=1, seed=5)
    evals, evecs = np.linalg.eigh(H)
    # The kernel eigenvector is the orbit tangent.
    kernel_idx = np.argmin(evals)
    orbit_basis = evecs[:, kernel_idx:kernel_idx + 1]

    datum = reduce_local_chart(
        H,
        orbit_tangent_basis=orbit_basis,
        log_orbit_volume=1.0,
        log_slice_jacobian=0.0,
    )
    assert datum.active_dim == 3
    assert datum.h_active.shape == (3, 3)
    assert datum.orbit is not None
    assert datum.orbit.orbit_dim == 1
    assert datum.orbit.log_orbit_volume == pytest.approx(1.0)


def test_orbit_removal_preserves_nonzero_eigenvalues() -> None:
    """The transverse Hessian should have the same nonzero eigenvalues
    as the original, minus the zero eigenvalue from the orbit."""
    H = _psd_with_kernel(4, kernel_dim=1, seed=6)
    evals_orig = np.sort(np.linalg.eigvalsh(H))[::-1]
    positive_evals = evals_orig[evals_orig > 1e-10]

    kernel_idx = np.argmin(np.linalg.eigvalsh(H))
    evecs = np.linalg.eigh(H)[1]
    orbit_basis = evecs[:, kernel_idx:kernel_idx + 1]

    datum = reduce_local_chart(H, orbit_tangent_basis=orbit_basis)
    evals_transverse = np.sort(np.linalg.eigvalsh(datum.h_active))[::-1]

    assert len(evals_transverse) == len(positive_evals)
    assert np.allclose(evals_transverse, positive_evals, atol=1e-10)


def test_orbit_removal_invariant_under_basis_rotation() -> None:
    """Different orthonormal bases for the same orbit tangent space
    should give the same transverse Hessian eigenvalues."""
    H = _psd_with_kernel(5, kernel_dim=2, seed=7)
    evals, evecs = np.linalg.eigh(H)
    kernel_mask = evals < 1e-10
    V1 = evecs[:, kernel_mask]

    # Rotate the orbit basis.
    rng = np.random.default_rng(42)
    R = np.linalg.qr(rng.normal(size=(2, 2)))[0]
    V2 = V1 @ R

    datum1 = reduce_local_chart(H, orbit_tangent_basis=V1)
    datum2 = reduce_local_chart(H, orbit_tangent_basis=V2)

    evals1 = np.sort(np.linalg.eigvalsh(datum1.h_active))
    evals2 = np.sort(np.linalg.eigvalsh(datum2.h_active))
    assert np.allclose(evals1, evals2, atol=1e-10)


def test_orbit_removal_rejects_non_kernel_tangent() -> None:
    """If the orbit tangent is not in ker(H), raise an error."""
    H = _spd(3, seed=8)  # SPD → no kernel
    v = np.array([[1.0], [0.0], [0.0]])  # not in ker(H)
    with pytest.raises(InputValidationError, match="kernel"):
        reduce_local_chart(H, orbit_tangent_basis=v)


def test_orbit_removal_rejects_non_orthonormal() -> None:
    """Non-orthonormal orbit basis should be rejected."""
    H = _psd_with_kernel(3, kernel_dim=1, seed=9)
    v = np.array([[2.0], [0.0], [0.0]])  # not unit norm
    with pytest.raises(InputValidationError, match="orthonormal"):
        reduce_local_chart(H, orbit_tangent_basis=v)


# ═══════════════════════════════════════════════════════════════════════
# Combined face + orbit reduction
# ═══════════════════════════════════════════════════════════════════════

def test_face_then_orbit() -> None:
    """Face restriction followed by orbit removal."""
    d = 6
    H_full = _psd_with_kernel(d, kernel_dim=1, seed=10)

    # Restrict to face {0,1,2,3} (4-dim), then remove 1-dim orbit.
    face_idx = np.array([0, 1, 2, 3])
    H_face = H_full[np.ix_(face_idx, face_idx)]

    # Find kernel of H_face.
    evals, evecs = np.linalg.eigh(H_face)
    kernel_idx = np.argmin(evals)
    orbit_basis = evecs[:, kernel_idx:kernel_idx + 1]

    # Only proceed if the kernel direction is actually near-zero.
    if evals[kernel_idx] < 1e-10:
        datum = reduce_local_chart(
            H_full,
            active_face_indices=face_idx,
            orbit_tangent_basis=orbit_basis,
        )
        assert datum.active_dim == 3
        assert datum.orbit is not None
        assert datum.orbit.orbit_dim == 1


# ═══════════════════════════════════════════════════════════════════════
# Transverse complement utility
# ═══════════════════════════════════════════════════════════════════════

def test_transverse_complement_dimensions() -> None:
    rng = np.random.default_rng(11)
    V = np.linalg.qr(rng.normal(size=(5, 2)))[0][:, :2]
    Q = transverse_complement(V, 5)
    assert Q.shape == (5, 3)
    # Verify orthonormality of Q.
    assert np.allclose(Q.T @ Q, np.eye(3), atol=1e-10)
    # Verify orthogonality with V.
    assert np.allclose(V.T @ Q, 0, atol=1e-10)


def test_transverse_complement_spans_full_space() -> None:
    rng = np.random.default_rng(12)
    V = np.linalg.qr(rng.normal(size=(4, 1)))[0][:, :1]
    Q = transverse_complement(V, 4)
    combined = np.hstack([V, Q])
    assert np.linalg.matrix_rank(combined) == 4


# ═══════════════════════════════════════════════════════════════════════
# Cone handling in reduction
# ═══════════════════════════════════════════════════════════════════════

def test_cone_restriction_to_face() -> None:
    """Orthant cone restricted to a face."""
    H = _spd(4, seed=13)
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=4)
    idx = np.array([0, 2])
    datum = reduce_local_chart(H, active_face_indices=idx, cone=cone)
    assert datum.active_dim == 2
    # The restricted cone should still be an orthant on the active face.
    assert datum.cone.kind == ConeKind.ORTHANT


def test_cone_transverse_becomes_polyhedral() -> None:
    """Orthant cone projected to transverse space becomes polyhedral."""
    H = _psd_with_kernel(3, kernel_dim=1, seed=14)
    evals, evecs = np.linalg.eigh(H)
    kernel_idx = np.argmin(evals)
    orbit_basis = evecs[:, kernel_idx:kernel_idx + 1]

    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)
    datum = reduce_local_chart(H, orbit_tangent_basis=orbit_basis, cone=cone)
    assert datum.active_dim == 2
    # The transverse cone should be polyhedral (projected orthant).
    assert datum.cone.kind == ConeKind.POLYHEDRAL


# ═══════════════════════════════════════════════════════════════════════
# Input validation
# ═══════════════════════════════════════════════════════════════════════

def test_rejects_non_square() -> None:
    H = np.ones((3, 4))
    with pytest.raises(InputValidationError, match="square"):
        reduce_local_chart(H)


def test_rejects_non_symmetric() -> None:
    H = np.array([[1.0, 2.0], [0.0, 1.0]])
    with pytest.raises(InputValidationError, match="symmetric"):
        reduce_local_chart(H)


def test_rejects_face_indices_out_of_range() -> None:
    H = _spd(3, seed=15)
    with pytest.raises(InputValidationError, match="range"):
        reduce_local_chart(H, active_face_indices=np.array([0, 5]))


def test_rejects_orbit_dim_too_large() -> None:
    """orbit_dim >= active_dim should fail."""
    H = np.zeros((2, 2))
    V = np.eye(2)
    with pytest.raises(InputValidationError, match="transverse"):
        reduce_local_chart(H, orbit_tangent_basis=V)
