"""Tests for the quadratic regime detector (regime.py, 0.3.3)."""
import numpy as np
import pytest

from nomogeo.regime import classify_regime, classify_from_hessian
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    OrbitSpec,
    ReducedLocalDatum,
    RegimeKind,
    RegularCone,
    RegularInterior,
    RegularQuotient,
    RegularQuotientCone,
    IndefiniteStationaryPoint,
    UnresolvedKernel,
    UnsupportedConeReduction,
)
from nomogeo.exceptions import InputValidationError


def _spd(d: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(d, d))
    return raw.T @ raw + 0.5 * np.eye(d)


def _psd_with_kernel(d: int, kernel_dim: int, seed: int = 0) -> np.ndarray:
    """Create a PSD matrix of size d with exactly kernel_dim zero eigenvalues."""
    rng = np.random.default_rng(seed)
    # Build in the eigenbasis.
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    eigs = np.zeros(d)
    eigs[:d - kernel_dim] = rng.uniform(0.5, 3.0, size=d - kernel_dim)
    return Q @ np.diag(eigs) @ Q.T


def _datum(H, cone=None, orbit=None):
    d = H.shape[0]
    if cone is None:
        cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d)
    return ReducedLocalDatum(active_dim=d, h_active=H, cone=cone, orbit=orbit)


# ═══════════════════════════════════════════════════════════════════════
# Regular Interior
# ═══════════════════════════════════════════════════════════════════════

def test_regular_interior_from_spd() -> None:
    H = _spd(3, seed=1)
    result = classify_regime(_datum(H))
    assert isinstance(result, RegularInterior)
    assert result.kind == RegimeKind.REGULAR_INTERIOR
    expected_logdet = float(np.linalg.slogdet(H)[1])
    assert result.log_det_h_active == pytest.approx(expected_logdet)


def test_regular_interior_1d() -> None:
    H = np.array([[2.5]])
    result = classify_from_hessian(H)
    assert isinstance(result, RegularInterior)
    assert result.log_det_h_active == pytest.approx(np.log(2.5))


# ═══════════════════════════════════════════════════════════════════════
# Regular Cone
# ═══════════════════════════════════════════════════════════════════════

def test_regular_cone_orthant() -> None:
    H = _spd(3, seed=2)
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, RegularCone)
    assert result.kind == RegimeKind.REGULAR_CONE


def test_regular_cone_abstract() -> None:
    H = _spd(2, seed=3)
    cone = ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=2)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, RegularCone)


def test_regular_cone_full_dim_polyhedral() -> None:
    """A polyhedral cone {x : x_1 >= 0} in R^2 is full-dimensional."""
    H = _spd(2, seed=4)
    A = np.array([[1.0, 0.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=2, halfspace_normals=A)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, RegularCone)


# ═══════════════════════════════════════════════════════════════════════
# Regular Quotient
# ═══════════════════════════════════════════════════════════════════════

def test_regular_quotient() -> None:
    """Quotient with genuinely reduced transverse datum.
    If the ambient space is 4-dim and orbit_dim=1, the transverse
    Hessian is 3×3 (already slice-reduced by the caller)."""
    H_transverse = _spd(3, seed=5)  # 3×3 transverse Hessian
    orbit = OrbitSpec(orbit_dim=1, jacobian_convention=JacobianConvention.NONE)
    result = classify_regime(_datum(H_transverse, orbit=orbit))
    assert isinstance(result, RegularQuotient)
    assert result.orbit_dim == 1
    assert result.datum.active_dim == 3  # transverse dim, not ambient


# ═══════════════════════════════════════════════════════════════════════
# Regular Quotient-Cone
# ═══════════════════════════════════════════════════════════════════════

def test_regular_quotient_cone() -> None:
    """Quotient-cone with genuinely reduced transverse datum.
    If ambient is 5-dim, orbit_dim=2, the transverse Hessian is 3×3."""
    H_transverse = _spd(3, seed=6)  # 3×3 transverse Hessian
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)  # cone on transverse space
    orbit = OrbitSpec(orbit_dim=2, jacobian_convention=JacobianConvention.SLICE_LEBESGUE)
    result = classify_regime(_datum(H_transverse, cone=cone, orbit=orbit))
    assert isinstance(result, RegularQuotientCone)
    assert result.orbit_dim == 2
    assert result.datum.active_dim == 3  # transverse dim


# ═══════════════════════════════════════════════════════════════════════
# Unresolved Kernel
# ═══════════════════════════════════════════════════════════════════════

def test_unresolved_kernel_rank_deficient() -> None:
    H = _psd_with_kernel(4, kernel_dim=1, seed=7)
    result = classify_regime(_datum(H))
    assert isinstance(result, UnresolvedKernel)
    assert result.kernel.kernel_dim == 1
    assert result.kernel.positive_normal_dim == 3
    assert result.kernel.kernel_basis.shape == (4, 1)
    assert result.kernel.positive_normal_basis.shape == (4, 3)
    assert np.all(result.kernel.positive_normal_eigenvalues > 0)


def test_unresolved_kernel_dim_2() -> None:
    H = _psd_with_kernel(5, kernel_dim=2, seed=8)
    result = classify_regime(_datum(H))
    assert isinstance(result, UnresolvedKernel)
    assert result.kernel.kernel_dim == 2
    assert result.kernel.positive_normal_dim == 3


def test_unresolved_kernel_normal_prefactor() -> None:
    """The normal prefactor should match the analytic formula."""
    H = _psd_with_kernel(4, kernel_dim=1, seed=9)
    result = classify_regime(_datum(H))
    assert isinstance(result, UnresolvedKernel)
    lam = result.kernel.positive_normal_eigenvalues
    expected = 0.5 * len(lam) * np.log(2.0 * np.pi) - 0.5 * float(np.sum(np.log(lam)))
    assert result.kernel.log_normal_prefactor == pytest.approx(expected)


def test_zero_matrix_is_fully_kernel() -> None:
    H = np.zeros((3, 3))
    result = classify_regime(_datum(H))
    assert isinstance(result, UnresolvedKernel)
    assert result.kernel.kernel_dim == 3
    assert result.kernel.positive_normal_dim == 0


# ═══════════════════════════════════════════════════════════════════════
# Unsupported Cone Reduction
# ═══════════════════════════════════════════════════════════════════════

def test_non_full_span_polyhedral_returns_unsupported() -> None:
    """A polyhedral cone {x : x_1 >= 0, -x_1 >= 0} = {x : x_1 = 0}
    has span dimension d-1, which is not supported."""
    H = _spd(3, seed=10)
    A = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=3, halfspace_normals=A)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, UnsupportedConeReduction)
    assert result.detected_span_dim == 2
    assert result.active_dim == 3


def test_counterexample_rank_vs_span() -> None:
    """The rank(A) bug counterexample: A=[[1,0]] has rank 1, but
    the cone {x : x_1 >= 0} is full-dimensional in R^2."""
    H = _spd(2, seed=11)
    A = np.array([[1.0, 0.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=2, halfspace_normals=A)
    result = classify_regime(_datum(H, cone=cone))
    # This MUST be RegularCone, not UnsupportedConeReduction.
    assert isinstance(result, RegularCone)


# ═══════════════════════════════════════════════════════════════════════
# Input validation
# ═══════════════════════════════════════════════════════════════════════

def test_rejects_non_symmetric() -> None:
    H = np.array([[1.0, 2.0], [0.0, 1.0]])
    with pytest.raises(InputValidationError, match="symmetric"):
        classify_regime(_datum(H))


def test_rejects_negative_definite() -> None:
    """Negative-definite Hessian → IndefiniteStationaryPoint (not an exception).

    The three-way spectral pre-check catches this before regime
    classification and returns a typed result instead of raising.
    """
    H = -np.eye(2)
    result = classify_regime(_datum(H))
    assert isinstance(result, IndefiniteStationaryPoint)
    assert result.kind == RegimeKind.INDEFINITE_STATIONARY
    assert result.negative_dim == 2
    assert result.positive_dim == 0
    assert result.min_eigenvalue < 0
    assert result.suggested_action in ("refit", "discard")


def test_rejects_mismatched_dim() -> None:
    H = np.eye(3)
    datum = ReducedLocalDatum(
        active_dim=2, h_active=H,
        cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=3),
    )
    with pytest.raises(InputValidationError, match="active_dim"):
        classify_regime(datum)


# ═══════════════════════════════════════════════════════════════════════
# Polyhedral counterexamples
# ═══════════════════════════════════════════════════════════════════════

def test_polyhedral_full_span_multiple_constraints() -> None:
    """A = [[1,0],[0,1]] in R^2 — positive orthant, full-dimensional."""
    H = _spd(2, seed=12)
    A = np.array([[1.0, 0.0], [0.0, 1.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=2, halfspace_normals=A)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, RegularCone)


def test_polyhedral_one_implicit_equality() -> None:
    """A = [[1,0,0],[0,1,0],[-1,0,0]] in R^3 — x_1=0 is an implicit
    equality, so span dim = 2."""
    H = _spd(3, seed=13)
    A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=3, halfspace_normals=A)
    result = classify_regime(_datum(H, cone=cone))
    assert isinstance(result, UnsupportedConeReduction)
    assert result.detected_span_dim == 2


def test_polyhedral_missing_normals_raises() -> None:
    H = _spd(2, seed=14)
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=2)
    with pytest.raises(InputValidationError, match="halfspace_normals"):
        classify_regime(_datum(H, cone=cone))


# ═══════════════════════════════════════════════════════════════════════
# classify_from_hessian convenience
# ═══════════════════════════════════════════════════════════════════════

def test_classify_from_hessian_defaults_to_full_space() -> None:
    H = _spd(3, seed=15)
    result = classify_from_hessian(H)
    assert isinstance(result, RegularInterior)


def test_classify_from_hessian_with_cone() -> None:
    H = _spd(3, seed=16)
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)
    result = classify_from_hessian(H, cone=cone)
    assert isinstance(result, RegularCone)


def test_classify_from_hessian_with_orbit() -> None:
    """When using classify_from_hessian with orbit, H must be the
    transverse Hessian (already slice-reduced)."""
    H_transverse = _spd(3, seed=17)  # 3×3 = transverse Hessian
    orbit = OrbitSpec(orbit_dim=1, jacobian_convention=JacobianConvention.NONE)
    result = classify_from_hessian(H_transverse, orbit=orbit)
    assert isinstance(result, RegularQuotient)
    assert result.datum.active_dim == 3
