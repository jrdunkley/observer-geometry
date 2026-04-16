"""Tests for the typed evidence dispatcher (evidence.py, 0.3.3)."""
import numpy as np
import pytest

from nomogeo.evidence import dispatch_evidence, local_evidence, local_evidence_from_hessian
from nomogeo.regime import classify_regime, classify_from_hessian
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    OrbitSpec,
    QuadraticKernelSeed,
    ReducedLocalDatum,
    RegimeKind,
    RegularCone,
    RegularEvidenceTemplate,
    RegularInterior,
    RegularQuotient,
    RegularQuotientCone,
    UnresolvedEvidenceRefusal,
    UnresolvedKernel,
    UnsupportedConeReduction,
    UnsupportedEvidenceRefusal,
)


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


def _datum(H, cone=None, orbit=None):
    d = H.shape[0]
    if cone is None:
        cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d)
    return ReducedLocalDatum(active_dim=d, h_active=H, cone=cone, orbit=orbit)


# ═══════════════════════════════════════════════════════════════════════
# dispatch_evidence on each regime variant
# ═══════════════════════════════════════════════════════════════════════

def test_regular_interior_evidence() -> None:
    H = _spd(3, seed=100)
    cls = classify_regime(_datum(H))
    assert isinstance(cls, RegularInterior)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_INTERIOR
    assert ev.active_dim == 3

    # Check the log-det term: (m/2) log(2π) - (1/2) log det(H)
    expected = 0.5 * 3 * np.log(2.0 * np.pi) - 0.5 * np.linalg.slogdet(H)[1]
    assert ev.log_quadratic_det_term == pytest.approx(expected)
    assert ev.log_cone_mass is None
    assert ev.log_orbit_volume is None


def test_regular_cone_evidence() -> None:
    H = _spd(3, seed=101)
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)
    cls = classify_regime(_datum(H, cone=cone))
    assert isinstance(cls, RegularCone)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_CONE
    assert ev.log_cone_mass is None  # deferred in 0.3.3
    assert ev.cone_mass_known_exact is True
    assert ev.log_orbit_volume is None


def test_regular_quotient_evidence() -> None:
    """Quotient evidence with genuinely reduced transverse datum.
    Ambient 4-dim, orbit_dim=1, so transverse H is 3×3."""
    H_transverse = _spd(3, seed=102)  # 3×3 transverse Hessian
    orbit = OrbitSpec(
        orbit_dim=1,
        jacobian_convention=JacobianConvention.SLICE_LEBESGUE,
        log_orbit_volume=1.5,
        log_slice_jacobian=0.3,
    )
    cls = classify_regime(_datum(H_transverse, orbit=orbit))
    assert isinstance(cls, RegularQuotient)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_QUOTIENT
    assert ev.active_dim == 3  # transverse dim, not ambient
    assert ev.log_orbit_volume == pytest.approx(1.5)
    assert ev.log_slice_jacobian == pytest.approx(0.3)
    assert ev.log_cone_mass is None

    # Check the Gaussian dimension is the transverse dimension.
    # (m/2) log(2π) - (1/2) log det(H_⊥) where m = 3.
    expected = 0.5 * 3 * np.log(2.0 * np.pi) - 0.5 * np.linalg.slogdet(H_transverse)[1]
    assert ev.log_quadratic_det_term == pytest.approx(expected)


def test_regular_quotient_cone_evidence() -> None:
    """Quotient-cone evidence with genuinely reduced transverse datum.
    Ambient 5-dim, orbit_dim=2, so transverse H is 3×3."""
    H_transverse = _spd(3, seed=103)  # 3×3 transverse Hessian
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)  # cone on transverse space
    orbit = OrbitSpec(
        orbit_dim=2,
        jacobian_convention=JacobianConvention.SLICE_RIEMANNIAN,
        log_orbit_volume=2.0,
        log_slice_jacobian=-0.1,
    )
    cls = classify_regime(_datum(H_transverse, cone=cone, orbit=orbit))
    assert isinstance(cls, RegularQuotientCone)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_QUOTIENT_CONE
    assert ev.active_dim == 3  # transverse dim
    assert ev.log_cone_mass is None  # deferred
    assert ev.log_orbit_volume == pytest.approx(2.0)
    assert ev.log_slice_jacobian == pytest.approx(-0.1)


def test_unresolved_kernel_evidence_refusal() -> None:
    H = _psd_with_kernel(4, kernel_dim=1, seed=104)
    cls = classify_regime(_datum(H))
    assert isinstance(cls, UnresolvedKernel)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, UnresolvedEvidenceRefusal)
    assert ev.kernel_dim == 1
    assert "kernel" in ev.reason.lower()


def test_unsupported_cone_evidence_refusal() -> None:
    H = _spd(3, seed=105)
    # Cone with implicit equality: x_1 >= 0 AND -x_1 >= 0 → x_1 = 0
    A = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    cone = ConeSpec(kind=ConeKind.POLYHEDRAL, ambient_dim=3, halfspace_normals=A)
    cls = classify_regime(_datum(H, cone=cone))
    assert isinstance(cls, UnsupportedConeReduction)
    ev = dispatch_evidence(cls)
    assert isinstance(ev, UnsupportedEvidenceRefusal)
    assert ev.cone_kind == ConeKind.POLYHEDRAL


# ═══════════════════════════════════════════════════════════════════════
# Convenience entry points
# ═══════════════════════════════════════════════════════════════════════

def test_local_evidence_composites_classify_dispatch() -> None:
    H = _spd(3, seed=106)
    datum = _datum(H)
    ev = local_evidence(datum)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_INTERIOR


def test_local_evidence_from_hessian_composites() -> None:
    H = _spd(3, seed=107)
    ev = local_evidence_from_hessian(H)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_INTERIOR


def test_local_evidence_from_hessian_with_cone_and_orbit() -> None:
    """Convenience wrapper with cone + orbit: H must be transverse.
    Ambient 4-dim, orbit_dim=1 → transverse H is 3×3."""
    H_transverse = _spd(3, seed=108)  # 3×3 transverse Hessian
    cone = ConeSpec(kind=ConeKind.ORTHANT, ambient_dim=3)  # cone on transverse space
    orbit = OrbitSpec(orbit_dim=1, jacobian_convention=JacobianConvention.NONE)
    ev = local_evidence_from_hessian(H_transverse, cone=cone, orbit=orbit)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.regime == RegimeKind.REGULAR_QUOTIENT_CONE
    assert ev.active_dim == 3


# ═══════════════════════════════════════════════════════════════════════
# Evidence numerical correctness
# ═══════════════════════════════════════════════════════════════════════

def test_interior_evidence_matches_laplace_formula() -> None:
    """For a diagonal SPD matrix, check against the explicit Laplace formula."""
    d = 5
    diag_vals = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    H = np.diag(diag_vals)

    ev = local_evidence_from_hessian(H)
    assert isinstance(ev, RegularEvidenceTemplate)

    # (m/2) log(2π) - (1/2) Σ log(λ_i)
    expected = 0.5 * d * np.log(2.0 * np.pi) - 0.5 * np.sum(np.log(diag_vals))
    assert ev.log_quadratic_det_term == pytest.approx(expected)


def test_quotient_without_orbit_data_has_none_fields() -> None:
    """If OrbitSpec has no volume/Jacobian, evidence should carry None.
    Transverse H is 3×3 (ambient 4-dim, orbit_dim=1)."""
    H_transverse = _spd(3, seed=109)  # 3×3 transverse Hessian
    orbit = OrbitSpec(orbit_dim=1, jacobian_convention=JacobianConvention.NONE)
    ev = local_evidence_from_hessian(H_transverse, orbit=orbit)
    assert isinstance(ev, RegularEvidenceTemplate)
    assert ev.active_dim == 3
    assert ev.log_orbit_volume is None
    assert ev.log_slice_jacobian is None


def test_unknown_classification_type_raises() -> None:
    """dispatch_evidence should reject unrecognized classification types."""
    with pytest.raises(TypeError, match="Unknown"):
        dispatch_evidence("not_a_classification")  # type: ignore
