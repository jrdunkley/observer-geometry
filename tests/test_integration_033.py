"""End-to-end integration test: tower → regime → evidence (0.3.3 pipeline)."""
import numpy as np
import pytest

from nomogeo.affine import tower_affine_hidden_elimination
from nomogeo.evidence import dispatch_evidence, local_evidence
from nomogeo.regime import classify_regime
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    OrbitSpec,
    ReducedLocalDatum,
    RegimeKind,
    RegularEvidenceTemplate,
    RegularInterior,
    UnresolvedEvidenceRefusal,
    UnresolvedKernel,
)


def _random_spd(d: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(d, d))
    return raw.T @ raw + 0.5 * np.eye(d)


# ═══════════════════════════════════════════════════════════════════════
# Full pipeline: tower → datum → classify → evidence
# ═══════════════════════════════════════════════════════════════════════

def test_tower_to_interior_evidence_pipeline() -> None:
    """Construct a two-layer hidden model, tower-eliminate, classify, dispatch."""
    # A 6-dim hidden model that we eliminate in two 3-dim blocks.
    D = _random_spd(6, seed=200)
    J = np.array([0.1, -0.2, 0.3, 0.4, -0.5, 0.6])
    A = 1.0

    # Tower elimination.
    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("layer_1", [0, 1, 2]), ("layer_2", [0, 1, 2])],
    )

    # All hidden eliminated → final_precision is 0×0.
    assert tower.final_coupling.shape == (0,)
    assert tower.final_precision.shape == (0, 0)

    # Now suppose the visible Hessian after elimination is given by:
    visible_H = _random_spd(3, seed=201)  # 3-dim visible sector

    # Build datum → classify → evidence.
    datum = ReducedLocalDatum(
        active_dim=3,
        h_active=visible_H,
        cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=3),
    )
    classification = classify_regime(datum)
    assert isinstance(classification, RegularInterior)

    evidence = dispatch_evidence(classification)
    assert isinstance(evidence, RegularEvidenceTemplate)
    assert evidence.regime == RegimeKind.REGULAR_INTERIOR
    assert evidence.active_dim == 3

    # The total log evidence would be:
    # log Z ≈ -tower.final_action + evidence.log_quadratic_det_term
    # (combining hidden elimination with visible Laplace approximation).
    total_log_evidence = -tower.final_action + evidence.log_quadratic_det_term
    assert np.isfinite(total_log_evidence)


def test_tower_partial_to_kernel_regime() -> None:
    """Tower leaves residual hidden → build PSD datum → kernel regime."""
    D = _random_spd(4, seed=202)
    J = np.zeros(4)
    A = 0.0

    # Only eliminate 2 of 4 hidden variables.
    tower = tower_affine_hidden_elimination(
        A, J, D, stages=[("first_half", [0, 1])]
    )

    assert tower.final_precision.shape == (2, 2)

    # Suppose the visible Hessian is singular (PSD with 1-dim kernel).
    v = np.array([1.0, 0.5])
    visible_H = np.outer(v, v)  # rank 1, kernel dim 1

    datum = ReducedLocalDatum(
        active_dim=2,
        h_active=visible_H,
        cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=2),
    )
    classification = classify_regime(datum)
    assert isinstance(classification, UnresolvedKernel)
    assert classification.kernel.kernel_dim == 1

    evidence = dispatch_evidence(classification)
    assert isinstance(evidence, UnresolvedEvidenceRefusal)
    assert evidence.kernel_dim == 1


def test_tower_with_orbit_quotient_evidence() -> None:
    """Tower + quotient orbit → quotient evidence with volume data.

    Convention check: the visible sector is 4-dim ambient, with orbit_dim=1.
    After exact slice reduction, the transverse Hessian is 3×3.  The
    ReducedLocalDatum carries this 3×3 transverse Hessian, and the orbit
    volume multiplies the evidence externally.  The Gaussian dimension
    in the evidence template must be 3 (transverse), not 4 (ambient).
    """
    D = _random_spd(4, seed=203)
    J = np.ones(4) * 0.1
    A = 0.5

    tower = tower_affine_hidden_elimination(
        A, J, D, stages=[("hidden", [0, 1, 2, 3])]
    )

    # Simulate a 4-dim visible space with orbit_dim=1:
    # After exact slice reduction, the transverse Hessian is 3×3.
    transverse_H = _random_spd(3, seed=204)
    orbit = OrbitSpec(
        orbit_dim=1,
        jacobian_convention=JacobianConvention.SLICE_LEBESGUE,
        log_orbit_volume=np.log(2 * np.pi),
        log_slice_jacobian=0.0,
    )
    datum = ReducedLocalDatum(
        active_dim=3,  # transverse dimension, NOT ambient
        h_active=transverse_H,  # 3×3 transverse Hessian
        cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=3),
        orbit=orbit,
    )
    evidence = local_evidence(datum)
    assert isinstance(evidence, RegularEvidenceTemplate)
    assert evidence.regime == RegimeKind.REGULAR_QUOTIENT
    assert evidence.active_dim == 3  # transverse dim
    assert evidence.log_orbit_volume == pytest.approx(np.log(2 * np.pi))

    # Verify the Gaussian prefactor uses transverse dimension.
    expected_det_term = (
        0.5 * 3 * np.log(2.0 * np.pi)
        - 0.5 * np.linalg.slogdet(transverse_H)[1]
    )
    assert evidence.log_quadratic_det_term == pytest.approx(expected_det_term)


# ═══════════════════════════════════════════════════════════════════════
# Consistency: tower half-logdet + remaining half-logdet = full half-logdet
# ═══════════════════════════════════════════════════════════════════════

def test_tower_half_logdet_plus_residual_equals_full_half_logdet() -> None:
    """When tower eliminates some variables, tower.total_half_log_det +
    (1/2)*logdet(residual) should equal (1/2)*logdet(full D).

    This is exact for J=0 (no coupling shift), since the determinant
    of the Schur complement times the eliminated block's determinant
    equals the full matrix determinant."""
    D = _random_spd(6, seed=205)
    J = np.zeros(6)
    A = 0.0

    tower = tower_affine_hidden_elimination(
        A, J, D, stages=[("block_a", [0, 1, 2])]
    )

    # tower.total_half_log_det = (1/2) logdet(D_ee)
    # The residual is the Schur complement S = D_kk - D_ke D_ee^{-1} D_ek
    # and det(D) = det(D_ee) * det(S)
    # so (1/2)*logdet(D) = tower.total_half_log_det + (1/2)*logdet(S)
    full_half_logdet = 0.5 * np.linalg.slogdet(D)[1]
    residual_half_logdet = 0.5 * np.linalg.slogdet(tower.final_precision)[1]
    assert tower.total_half_log_det + residual_half_logdet == pytest.approx(full_half_logdet)


def test_full_pipeline_scalar_consistency() -> None:
    """For J=0, tower_action + evidence_det_term should match the
    direct Laplace approximation of a purely quadratic model.

    The accounting:
      tower.final_action = total_half_log_det  (since A=0, J=0)
      evidence.log_quadratic_det_term = (m/2)log(2π) - (1/2)logdet(S)
      hidden_gaussian_prefactor = (k/2)log(2π)  [for k eliminated dims]

    Combined: -tower.final_action + evidence_det_term + hidden_prefactor
            = -(1/2)logdet(D_ee) + (m/2)log(2π) - (1/2)logdet(S) + (k/2)log(2π)
            = (d/2)log(2π) - (1/2)logdet(D)
            = direct Laplace
    """
    d = 4
    H = _random_spd(d, seed=206)

    # Direct Laplace for A=0, full quadratic:
    direct = 0.5 * d * np.log(2.0 * np.pi) - 0.5 * np.linalg.slogdet(H)[1]

    # Split: eliminate 2 hidden, use 2×2 Schur complement as visible.
    J = np.zeros(d)
    A = 0.0
    tower = tower_affine_hidden_elimination(
        A, J, H, stages=[("hidden", [0, 1])]
    )

    visible_H = tower.final_precision  # 2×2 Schur complement
    datum = ReducedLocalDatum(
        active_dim=2,
        h_active=visible_H,
        cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=2),
    )
    evidence = local_evidence(datum)
    assert isinstance(evidence, RegularEvidenceTemplate)

    # The tower handles the determinant but not the (2π)^{k/2} prefactor
    # for the eliminated dimensions.  That is a separate accounting step.
    hidden_dim = 2
    hidden_gaussian_prefactor = 0.5 * hidden_dim * np.log(2.0 * np.pi)
    combined = -tower.final_action + evidence.log_quadratic_det_term + hidden_gaussian_prefactor
    assert combined == pytest.approx(direct)
