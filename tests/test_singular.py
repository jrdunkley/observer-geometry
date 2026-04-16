"""Tests for the certified singular classifier (singular.py, 0.3.3)."""
import numpy as np
import pytest
from math import lgamma

from nomogeo.singular import classify_singular_kernel
from nomogeo.singular_types import (
    BranchChannelTemplate,
    CriticalBranchTemplate,
    SingularRegimeKind,
    UnresolvedKernelJetRefusal,
    WeightedHomogeneousTemplate,
)
from nomogeo.kernel_reduction import (
    PolynomialJet,
    ReducedKernelActionDatum,
)
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
)


def _make_datum(kernel_dim, jet, positive_normal_dim=0, log_normal_prefactor=0.0):
    """Build a ReducedKernelActionDatum for testing."""
    return ReducedKernelActionDatum(
        kernel_dim=kernel_dim,
        kernel_action=jet,
        kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=kernel_dim),
        positive_normal_dim=positive_normal_dim,
        log_normal_prefactor=log_normal_prefactor,
        total_dim=kernel_dim + positive_normal_dim,
    )


# ═══════════════════════════════════════════════════════════════════════
# Channel 1: Weighted-homogeneous — 1D cases
# ═══════════════════════════════════════════════════════════════════════

def test_quartic_kernel_1d() -> None:
    """Φ(u) = u⁴ → weights = (1/4,), exponent = 1/4."""
    jet = PolynomialJet(dim=1, terms={(4,): 1.0})
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.kind == SingularRegimeKind.WEIGHTED_HOMOGENEOUS
    assert result.weights == pytest.approx((0.25,))
    assert result.total_weight == pytest.approx(0.25)
    assert result.leading_monomial_degrees == (4,)


def test_sextic_kernel_1d() -> None:
    """Φ(u) = u⁶ → weights = (1/6,), exponent = 1/6."""
    jet = PolynomialJet(dim=1, terms={(6,): 1.0}, certified_order=6)
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.weights == pytest.approx((1.0 / 6.0,))
    assert result.total_weight == pytest.approx(1.0 / 6.0)


def test_quartic_with_positive_coefficient() -> None:
    """Φ(u) = 3u⁴ → same structure, different coefficient."""
    jet = PolynomialJet(dim=1, terms={(4,): 3.0})
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.weights == pytest.approx((0.25,))


# ═══════════════════════════════════════════════════════════════════════
# Channel 1: Weighted-homogeneous — 2D cases
# ═══════════════════════════════════════════════════════════════════════

def test_anisotropic_quartic_sextic() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁶ → ω = (1/4, 1/6), |ω| = 5/12."""
    jet = PolynomialJet(dim=2, terms={(4, 0): 1.0, (0, 6): 1.0}, certified_order=6)
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.weights == pytest.approx((0.25, 1.0 / 6.0))
    assert result.total_weight == pytest.approx(5.0 / 12.0)
    assert result.leading_monomial_degrees == (4, 6)


def test_product_quartic() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁴ → ω = (1/4, 1/4), |ω| = 1/2."""
    jet = PolynomialJet(dim=2, terms={(4, 0): 1.0, (0, 4): 1.0})
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.total_weight == pytest.approx(0.5)


def test_wh_rejects_same_weight_mixed_principal_part() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁴ + 0.5·u₁²u₂² → weighted degree of
    cross term = 2*(1/4) + 2*(1/4) = 1.

    The cross term sits on the Newton boundary (weighted degree = 1),
    which makes it part of the principal germ.  The separable product
    integral formula is NOT licensed for non-separable principal parts.
    The classifier must refuse."""
    jet = PolynomialJet(dim=2, terms={(4, 0): 1.0, (0, 4): 1.0, (2, 2): 0.5})
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert not isinstance(result, WeightedHomogeneousTemplate)


def test_wh_with_strictly_subordinate_cross_term() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁶ + u₁²u₂² → weighted degree of
    cross term = 2*(1/4) + 2*(1/6) = 5/6 + 1/3 = 10/12 + 4/12... wait:
    Actually 2/4 + 2/6 = 1/2 + 1/3 = 5/6 < 1.  Subordination fails.

    Use instead: u₁⁴ + u₂⁴ + u₁u₂³ → cross term weighted degree
    = 1/4 + 3/4 = 1.0.  Still on Newton boundary, refuse.

    For a truly subordinate cross term, need weighted degree > 1:
    u₁⁴ + u₂⁴ + u₁²u₂⁴ → weighted degree = 2/4 + 4/4 = 3/2 > 1.
    """
    jet = PolynomialJet(dim=2, terms={
        (4, 0): 1.0, (0, 4): 1.0, (2, 4): 0.5,
    })
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)


# ═══════════════════════════════════════════════════════════════════════
# Channel 1: principal integral correctness
# ═══════════════════════════════════════════════════════════════════════

def test_quartic_principal_integral() -> None:
    """For Φ(u) = u⁴ on R, the principal integral is
    ∫ e^{-|u|⁴} du = 2·Γ(1/4) / 4."""
    jet = PolynomialJet(dim=1, terms={(4,): 1.0})
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.log_principal_integral is not None

    # Expected: 2 · Γ(1/4) / 4
    from math import gamma
    expected = 2.0 * gamma(0.25) / 4.0
    assert np.exp(result.log_principal_integral) == pytest.approx(expected, rel=1e-6)


def test_sextic_principal_integral() -> None:
    """For Φ(u) = u⁶ on R, the integral is 2·Γ(1/6) / 6."""
    jet = PolynomialJet(dim=1, terms={(6,): 1.0}, certified_order=6)
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.log_principal_integral is not None

    from math import gamma
    expected = 2.0 * gamma(1.0 / 6.0) / 6.0
    assert np.exp(result.log_principal_integral) == pytest.approx(expected, rel=1e-6)


def test_2d_product_integral_factors() -> None:
    """For Φ(u₁, u₂) = u₁⁴ + u₂⁶, the integral factors into 1D integrals."""
    jet = PolynomialJet(dim=2, terms={(4, 0): 1.0, (0, 6): 1.0}, certified_order=6)
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    assert result.log_principal_integral is not None

    from math import gamma
    I1 = 2.0 * gamma(0.25) / 4.0
    I2 = 2.0 * gamma(1.0 / 6.0) / 6.0
    expected = I1 * I2
    assert np.exp(result.log_principal_integral) == pytest.approx(expected, rel=1e-6)


# ═══════════════════════════════════════════════════════════════════════
# Channel 1: Weighted-homogeneous scaling verification
# ═══════════════════════════════════════════════════════════════════════

def test_quartic_scaling_law_numerical() -> None:
    """Numerically verify n^{-1/4} scaling for Φ(u) = u⁴."""
    from scipy import integrate

    ns = [100, 400, 1600]
    integrals = []
    for n_val in ns:
        val, _ = integrate.quad(lambda u: np.exp(-n_val * u**4), -5, 5)
        integrals.append(val)

    # Check that I(n) ~ n^{-1/4}: log(I) ~ -1/4 * log(n).
    log_n = np.log(ns)
    log_I = np.log(integrals)
    # Fit slope.
    slope = (log_I[-1] - log_I[0]) / (log_n[-1] - log_n[0])
    assert slope == pytest.approx(-0.25, abs=0.02)


# ═══════════════════════════════════════════════════════════════════════
# Channel 2: Branch-restiffening — critical quadratic branch
# ═══════════════════════════════════════════════════════════════════════

def test_critical_quadratic_branch() -> None:
    """Φ(u, v) = u⁴ + u²v² + v⁴.

    Branch order c = 2.
    (1/2)∂²ᵤΦ(0,v) = 0 + v² → H₀ = 1.0, stiffness_ord = 2.
    Φ(0,v) = v⁴ → d = 4, B₀ = 1.0.
    ∂ᵤ⁴Φ(0,0) = 4! · 1 = 24 → C₀ = 1.0, a = 4.
    κ = (4(4-2) - 2·4) / (2·4) = (8 - 8) / 8 = 0.

    But κ = 0 fails the positivity condition.  Let's use a different example.
    """
    pass  # See next test


def test_critical_quadratic_branch_positive_kappa() -> None:
    """Φ(u, v) = u⁶ + u²v² + v⁴.

    Stiffness: (1/2)∂²ᵤΦ(0,v) at u=0 is the coefficient of u² as f(v):
        The u² terms: u²v² → coefficient of u² is v².
        So (1/2)∂²ᵤΦ(0,v) = v².  → H₀ = 1.0, stiffness_ord = 2.
    Note: the u⁶ contributes a u⁶ term, not a u² term.

    Φ(0,v) = v⁴ → d = 4, B₀ = 1.0.
    Higher transverse at v=0: u⁶ → a = 6, C₀ = 1.0.
    κ = (4(6-2) - 2·6) / (2·6) = (16 - 12) / 12 = 4/12 = 1/3.
    """
    jet = PolynomialJet(
        dim=2,
        terms={(6, 0): 1.0, (2, 2): 1.0, (0, 4): 1.0},
        certified_order=6,
    )
    datum = _make_datum(2, jet, positive_normal_dim=3, log_normal_prefactor=1.5)
    result = classify_singular_kernel(datum, branch_straightened=True)
    assert isinstance(result, CriticalBranchTemplate)
    assert result.kind == SingularRegimeKind.CRITICAL_QUADRATIC_BRANCH
    assert result.kappa == pytest.approx(1.0 / 3.0)
    assert result.branch_datum.branch_order == 2
    assert result.branch_datum.stiffness_vanishing_order == 2
    assert result.branch_datum.branch_action_degree == 4
    assert result.branch_datum.higher_transverse_degree == 6
    assert result.branch_datum.kappa == pytest.approx(1.0 / 3.0)
    assert result.branch_datum.stiffness_leading_coeff == pytest.approx(1.0)
    assert result.branch_datum.branch_leading_coeff == pytest.approx(1.0)
    assert result.branch_datum.higher_transverse_coeff == pytest.approx(1.0)
    # Normal prefactor passes through.
    assert result.log_normal_prefactor == pytest.approx(1.5)
    assert result.positive_normal_dim == 3


def test_critical_branch_log_coefficient() -> None:
    """Verify the log coefficient for the critical quadratic branch.

    log_coeff = log(√π / d · H₀^{-1/2} · κ)
    For d=4, H₀=1, κ=1/3:
      log(√π / 4 · 1 · 1/3) = log(√π / 12)
    """
    jet = PolynomialJet(
        dim=2,
        terms={(6, 0): 1.0, (2, 2): 1.0, (0, 4): 1.0},
        certified_order=6,
    )
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum, branch_straightened=True)
    assert isinstance(result, CriticalBranchTemplate)

    expected = np.log(np.sqrt(np.pi) / 4.0 * 1.0 * (1.0 / 3.0))
    # = log(√π / 12)
    assert result.log_coefficient == pytest.approx(expected)


# ═══════════════════════════════════════════════════════════════════════
# Channel 2: non-critical branch channel
# ═══════════════════════════════════════════════════════════════════════

def test_non_critical_branch_channel() -> None:
    """Φ(u, v) = u⁶ + u²v³ + v⁶ (stiffness vanishing order = 3, not 2).

    WH fails: pure powers u⁶, v⁶ give weights (1/6, 1/6), but the
    cross term u²v³ has weighted degree 2/6 + 3/6 = 5/6 < 1, so
    subordination fails.

    Branch analysis with the corrected branch-order detection:
      c = 2  (u² coefficient as polynomial in v is v³, not identically zero).
      H_v_coeffs = {3: 1.0}.
      h_at_zero = 0 → stiffness vanishes at origin.
      stiffness_ord = 3, H₀ = 1.0.
      Since stiffness_ord ≠ 2 → general BranchChannelTemplate.
    """
    jet = PolynomialJet(
        dim=2,
        terms={(6, 0): 1.0, (2, 3): 1.0, (0, 6): 1.0},
        certified_order=6,
    )
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum, branch_straightened=True)
    # This should be a general BranchChannelTemplate (not critical).
    assert isinstance(result, BranchChannelTemplate)
    assert result.kind == SingularRegimeKind.BRANCH_RESTIFFENING
    assert result.branch_datum.stiffness_vanishing_order == 3
    assert result.branch_datum.branch_order == 2
    assert result.branch_datum.branch_channel_exponent == pytest.approx(1.5)


# ═══════════════════════════════════════════════════════════════════════
# Refusals
# ═══════════════════════════════════════════════════════════════════════

def test_refusal_no_pure_powers() -> None:
    """Φ(u₁, u₂) = u₁²u₂² (no pure power in either variable)."""
    jet = PolynomialJet(dim=2, terms={(2, 2): 1.0})
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, UnresolvedKernelJetRefusal)
    assert result.kind == SingularRegimeKind.UNRESOLVED_KERNEL_JET


def test_refusal_negative_coefficient() -> None:
    """Φ(u) = -u⁴ — negative coefficient, not positive."""
    jet = PolynomialJet(dim=1, terms={(4,): -1.0})
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, UnresolvedKernelJetRefusal)


def test_refusal_3d_kernel_no_branch() -> None:
    """3D kernel germ without separable principal part → refusal."""
    jet = PolynomialJet(
        dim=3,
        terms={(2, 2, 0): 1.0, (0, 2, 2): 1.0},  # only cross terms
    )
    datum = _make_datum(3, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, UnresolvedKernelJetRefusal)
    assert "branch" in result.reason or "pure-power" in result.reason


def test_refusal_subordination_fails() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁴ + u₁u₂².

    The cross term u₁u₂² has weighted degree 1*(1/4) + 2*(1/4) = 3/4 < 1.
    Subordination fails → cannot certify as weighted-homogeneous.
    Falls through to branch channel attempt (2D), which may also fail.
    """
    jet = PolynomialJet(
        dim=2,
        terms={(4, 0): 1.0, (0, 4): 1.0, (1, 2): 1.0},
    )
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    # Should not be WeightedHomogeneousTemplate.
    assert not isinstance(result, WeightedHomogeneousTemplate)


def test_refusal_restricted_cone() -> None:
    """Weighted-homogeneous on non-full-space cone → refusal
    (not yet supported for certification)."""
    jet = PolynomialJet(dim=1, terms={(4,): 1.0})
    datum = ReducedKernelActionDatum(
        kernel_dim=1,
        kernel_action=jet,
        kernel_cone=ConeSpec(
            kind=ConeKind.POLYHEDRAL,
            ambient_dim=1,
            halfspace_normals=np.array([[1.0]]),
        ),
        positive_normal_dim=0,
        log_normal_prefactor=0.0,
        total_dim=1,
    )
    result = classify_singular_kernel(datum)
    assert isinstance(result, UnresolvedKernelJetRefusal)


# ═══════════════════════════════════════════════════════════════════════
# End-to-end: asymptotic exponent verification
# ═══════════════════════════════════════════════════════════════════════

def test_quartic_exponent_from_classifier_matches_numerical() -> None:
    """The classifier says exponent = 1/4.  Numerical integration confirms."""
    from scipy import integrate

    jet = PolynomialJet(dim=1, terms={(4,): 1.0})
    datum = _make_datum(1, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    predicted_exponent = result.total_weight  # 0.25

    # Verify numerically.
    ns = [100, 1000, 10000]
    log_n = np.log(ns)
    log_I = []
    for n_val in ns:
        val, _ = integrate.quad(lambda u: np.exp(-n_val * u**4), -5, 5)
        log_I.append(np.log(val))

    slope = np.polyfit(log_n, log_I, 1)[0]
    assert slope == pytest.approx(-predicted_exponent, abs=0.02)


def test_2d_anisotropic_exponent_numerical() -> None:
    """Φ(u₁, u₂) = u₁⁴ + u₂⁶.  Predicted exponent = 1/4 + 1/6 = 5/12."""
    from scipy import integrate

    jet = PolynomialJet(dim=2, terms={(4, 0): 1.0, (0, 6): 1.0}, certified_order=6)
    datum = _make_datum(2, jet)
    result = classify_singular_kernel(datum)
    assert isinstance(result, WeightedHomogeneousTemplate)
    predicted = result.total_weight  # 5/12

    ns = [100, 400, 1600]
    log_n = np.log(ns)
    log_I = []
    for n_val in ns:
        val, _ = integrate.dblquad(
            lambda u2, u1: np.exp(-n_val * (u1**4 + u2**6)),
            -3, 3, -3, 3,
        )
        log_I.append(np.log(val))

    slope = np.polyfit(log_n, log_I, 1)[0]
    assert slope == pytest.approx(-predicted, abs=0.03)
