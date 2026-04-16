"""Tests for kernel reduction (kernel_reduction.py, 0.3.3 Prop 9.1)."""
import numpy as np
import pytest

from nomogeo.kernel_reduction import (
    PolynomialJet,
    ReducedKernelActionDatum,
    reduce_kernel_action,
)
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    QuadraticKernelSeed,
)
from nomogeo.exceptions import InputValidationError


# ═══════════════════════════════════════════════════════════════════════
# PolynomialJet basics
# ═══════════════════════════════════════════════════════════════════════

def test_polynomial_jet_evaluate() -> None:
    """Evaluate a simple polynomial."""
    # Ψ(x₁, x₂) = 3x₁² + 2x₁x₂ + x₂⁴
    jet = PolynomialJet(
        dim=2,
        terms={(2, 0): 3.0, (1, 1): 2.0, (0, 4): 1.0},
    )
    # At (1, 1): 3 + 2 + 1 = 6
    assert jet.evaluate(np.array([1.0, 1.0])) == pytest.approx(6.0)
    # At (0, 0): 0
    assert jet.evaluate(np.array([0.0, 0.0])) == pytest.approx(0.0)
    # At (2, 0): 3*4 = 12
    assert jet.evaluate(np.array([2.0, 0.0])) == pytest.approx(12.0)


def test_polynomial_jet_gradient() -> None:
    """Gradient of Ψ(x₁, x₂) = x₁² + x₂² at (1, 2)."""
    jet = PolynomialJet(
        dim=2,
        terms={(2, 0): 1.0, (0, 2): 1.0},
    )
    grad = jet.gradient(np.array([1.0, 2.0]))
    assert grad[0] == pytest.approx(2.0)  # 2x₁
    assert grad[1] == pytest.approx(4.0)  # 2x₂


def test_polynomial_jet_hessian_at_origin() -> None:
    """Hessian of Ψ = 3x₁² + 2x₁x₂ + 5x₂² + x₁⁴ at origin."""
    jet = PolynomialJet(
        dim=2,
        terms={(2, 0): 3.0, (1, 1): 2.0, (0, 2): 5.0, (4, 0): 1.0},
    )
    H = jet.hessian_at_origin()
    assert H[0, 0] == pytest.approx(6.0)   # 2 * 3
    assert H[0, 1] == pytest.approx(2.0)
    assert H[1, 0] == pytest.approx(2.0)
    assert H[1, 1] == pytest.approx(10.0)  # 2 * 5


def test_polynomial_jet_restrict_to_subspace() -> None:
    """Restrict Ψ(x₁, x₂) = x₁² + x₂² to the x₁-axis."""
    jet = PolynomialJet(
        dim=2,
        terms={(2, 0): 1.0, (0, 2): 1.0},
    )
    # basis = [[1], [0]] — the x₁-axis
    basis = np.array([[1.0], [0.0]])
    restricted = jet.restrict_to_subspace(basis)
    assert restricted.dim == 1
    # Should be u² (only x₁² survives)
    assert restricted.evaluate(np.array([2.0])) == pytest.approx(4.0)


def test_polynomial_jet_restrict_to_diagonal() -> None:
    """Restrict Ψ(x₁, x₂) = x₁² + x₂² to the diagonal {(t,t)}."""
    jet = PolynomialJet(
        dim=2,
        terms={(2, 0): 1.0, (0, 2): 1.0},
    )
    # basis = [[1/√2], [1/√2]]
    s = 1.0 / np.sqrt(2.0)
    basis = np.array([[s], [s]])
    restricted = jet.restrict_to_subspace(basis)
    assert restricted.dim == 1
    # Ψ(t/√2, t/√2) = t²/2 + t²/2 = t²
    assert restricted.evaluate(np.array([1.0])) == pytest.approx(1.0)
    assert restricted.evaluate(np.array([3.0])) == pytest.approx(9.0)


# ═══════════════════════════════════════════════════════════════════════
# Kernel reduction: pure-quadratic-normal case
# ═══════════════════════════════════════════════════════════════════════

def _make_seed(kernel_dim, normal_dim, normal_eigenvalues, seed=0):
    """Build a QuadraticKernelSeed with standard basis alignment."""
    d = kernel_dim + normal_dim
    kernel_basis = np.eye(d, kernel_dim)     # first r cols
    normal_basis = np.eye(d, normal_dim, kernel_dim)  # next p cols
    # Fix: np.eye(d, normal_dim, kernel_dim) doesn't work as expected.
    # Use explicit construction.
    kernel_basis = np.zeros((d, kernel_dim))
    for i in range(kernel_dim):
        kernel_basis[i, i] = 1.0
    normal_basis = np.zeros((d, normal_dim))
    for i in range(normal_dim):
        normal_basis[kernel_dim + i, i] = 1.0

    eigs = np.array(normal_eigenvalues, dtype=float)
    log_prefactor = (
        0.5 * normal_dim * np.log(2.0 * np.pi)
        - 0.5 * np.sum(np.log(eigs))
    )

    return QuadraticKernelSeed(
        kernel_dim=kernel_dim,
        kernel_basis=kernel_basis,
        kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=kernel_dim),
        positive_normal_dim=normal_dim,
        positive_normal_basis=normal_basis,
        positive_normal_eigenvalues=eigs,
        log_normal_prefactor=log_prefactor,
    )


def test_pure_quartic_kernel_reduction() -> None:
    """Ψ(u, y) = u⁴ + y² (1-dim kernel, 1-dim normal).

    H_N = [[2]], kernel action Φ(u) = u⁴.
    No mixed cubic terms, so no quartic correction.
    """
    seed = _make_seed(kernel_dim=1, normal_dim=1, normal_eigenvalues=[2.0])
    jet = PolynomialJet(
        dim=2,
        terms={(4, 0): 1.0, (0, 2): 1.0},
    )
    result = reduce_kernel_action(seed, jet)
    assert isinstance(result, ReducedKernelActionDatum)
    assert result.kernel_dim == 1
    assert result.positive_normal_dim == 1

    # The reduced action should be Φ(u) = u⁴.
    phi = result.kernel_action
    assert phi.dim == 1
    assert phi.evaluate(np.array([2.0])) == pytest.approx(16.0)


def test_quartic_plus_sextic_kernel_reduction() -> None:
    """Ψ(u₁, u₂, y) = u₁⁴ + u₂⁶ + y² (2-dim kernel, 1-dim normal)."""
    seed = _make_seed(kernel_dim=2, normal_dim=1, normal_eigenvalues=[2.0])
    jet = PolynomialJet(
        dim=3,
        terms={(4, 0, 0): 1.0, (0, 6, 0): 1.0, (0, 0, 2): 1.0},
    )
    result = reduce_kernel_action(seed, jet)
    assert result.kernel_dim == 2

    phi = result.kernel_action
    # Φ(u₁, u₂) = u₁⁴ + u₂⁶
    assert phi.evaluate(np.array([1.0, 0.0])) == pytest.approx(1.0)
    assert phi.evaluate(np.array([0.0, 1.0])) == pytest.approx(1.0)
    assert phi.evaluate(np.array([1.0, 1.0])) == pytest.approx(2.0)


def test_kernel_reduction_with_cubic_correction() -> None:
    """Ψ(u, y) = u⁴ + y² + u²y (mixed cubic term).

    The cubic coupling V(u) = u² (coefficient of y in the cubic part).
    H_NN = 2, so H_NN^{-1} = 0.5.
    The quartic correction: ΔΦ₄ = -0.5 * 0.5 * u² * u² = -0.25 u⁴.
    So Φ(u) = u⁴ - 0.25 u⁴ = 0.75 u⁴.
    """
    seed = _make_seed(kernel_dim=1, normal_dim=1, normal_eigenvalues=[2.0])
    jet = PolynomialJet(
        dim=2,
        terms={
            (4, 0): 1.0,    # u⁴
            (0, 2): 1.0,    # y²
            (2, 1): 1.0,    # u²y — cubic mixed term
        },
        certified_order=4,
    )
    result = reduce_kernel_action(seed, jet)
    phi = result.kernel_action

    # Φ(u) should be 0.75 u⁴.
    assert phi.evaluate(np.array([1.0])) == pytest.approx(0.75)
    assert phi.evaluate(np.array([2.0])) == pytest.approx(0.75 * 16.0)


def test_kernel_reduction_preserves_log_prefactor() -> None:
    """The normal prefactor from the seed should pass through."""
    eigs = [1.0, 3.0]
    seed = _make_seed(kernel_dim=1, normal_dim=2, normal_eigenvalues=eigs)
    jet = PolynomialJet(
        dim=3,
        terms={(4, 0, 0): 1.0, (0, 2, 0): 0.5, (0, 0, 2): 1.5},
    )
    result = reduce_kernel_action(seed, jet)
    assert result.log_normal_prefactor == pytest.approx(seed.log_normal_prefactor)


# ═══════════════════════════════════════════════════════════════════════
# Numerical verification against direct integration
# ═══════════════════════════════════════════════════════════════════════

def test_reduction_vs_direct_integration_1d() -> None:
    """Compare reduced kernel integral against full 2D numerical integral.

    Ψ(u, y) = u⁴ + 2y² (no mixed terms).
    Full integral: ∫∫ e^{-n(u⁴ + 2y²)} du dy
                 = ∫ e^{-nu⁴} du · ∫ e^{-2ny²} dy
                 = ∫ e^{-nu⁴} du · √(π/(2n))

    Kernel integral: ∫ e^{-nΦ(u)} du = ∫ e^{-nu⁴} du

    Ratio: n^{-1/2} · √(π/2) (the normal prefactor).
    """
    from scipy import integrate

    n = 100.0
    # Full 2D integral.
    full, _ = integrate.dblquad(
        lambda y, u: np.exp(-n * (u**4 + 2 * y**2)),
        -3, 3, -3, 3,
    )
    # 1D kernel integral.
    kernel, _ = integrate.quad(
        lambda u: np.exp(-n * u**4), -3, 3,
    )
    # Normal prefactor.
    normal_factor = np.sqrt(np.pi / (2 * n))

    # full ≈ kernel * normal_factor
    assert full == pytest.approx(kernel * normal_factor, rel=1e-4)


# ═══════════════════════════════════════════════════════════════════════
# Input validation
# ═══════════════════════════════════════════════════════════════════════

def test_rejects_dimension_mismatch() -> None:
    seed = _make_seed(kernel_dim=1, normal_dim=1, normal_eigenvalues=[1.0])
    jet = PolynomialJet(dim=3, terms={(2, 0, 0): 1.0})
    with pytest.raises(InputValidationError, match="dim"):
        reduce_kernel_action(seed, jet)
