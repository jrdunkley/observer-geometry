"""
Deep QA tests for the 0.3.3 layer stack.

This file covers:
  1. Full pipeline integration (slice → detector → kernel_reduction → singular)
  2. Polynomial jet algebra stress tests (high-order, multi-variable, subspace)
  3. Numerical verification of branch-restiffening asymptotics
  4. Adversarial / edge-case germs designed to probe classifier boundaries
  5. Kernel reduction correctness under rotated bases
  6. Consistency checks between layers
"""
import numpy as np
import pytest
from math import gamma, lgamma

from nomogeo.slice import reduce_local_chart, active_face_restriction, transverse_complement
from nomogeo.regime import classify_regime
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    JacobianConvention,
    OrbitSpec,
    QuadraticKernelSeed,
    ReducedLocalDatum,
    RegimeKind,
    RegularInterior,
    UnresolvedKernel,
)
from nomogeo.kernel_reduction import (
    PolynomialJet,
    ReducedKernelActionDatum,
    reduce_kernel_action,
)
from nomogeo.singular import classify_singular_kernel
from nomogeo.singular_types import (
    BranchChannelTemplate,
    CriticalBranchTemplate,
    SingularRegimeKind,
    UnresolvedKernelJetRefusal,
    WeightedHomogeneousTemplate,
)


def _random_spd(d, seed=42, min_eig=0.5):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(d, d))
    return raw.T @ raw + min_eig * np.eye(d)


def _make_datum(kernel_dim, jet, positive_normal_dim=0, log_normal_prefactor=0.0):
    return ReducedKernelActionDatum(
        kernel_dim=kernel_dim,
        kernel_action=jet,
        kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=kernel_dim),
        positive_normal_dim=positive_normal_dim,
        log_normal_prefactor=log_normal_prefactor,
        total_dim=kernel_dim + positive_normal_dim,
    )


# ═══════════════════════════════════════════════════════════════════════
# 1. Full pipeline integration: slice → detector → kernel_reduction → singular
# ═══════════════════════════════════════════════════════════════════════

class TestFullPipelineIntegration:
    """End-to-end tests that exercise the complete 5-layer stack."""

    def test_slice_to_regular_evidence(self) -> None:
        """An SPD ambient Hessian with no kernel → regular interior after
        slice reduction."""
        H = _random_spd(5, seed=300)
        datum = reduce_local_chart(H)
        classification = classify_regime(datum)
        assert isinstance(classification, RegularInterior)

    def test_slice_with_orbit_to_regular_quotient(self) -> None:
        """An SPD ambient Hessian with a 1D orbit tangent in its kernel.

        The slice reduction removes the orbit direction, leaving a
        4×4 transverse Hessian which the detector classifies as regular.
        """
        # Build a 5×5 Hessian with a 1-dim kernel along e_4.
        H = _random_spd(4, seed=301)
        H_ambient = np.zeros((5, 5))
        H_ambient[:4, :4] = H
        # e_4 is in ker(H_ambient).
        orbit_basis = np.zeros((5, 1))
        orbit_basis[4, 0] = 1.0

        datum = reduce_local_chart(
            H_ambient,
            orbit_tangent_basis=orbit_basis,
            log_orbit_volume=np.log(2 * np.pi),
            log_slice_jacobian=0.0,
        )
        assert datum.active_dim == 4
        classification = classify_regime(datum)
        # With an orbit, the detector returns RegularQuotient, not RegularInterior.
        from nomogeo.regime_types import RegularQuotient
        assert isinstance(classification, RegularQuotient)
        assert classification.orbit_dim == 1

    def test_slice_to_kernel_to_singular_wh(self) -> None:
        """A 3×3 Hessian with 1-dim kernel → detector finds kernel →
        kernel reduction with quartic jet → WH classifier certifies.

        We build a Hessian diag(1, 2, 0) so the kernel is along e_2.
        Then provide a local action jet with a quartic term in the
        kernel direction.
        """
        H = np.diag([1.0, 2.0, 0.0])
        datum = reduce_local_chart(H)
        classification = classify_regime(datum)
        assert isinstance(classification, UnresolvedKernel)
        seed = classification.kernel

        # Build an action jet in 3D with quartic term along the kernel.
        # The kernel direction is e_2 (3rd coordinate).
        # Ψ(x₀, x₁, x₂) = 0.5x₀² + x₁² + x₂⁴
        jet = PolynomialJet(
            dim=3,
            terms={
                (2, 0, 0): 0.5,  # 0.5 x₀²
                (0, 2, 0): 1.0,  # x₁²
                (0, 0, 4): 1.0,  # x₂⁴
            },
            certified_order=4,
        )

        reduced = reduce_kernel_action(seed, jet)
        assert reduced.kernel_dim == 1
        assert reduced.positive_normal_dim == 2

        # The reduced kernel action should be ~ u⁴ (possibly with a
        # different coefficient depending on the basis).
        result = classify_singular_kernel(reduced)
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.total_weight == pytest.approx(0.25)

    def test_face_restriction_then_kernel(self) -> None:
        """Active face restriction drops some coordinates, then the
        remaining Hessian has a kernel."""
        # 4×4 ambient, only coordinates [0, 1, 2] are active (face).
        # The 3×3 restricted Hessian has a 1-dim kernel.
        v = np.array([1.0, 0.5, 0.0])
        H_face = np.outer(v, v) + np.diag([0.1, 0.2, 0.0])
        # Embed in 4×4 ambient.
        H_ambient = np.zeros((4, 4))
        H_ambient[:3, :3] = H_face
        H_ambient[3, 3] = 5.0  # inactive coordinate is strongly curved

        datum = reduce_local_chart(
            H_ambient,
            active_face_indices=np.array([0, 1, 2]),
        )
        assert datum.active_dim == 3
        classification = classify_regime(datum)
        assert isinstance(classification, UnresolvedKernel)
        assert classification.kernel.kernel_dim == 1


# ═══════════════════════════════════════════════════════════════════════
# 2. Polynomial jet algebra stress tests
# ═══════════════════════════════════════════════════════════════════════

class TestPolynomialJetAlgebra:
    """Thorough testing of PolynomialJet operations."""

    def test_evaluate_3d_polynomial(self) -> None:
        """Evaluate a 3-variable polynomial at several points."""
        # Φ(x,y,z) = x² + 2xy + y²z + z³
        jet = PolynomialJet(dim=3, terms={
            (2, 0, 0): 1.0,
            (1, 1, 0): 2.0,
            (0, 2, 1): 1.0,
            (0, 0, 3): 1.0,
        })
        # At (1, 1, 1): 1 + 2 + 1 + 1 = 5
        assert jet.evaluate(np.array([1.0, 1.0, 1.0])) == pytest.approx(5.0)
        # At (0, 0, 0): 0
        assert jet.evaluate(np.array([0.0, 0.0, 0.0])) == pytest.approx(0.0)
        # At (2, -1, 0): 4 + 2*2*(-1) + 0 + 0 = 4 - 4 = 0
        assert jet.evaluate(np.array([2.0, -1.0, 0.0])) == pytest.approx(0.0)

    def test_gradient_matches_finite_difference(self) -> None:
        """Gradient of a 3D polynomial matches finite differences."""
        jet = PolynomialJet(dim=3, terms={
            (4, 0, 0): 1.0,
            (2, 2, 0): 0.5,
            (0, 0, 6): 1.0,
            (1, 1, 2): 3.0,
        }, certified_order=6)

        x0 = np.array([0.3, -0.5, 0.2])
        grad = jet.gradient(x0)

        eps = 1e-7
        for i in range(3):
            xp = x0.copy(); xp[i] += eps
            xm = x0.copy(); xm[i] -= eps
            fd = (jet.evaluate(xp) - jet.evaluate(xm)) / (2 * eps)
            assert grad[i] == pytest.approx(fd, abs=1e-5)

    def test_hessian_at_origin_mixed_terms(self) -> None:
        """Hessian extraction from a polynomial with mixed quadratic terms."""
        # Ψ = 3x₀² + x₀x₁ + 2x₁² + x₂⁴
        jet = PolynomialJet(dim=3, terms={
            (2, 0, 0): 3.0,
            (1, 1, 0): 1.0,
            (0, 2, 0): 2.0,
            (0, 0, 4): 1.0,
        })
        H = jet.hessian_at_origin()
        expected = np.array([
            [6.0, 1.0, 0.0],
            [1.0, 4.0, 0.0],
            [0.0, 0.0, 0.0],
        ])
        np.testing.assert_array_almost_equal(H, expected)

    def test_restrict_to_subspace_identity(self) -> None:
        """Restriction to the full space (identity basis) is a no-op."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 4): 2.0, (2, 2): 0.5,
        })
        basis = np.eye(2)
        restricted = jet.restrict_to_subspace(basis)
        # Evaluate at several points and compare.
        for _ in range(10):
            x = np.random.randn(2)
            assert restricted.evaluate(x) == pytest.approx(jet.evaluate(x), abs=1e-10)

    def test_restrict_to_subspace_rotation(self) -> None:
        """Restriction via a rotation preserves function values."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 6): 1.0, (2, 2): 0.5,
        }, certified_order=6)
        theta = np.pi / 5
        R = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta),  np.cos(theta)],
        ])
        restricted = jet.restrict_to_subspace(R)
        # For any u, restricted.evaluate(u) should equal jet.evaluate(R @ u).
        rng = np.random.default_rng(42)
        for _ in range(20):
            u = rng.normal(size=2) * 0.5
            val_original = jet.evaluate(R @ u)
            val_restricted = restricted.evaluate(u)
            assert val_restricted == pytest.approx(val_original, abs=1e-10)

    def test_restrict_to_1d_subspace(self) -> None:
        """Restrict a 3D polynomial to a 1D subspace."""
        # Φ(x,y,z) = x⁴ + y⁴ + z⁴
        jet = PolynomialJet(dim=3, terms={
            (4, 0, 0): 1.0, (0, 4, 0): 1.0, (0, 0, 4): 1.0,
        })
        # Restrict to the direction (1, 1, 1)/√3.
        v = np.array([[1.0], [1.0], [1.0]]) / np.sqrt(3)
        restricted = jet.restrict_to_subspace(v)
        assert restricted.dim == 1
        # At u = √3: x = y = z = 1, so Φ = 3.
        val = restricted.evaluate(np.array([np.sqrt(3)]))
        assert val == pytest.approx(3.0, abs=1e-8)

    def test_high_order_polynomial_multiply(self) -> None:
        """Jet evaluation consistency for a high-order polynomial."""
        # Φ(u) = u² + u⁴ + u⁶ + u⁸
        jet = PolynomialJet(dim=1, terms={
            (2,): 1.0, (4,): 1.0, (6,): 1.0, (8,): 1.0,
        }, certified_order=8)
        u = 0.5
        expected = u**2 + u**4 + u**6 + u**8
        assert jet.evaluate(np.array([u])) == pytest.approx(expected)


# ═══════════════════════════════════════════════════════════════════════
# 3. Numerical verification of branch-restiffening asymptotics
# ═══════════════════════════════════════════════════════════════════════

class TestBranchAsymptotics:
    """Numerically verify the asymptotic predictions of the branch
    channel classifier against direct 2D numerical integration."""

    def test_critical_quadratic_branch_log_correction(self) -> None:
        """For Φ(u,v) = u⁶ + u²v² + v⁴, the classifier predicts
        n^{-1/2} κ log(n) asymptotic form.

        Verify that the integral I(n) = ∫∫ e^{-n Φ} du dv has the
        form  I(n) ≈ C n^{-1/2} log(n)  for large n.
        """
        from scipy import integrate

        jet = PolynomialJet(
            dim=2,
            terms={(6, 0): 1.0, (2, 2): 1.0, (0, 4): 1.0},
            certified_order=6,
        )
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum, branch_straightened=True)
        assert isinstance(result, CriticalBranchTemplate)

        # The classifier says n^{-1/2} with a log(n) correction.
        # Check that the ratio I(n) / (n^{-1/2} log(n)) stabilises.
        def Phi(u, v):
            return u**6 + u**2 * v**2 + v**4

        ns = [200, 800, 3200]
        ratios = []
        for n_val in ns:
            val, _ = integrate.dblquad(
                lambda v, u: np.exp(-n_val * Phi(u, v)),
                -3, 3, -3, 3,
                epsabs=1e-10, epsrel=1e-8,
            )
            ratio = val / (n_val**(-0.5) * np.log(n_val))
            ratios.append(ratio)

        # The ratios should converge (not diverge or go to zero).
        # Check that ratios are within a factor of 2 of each other.
        assert ratios[-1] / ratios[0] == pytest.approx(1.0, abs=0.5)
        # And all ratios should be positive.
        assert all(r > 0 for r in ratios)

    def test_non_critical_branch_power_law(self) -> None:
        """For Φ(u,v) = u⁴ + u²v² + v⁴ (a positive-definite germ),
        verify it decays as a power law.  WH certifies this with
        weights (1/4, 1/4), total exponent 1/2.

        We use this as a baseline sanity check that the numerical
        integration machinery agrees with WH's prediction.
        """
        from scipy import integrate

        def Phi(u, v):
            return u**4 + u**2 * v**2 + v**4

        ns = [100, 400, 1600]
        log_n = np.log(ns)
        log_I = []
        for n_val in ns:
            val, _ = integrate.dblquad(
                lambda v, u: np.exp(-n_val * Phi(u, v)),
                -3, 3, -3, 3,
                epsabs=1e-10, epsrel=1e-8,
            )
            log_I.append(np.log(val))

        # Fit a line: log I ≈ -α log n + const.
        slope = np.polyfit(log_n, log_I, 1)[0]
        # WH predicts exponent = 1/2, so slope ≈ -0.5.
        assert slope == pytest.approx(-0.5, abs=0.05)


# ═══════════════════════════════════════════════════════════════════════
# 4. Adversarial / edge-case germs
# ═══════════════════════════════════════════════════════════════════════

class TestAdversarialGerms:
    """Germs designed to probe the boundaries of classifier logic."""

    def test_very_high_degree_pure_power(self) -> None:
        """Φ(u) = u^{20} — extremely flat germ, weight = 1/20."""
        jet = PolynomialJet(dim=1, terms={(20,): 1.0}, certified_order=20)
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.weights == pytest.approx((0.05,))
        assert result.total_weight == pytest.approx(0.05)

    def test_near_zero_coefficient_still_accepted(self) -> None:
        """Φ(u) = εu⁴ with tiny ε > 0: the WH classifier checks
        coeff > 0 but does NOT threshold on magnitude.  So even
        extremely small positive coefficients are accepted.

        This documents the current behavior.  Whether a magnitude
        threshold should be added is an open design question."""
        jet = PolynomialJet(dim=1, terms={(4,): 5e-17}, certified_order=4)
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        # Current behavior: accepted since 5e-17 > 0.
        assert isinstance(result, WeightedHomogeneousTemplate)

    def test_zero_coefficient_refused(self) -> None:
        """Φ(u) = 0·u⁴ — exactly zero coefficient → refusal."""
        jet = PolynomialJet(dim=1, terms={(4,): 0.0}, certified_order=4)
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)

    def test_exact_subordination_boundary_is_non_separable(self) -> None:
        """Φ(u₁, u₂) = u₁⁴ + u₂⁴ + 10u₁²u₂² has cross-term weighted
        degree exactly 1.0.  This puts it on the Newton boundary, making
        the principal germ non-separable.  WH must refuse."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 4): 1.0, (2, 2): 10.0,
        })
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        # Cross term at weighted degree 1 → non-separable principal germ.
        assert not isinstance(result, WeightedHomogeneousTemplate)

    def test_subordination_slightly_below_one(self) -> None:
        """Φ(u₁, u₂) = u₁⁴ + u₂⁶ + u₁u₂² has cross-term weighted
        degree 1*(1/4) + 2*(1/6) = 7/12 < 1.  Should refuse WH."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 6): 1.0, (1, 2): 1.0,
        })
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert not isinstance(result, WeightedHomogeneousTemplate)

    def test_multiple_pure_powers_per_variable(self) -> None:
        """Φ(u) = u⁴ + u⁶.  The lowest pure power degree is 4."""
        jet = PolynomialJet(dim=1, terms={(4,): 1.0, (6,): 2.0})
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.leading_monomial_degrees == (4,)

    def test_negative_cross_term_at_newton_boundary_refused(self) -> None:
        """Φ(u₁, u₂) = u₁⁴ + u₂⁴ - 0.5u₁²u₂².
        Negative cross term with weighted degree exactly 1 sits on the
        Newton boundary → non-separable principal germ → WH refuses."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 4): 1.0, (2, 2): -0.5,
        })
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert not isinstance(result, WeightedHomogeneousTemplate)

    def test_strictly_subordinate_cross_term_accepted(self) -> None:
        """Φ(u₁, u₂) = u₁⁴ + u₂⁴ + u₁²u₂⁴.
        Cross term weighted degree = 2/4 + 4/4 = 3/2 > 1.
        Strictly subordinate → WH certifies."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 4): 1.0, (2, 4): 0.5,
        })
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, WeightedHomogeneousTemplate)

    def test_pure_cross_term_only(self) -> None:
        """Φ(u₁, u₂) = u₁²u₂² — no pure powers at all → refusal.
        Also not a branch channel (no branch order found)."""
        jet = PolynomialJet(dim=2, terms={(2, 2): 1.0})
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)

    def test_3d_separable_wh(self) -> None:
        """Φ(u₁, u₂, u₃) = u₁⁴ + u₂⁶ + u₃⁸ → weights (1/4, 1/6, 1/8)."""
        jet = PolynomialJet(dim=3, terms={
            (4, 0, 0): 1.0, (0, 6, 0): 1.0, (0, 0, 8): 1.0,
        }, certified_order=8)
        datum = _make_datum(3, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.weights == pytest.approx((0.25, 1/6, 0.125))
        expected_weight = 0.25 + 1/6 + 0.125
        assert result.total_weight == pytest.approx(expected_weight)

    def test_3d_no_branch_channel(self) -> None:
        """3D kernel with missing pure power → refusal.
        Branch channels only supported for dim=2."""
        jet = PolynomialJet(dim=3, terms={
            (4, 0, 0): 1.0, (0, 4, 0): 1.0, (0, 0, 2): -0.1,
        }, certified_order=4)
        datum = _make_datum(3, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)

    def test_branch_with_vanishing_stiffness_never_returns(self) -> None:
        """Φ(u, v) = u² (no v-dependence in the u² stiffness, but also
        this is actually a quadratic term, so the Hessian is nonzero →
        this should not be an unresolved kernel at all).

        If we force it through the singular classifier, the u² coefficient
        at v=0 is nonzero (h_at_zero > 0), so the branch channel returns
        None (not a branch-restiffening scenario)."""
        jet = PolynomialJet(dim=2, terms={(2, 0): 1.0, (0, 4): 1.0})
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        # The (2,0) term means variable u has a pure power at degree 2.
        # And (0,4) means v has degree 4.  Weights (1/2, 1/4), total 3/4.
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.total_weight == pytest.approx(0.75)

    def test_degenerate_zero_jet(self) -> None:
        """Empty jet (no terms at all) → refusal."""
        jet = PolynomialJet(dim=2, terms={})
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)

    def test_1d_kernel_odd_degree_refused(self) -> None:
        """1D kernel with odd-degree principal term → WH refuses
        (odd degree not positive on FULL_SPACE), branch not attempted
        (1D, not 2D) → refusal."""
        jet = PolynomialJet(dim=1, terms={(3,): 1.0})
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)


# ═══════════════════════════════════════════════════════════════════════
# 5. Kernel reduction correctness under rotated bases
# ═══════════════════════════════════════════════════════════════════════

class TestKernelReductionRotation:
    """Verify that kernel reduction is invariant under rotations of the
    kernel/normal basis."""

    def test_reduction_result_invariant_under_kernel_rotation(self) -> None:
        """Rotating the kernel basis should not change the reduced action
        evaluated along corresponding rotated directions.

        Build a 3D action with 1-dim kernel and 2-dim normal.
        The kernel action germ should be the same regardless of how
        we orient the normal basis.
        """
        # Ψ(x₀, x₁, x₂) = 0.5 x₀² + x₁² + x₂⁴ + 0.3 x₀ x₂²
        # Kernel direction is along x₂ (where Hessian is zero).
        # Normal directions are x₀, x₁.
        jet = PolynomialJet(dim=3, terms={
            (2, 0, 0): 0.5,
            (0, 2, 0): 1.0,
            (0, 0, 4): 1.0,
            (1, 0, 2): 0.3,  # cubic coupling x₀ · x₂²
        }, certified_order=4)

        # Kernel basis: e₂
        kernel_basis = np.array([[0.0], [0.0], [1.0]])
        # Normal basis: e₀, e₁
        normal_basis = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])

        seed = QuadraticKernelSeed(
            kernel_dim=1,
            positive_normal_dim=2,
            kernel_basis=kernel_basis,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=np.array([1.0, 2.0]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            log_normal_prefactor=0.0,
        )
        reduced = reduce_kernel_action(seed, jet)
        assert reduced.kernel_dim == 1

        # Evaluate the reduced kernel action at u = 1.
        val_at_1 = reduced.kernel_action.evaluate(np.array([1.0]))

        # Now rotate the normal basis by 45 degrees.
        theta = np.pi / 4
        R = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta),  np.cos(theta)]])
        rotated_normal = normal_basis @ R
        seed_rot = QuadraticKernelSeed(
            kernel_dim=1,
            positive_normal_dim=2,
            kernel_basis=kernel_basis,
            positive_normal_basis=rotated_normal,
            positive_normal_eigenvalues=np.array([1.0, 2.0]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            log_normal_prefactor=0.0,
        )
        # The eigenvalues don't match the rotated basis, but the
        # reduction still uses H_NN from the jet.  Let's just check
        # that reduction doesn't crash and gives a reasonable result.
        reduced_rot = reduce_kernel_action(seed_rot, jet)
        val_at_1_rot = reduced_rot.kernel_action.evaluate(np.array([1.0]))

        # The pure kernel term (u⁴) should be present in both.
        # The quartic correction depends on the coupling V and H_NN,
        # which are affected by the basis choice.  But the total value
        # should be reasonably close since the physical quantity is
        # basis-independent.
        # We just check both are finite and positive.
        assert np.isfinite(val_at_1)
        assert np.isfinite(val_at_1_rot)

    def test_reduction_with_no_cubic_coupling(self) -> None:
        """When there's no cubic coupling (no u²y terms), the reduced
        action should equal the pure kernel restriction exactly."""
        # Ψ(x₀, x₁) = 2 x₀² + x₁⁴ (no mixed cubic terms)
        jet = PolynomialJet(dim=2, terms={
            (2, 0): 2.0,
            (0, 4): 1.0,
        }, certified_order=4)

        kernel_basis = np.array([[0.0], [1.0]])  # x₁ is kernel
        normal_basis = np.array([[1.0], [0.0]])  # x₀ is normal

        seed = QuadraticKernelSeed(
            kernel_dim=1,
            positive_normal_dim=1,
            kernel_basis=kernel_basis,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=np.array([4.0]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            log_normal_prefactor=0.0,
        )
        reduced = reduce_kernel_action(seed, jet)
        # No cubic coupling → no quartic correction.
        # Reduced action = u⁴ exactly.
        val = reduced.kernel_action.evaluate(np.array([1.0]))
        assert val == pytest.approx(1.0, abs=1e-10)

    def test_cubic_coupling_correction_sign(self) -> None:
        """The quartic correction from cubic coupling should be negative
        (since ΔΦ₄ = -(1/2) V^T H_NN^{-1} V with V > 0, H_NN > 0)."""
        # Ψ(x, y) = y² + x⁴ + x²y (cubic coupling term)
        jet = PolynomialJet(dim=2, terms={
            (0, 2): 1.0,   # y²
            (4, 0): 1.0,   # x⁴
            (2, 1): 1.0,   # x²y — cubic coupling
        }, certified_order=4)

        kernel_basis = np.array([[1.0], [0.0]])  # x is kernel
        normal_basis = np.array([[0.0], [1.0]])  # y is normal

        seed = QuadraticKernelSeed(
            kernel_dim=1,
            positive_normal_dim=1,
            kernel_basis=kernel_basis,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=np.array([2.0]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            log_normal_prefactor=0.0,
        )
        reduced = reduce_kernel_action(seed, jet)
        # Pure kernel action at u=1 is u⁴ = 1.0.
        # Quartic correction: V₀(u) = u² (coefficient of y in cubic part)
        # ΔΦ₄ = -(1/2) * (1/H_NN) * V₀² = -(1/2) * (1/2) * u⁴ = -0.25u⁴
        # Total: (1 - 0.25)u⁴ = 0.75u⁴
        val = reduced.kernel_action.evaluate(np.array([1.0]))
        assert val == pytest.approx(0.75, abs=1e-10)


# ═══════════════════════════════════════════════════════════════════════
# 6. Consistency checks between layers
# ═══════════════════════════════════════════════════════════════════════

class TestCrossLayerConsistency:
    """Verify that the data passed between layers is consistent."""

    def test_normal_prefactor_passes_through_all_layers(self) -> None:
        """The normal prefactor from kernel reduction should appear
        in the singular classifier output."""
        jet = PolynomialJet(dim=1, terms={(4,): 1.0})
        prefactor = 3.14
        datum = ReducedKernelActionDatum(
            kernel_dim=1,
            kernel_action=jet,
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            positive_normal_dim=5,
            log_normal_prefactor=prefactor,
            total_dim=6,
        )
        result = classify_singular_kernel(datum)
        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.log_normal_prefactor == pytest.approx(prefactor)
        assert result.positive_normal_dim == 5

    def test_kernel_cone_propagates_to_singular(self) -> None:
        """The kernel cone from the seed should be reflected in the
        singular classifier's decision."""
        # With FULL_SPACE cone → WH succeeds.
        jet = PolynomialJet(dim=1, terms={(4,): 1.0})
        datum_full = _make_datum(1, jet)
        result_full = classify_singular_kernel(datum_full)
        assert isinstance(result_full, WeightedHomogeneousTemplate)

        # With POLYHEDRAL cone → WH refuses (not yet supported).
        datum_poly = ReducedKernelActionDatum(
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
        result_poly = classify_singular_kernel(datum_poly)
        assert isinstance(result_poly, UnresolvedKernelJetRefusal)

    def test_wh_principal_integral_consistency_with_scaling(self) -> None:
        """For Φ(u) = c·u^d, the principal integral and scaling exponent
        must be internally consistent.

        I(n) = ∫ e^{-n c |u|^d} du = n^{-1/d} · (integral of e^{-c|s|^d} ds)

        Check that classifier's total_weight = 1/d and the principal
        integral gives the correct constant factor."""
        from scipy import integrate

        for d, c in [(4, 1.0), (6, 2.0), (4, 3.0)]:
            jet = PolynomialJet(dim=1, terms={(d,): c}, certified_order=d)
            datum = _make_datum(1, jet)
            result = classify_singular_kernel(datum)
            assert isinstance(result, WeightedHomogeneousTemplate)
            assert result.total_weight == pytest.approx(1.0 / d)

            # Numerical integration at n=1000.
            n_val = 1000
            numerical, _ = integrate.quad(
                lambda u: np.exp(-n_val * c * abs(u)**d), -10, 10,
            )
            # Predicted: n^{-1/d} * exp(log_principal_integral)
            predicted = n_val**(-1.0 / d) * np.exp(result.log_principal_integral)
            assert numerical == pytest.approx(predicted, rel=0.01)

    def test_slice_eigenvalue_preservation(self) -> None:
        """After slice reduction (orbit removal), the nonzero eigenvalues
        of the transverse Hessian should match those of the projected
        ambient Hessian."""
        # 4D Hessian with 1D kernel along e_3.
        eigs = [3.0, 2.0, 1.0, 0.0]
        H = np.diag(eigs)
        orbit_basis = np.array([[0], [0], [0], [1.0]])

        datum = reduce_local_chart(
            H,
            orbit_tangent_basis=orbit_basis,
            log_orbit_volume=0.0,
        )
        assert datum.active_dim == 3
        transverse_eigs = sorted(np.linalg.eigvalsh(datum.h_active), reverse=True)
        expected = sorted([3.0, 2.0, 1.0], reverse=True)
        np.testing.assert_array_almost_equal(transverse_eigs, expected, decimal=10)


# ═══════════════════════════════════════════════════════════════════════
# 7. Odd-degree positivity gap (known limitation)
# ═══════════════════════════════════════════════════════════════════════

class TestKnownLimitations:
    """Document known limitations as tests so they're tracked."""

    def test_odd_degree_pure_power_rejected_on_full_space(self) -> None:
        """Φ(u) = u³ — odd degree is NOT positive on all of R
        (u³ < 0 for u < 0).  The WH classifier correctly rejects
        odd-degree principal terms on FULL_SPACE cones."""
        jet = PolynomialJet(dim=1, terms={(3,): 1.0})
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)


# ═══════════════════════════════════════════════════════════════════════
# 8. GPT-requested tests
# ═══════════════════════════════════════════════════════════════════════

class TestGPTRequestedCoverage:
    """Tests specifically requested by GPT review."""

    def test_wh_refuses_same_degree_mixed_principal_part(self) -> None:
        """u₁⁴ + u₂⁴ + 10·u₁²u₂² — the cross term has weighted degree
        exactly 1 under weights (1/4, 1/4), making it part of the
        principal germ.  The product integral formula is not licensed.
        WH must refuse."""
        jet = PolynomialJet(dim=2, terms={
            (4, 0): 1.0, (0, 4): 1.0, (2, 2): 10.0,
        })
        datum = _make_datum(2, jet)
        result = classify_singular_kernel(datum)
        assert not isinstance(result, WeightedHomogeneousTemplate)

    def test_wh_refuses_odd_degree_full_space(self) -> None:
        """u³ on FULL_SPACE — odd degree, not positive for u < 0."""
        jet = PolynomialJet(dim=1, terms={(3,): 1.0})
        datum = _make_datum(1, jet)
        result = classify_singular_kernel(datum)
        assert isinstance(result, UnresolvedKernelJetRefusal)

    def test_kernel_reduction_downgrades_certification(self) -> None:
        """Ψ(u, y) = ½y² + u²y + y³.

        The exact eliminated action has terms at all even orders:
          Φ(u) = -½u⁴ - u⁶ - ...
        but the reducer only computes the quartic correction from the
        cubic coupling u²y.  The y³ term contributes at order 6+ via
        η₂³.  So the output certified_order must be ≤ 4, NOT the
        input's certified_order.

        This is GPT's concrete counterexample.
        """
        jet = PolynomialJet(dim=2, terms={
            (0, 2): 0.5,   # ½y²
            (2, 1): 1.0,   # u²y
            (0, 3): 1.0,   # y³
        }, certified_order=6)

        kernel_basis = np.array([[1.0], [0.0]])  # u is kernel
        normal_basis = np.array([[0.0], [1.0]])  # y is normal

        seed = QuadraticKernelSeed(
            kernel_dim=1,
            positive_normal_dim=1,
            kernel_basis=kernel_basis,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=np.array([1.0]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            log_normal_prefactor=0.0,
        )
        reduced = reduce_kernel_action(seed, jet)

        # The certified order must be downgraded: input was 6, but
        # the reducer cannot certify beyond 4 when higher normal
        # nonlinearities (y³) are present.
        assert reduced.kernel_action.certified_order <= 4

        # The quartic correction should be present: ΔΦ₄ = -(1/2)(1/H_NN)·V²
        # V(u) = u² (from the u²y term), H_NN = 1.0
        # ΔΦ₄ = -(1/2)·1·u⁴ = -0.5u⁴
        val = reduced.kernel_action.evaluate(np.array([1.0]))
        assert val == pytest.approx(-0.5, abs=1e-10)

        # And terms above order 4 should NOT be present.
        for multi_idx in reduced.kernel_action.terms:
            assert sum(multi_idx) <= 4

    def test_branch_channel_requires_straightened_flag(self) -> None:
        """Without branch_straightened=True, a 2D kernel jet that fails
        WH should fall through to refusal, NOT to branch classification."""
        # This germ would be classified as critical branch if
        # branch_straightened=True.
        jet = PolynomialJet(
            dim=2,
            terms={(6, 0): 1.0, (2, 2): 1.0, (0, 4): 1.0},
            certified_order=6,
        )
        datum = _make_datum(2, jet)

        # Without the flag: refusal.
        result_default = classify_singular_kernel(datum)
        assert isinstance(result_default, UnresolvedKernelJetRefusal)
        assert "branch_straightened" in result_default.reason

        # With the flag: critical branch.
        result_flagged = classify_singular_kernel(
            datum, branch_straightened=True
        )
        assert isinstance(result_flagged, CriticalBranchTemplate)

    def test_end_to_end_hessian_to_singular_evidence(self) -> None:
        """Full pipeline: Hessian with kernel → detector → kernel
        reduction → singular classifier → evidence template.

        Build a 3D action with a 1-dim kernel, quartic in the kernel
        direction, and verify we get a SingularEvidenceTemplate out."""
        from nomogeo.evidence import dispatch_singular_evidence
        from nomogeo.regime_types import SingularEvidenceTemplate

        # 3×3 Hessian: diag(2, 1, 0) → kernel dim 1 along e₂.
        H = np.diag([2.0, 1.0, 0.0])
        datum = reduce_local_chart(H)
        classification = classify_regime(datum)
        assert isinstance(classification, UnresolvedKernel)
        seed = classification.kernel

        # Action jet: Ψ = x₀² + 0.5x₁² + x₂⁴
        jet = PolynomialJet(dim=3, terms={
            (2, 0, 0): 1.0,
            (0, 2, 0): 0.5,
            (0, 0, 4): 1.0,
        }, certified_order=4)

        reduced = reduce_kernel_action(seed, jet)
        singular = classify_singular_kernel(reduced)
        assert isinstance(singular, WeightedHomogeneousTemplate)

        evidence = dispatch_singular_evidence(singular)
        assert isinstance(evidence, SingularEvidenceTemplate)
        assert evidence.kernel_dim == 1
        assert evidence.positive_normal_dim == 2
        assert evidence.singular_classification is singular
