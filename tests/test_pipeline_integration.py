"""
Integration tests for the full Layer 3→4→5 pipeline (0.3.3).

Families:
  A: Spectral hygiene (indefinite → IndefiniteStationaryPoint)
  B: Coordinate correctness (active vs kn entry points)
  C: Sparse vs dense jet extraction
  D: Coefficient conventions (known analytic cases)
  E: Singular classifier path (end-to-end)
  F: Tolerance stress (eigenvalue crossing +ε → 0 → −ε under scaling)
"""
from __future__ import annotations

import numpy as np
import pytest

from nomogeo.regime import classify_regime, classify_from_hessian
from nomogeo.regime_types import (
    ConeKind,
    ConeSpec,
    IndefiniteStationaryPoint,
    ReducedLocalDatum,
    RegimeKind,
    RegularInterior,
    UnresolvedKernel,
)
from nomogeo.kernel_reduction import (
    PolynomialJet,
    ReducedKernelActionDatum,
    reduce_kernel_action,
    reduce_kernel_action_kn,
)
from nomogeo.empirical_kernel_jet import extract_sparse_kn_jet
from nomogeo.evidence import dispatch_evidence


def _datum(H, cone=None):
    d = H.shape[0]
    if cone is None:
        cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d)
    return ReducedLocalDatum(active_dim=d, h_active=H, cone=cone)


# ═══════════════════════════════════════════════════════════════════════
# Family A: Spectral hygiene
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyA:
    """Indefinite Hessians must be caught by the pre-detector, not
    conflated with genuine kernel cases."""

    def test_negative_definite_gives_indefinite(self):
        """H = -I (all eigenvalues negative)."""
        H = -np.eye(3)
        result = classify_regime(_datum(H))
        assert isinstance(result, IndefiniteStationaryPoint)
        assert result.negative_dim == 3
        assert result.positive_dim == 0
        assert result.min_eigenvalue < 0

    def test_mixed_indefinite_gives_indefinite(self):
        """H with both positive and negative eigenvalues — saddle point."""
        H = np.diag([2.0, -1.0, 3.0])
        result = classify_regime(_datum(H))
        assert isinstance(result, IndefiniteStationaryPoint)
        assert result.negative_dim == 1
        assert result.positive_dim == 2
        assert result.min_eigenvalue == pytest.approx(-1.0, abs=1e-10)

    def test_slightly_negative_still_caught(self):
        """H with one small negative eigenvalue — borderline case."""
        H = np.diag([1.0, 1.0, -0.01])
        result = classify_regime(_datum(H))
        assert isinstance(result, IndefiniteStationaryPoint)
        assert result.negative_dim == 1

    def test_psd_with_kernel_gives_unresolved(self):
        """H with one zero eigenvalue and all others positive — genuine kernel."""
        H = np.diag([2.0, 3.0, 0.0])
        result = classify_regime(_datum(H))
        assert isinstance(result, UnresolvedKernel)
        assert result.kernel.kernel_dim == 1

    def test_spd_gives_regular(self):
        """H positive definite — regular interior."""
        H = np.diag([1.0, 2.0, 3.0])
        result = classify_regime(_datum(H))
        assert isinstance(result, RegularInterior)

    def test_indefinite_evidence_is_refusal(self):
        """IndefiniteStationaryPoint dispatches to an evidence refusal."""
        H = np.diag([2.0, -1.0])
        result = classify_regime(_datum(H))
        assert isinstance(result, IndefiniteStationaryPoint)
        evidence = dispatch_evidence(result)
        assert "indefinite" in evidence.reason.lower()


# ═══════════════════════════════════════════════════════════════════════
# Family B: Coordinate correctness
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyB:
    """reduce_kernel_action (active coords) and reduce_kernel_action_kn
    (kn coords) must produce the same reduced germ when given
    equivalent inputs."""

    def _make_seed(self, r, p, eigenvalues):
        """Create a QuadraticKernelSeed with standard basis."""
        d = r + p
        from nomogeo.regime_types import QuadraticKernelSeed
        kernel_basis = np.eye(d)[:, :r]
        normal_basis = np.eye(d)[:, r:]
        log_prefactor = 0.5 * p * np.log(2 * np.pi) - 0.5 * np.sum(np.log(eigenvalues))
        return QuadraticKernelSeed(
            kernel_dim=r,
            kernel_basis=kernel_basis,
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=r),
            positive_normal_dim=p,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=eigenvalues,
            log_normal_prefactor=log_prefactor,
        )

    def test_active_vs_kn_standard_basis(self):
        """When the eigenbasis IS the standard basis, both entry points
        should give identical results since no rotation is needed."""
        r, p = 1, 2
        eigs = np.array([2.0, 5.0])
        seed = self._make_seed(r, p, eigs)

        # Jet in active coords (= kn coords when basis is identity)
        jet = PolynomialJet(
            dim=3,
            terms={(4,0,0): 1.0, (2,1,0): 3.0, (0,2,0): 1.0, (0,0,2): 2.5},
            certified_order=4,
        )

        result_active = reduce_kernel_action(seed, jet)
        result_kn = reduce_kernel_action_kn(seed, jet)

        # Both should produce the same reduced germ
        for idx in set(result_active.kernel_action.terms) | set(result_kn.kernel_action.terms):
            c_a = result_active.kernel_action.terms.get(idx, 0.0)
            c_k = result_kn.kernel_action.terms.get(idx, 0.0)
            assert abs(c_a - c_k) < 1e-10, f"Mismatch at {idx}: {c_a} vs {c_k}"

    def test_active_vs_kn_rotated_basis(self):
        """When eigenbasis differs from standard basis, the active-coords
        entry point must rotate correctly, giving the same result as a
        pre-rotated kn-coords jet passed to the kn entry point."""
        r, p = 1, 1
        eigs = np.array([4.0])

        # Non-trivial rotation: kernel = (1,1)/sqrt(2), normal = (1,-1)/sqrt(2)
        s2 = 1.0 / np.sqrt(2)
        kernel_basis = np.array([[s2], [s2]])
        normal_basis = np.array([[s2], [-s2]])

        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=r,
            kernel_basis=kernel_basis,
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=r),
            positive_normal_dim=p,
            positive_normal_basis=normal_basis,
            positive_normal_eigenvalues=eigs,
            log_normal_prefactor=0.5 * np.log(2 * np.pi) - 0.5 * np.log(4.0),
        )

        # Action in ACTIVE coords: Psi(x1, x2) = x1^4 + x2^4 + 2*(x1^2 + x2^2)
        jet_active = PolynomialJet(
            dim=2,
            terms={(4,0): 1.0, (0,4): 1.0, (2,0): 1.0, (0,2): 1.0},
            certified_order=4,
        )
        result_active = reduce_kernel_action(seed, jet_active)

        # Manually rotate to kn coords and use kn entry point
        jet_kn = jet_active.restrict_to_subspace(
            np.hstack([kernel_basis, normal_basis])
        )
        result_kn = reduce_kernel_action_kn(seed, jet_kn)

        for idx in set(result_active.kernel_action.terms) | set(result_kn.kernel_action.terms):
            c_a = result_active.kernel_action.terms.get(idx, 0.0)
            c_k = result_kn.kernel_action.terms.get(idx, 0.0)
            assert abs(c_a - c_k) < 1e-10, f"Mismatch at {idx}: {c_a} vs {c_k}"


# ═══════════════════════════════════════════════════════════════════════
# Family C: Sparse vs dense jet extraction
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyC:
    """The sparse kn-jet extractor must agree with a dense extraction
    on the terms it claims to compute."""

    def test_sparse_matches_known_polynomial(self):
        """Extract jet of a known polynomial and verify coefficients."""
        # Psi(u, y1, y2) = u^4 + 2*u^2*y1 + y1^2 + 3*y2^2
        def action(theta):
            u, y1, y2 = theta
            return u**4 + 2*u**2*y1 + y1**2 + 3*y2**2

        theta_star = np.zeros(3)
        K = np.array([[1],[0],[0]], dtype=float)
        N = np.array([[0,0],[1,0],[0,1]], dtype=float)

        jet = extract_sparse_kn_jet(action, theta_star, K, N, order=4, h=1e-4)

        expected = {
            (4, 0, 0): 1.0,
            (2, 1, 0): 2.0,
            (0, 2, 0): 1.0,
            (0, 0, 2): 3.0,
        }
        for idx, exp_val in expected.items():
            got = jet.terms.get(idx, 0.0)
            assert abs(got - exp_val) < 1e-3, f"At {idx}: expected {exp_val}, got {got}"

    def test_sparse_certified_order_is_honest(self):
        """Sparse extractor must set certified_order=2, not the probing order."""
        def action(theta):
            return theta[0]**4 + 2*theta[0]**2*theta[1] + theta[1]**2

        theta_star = np.zeros(2)
        K = np.array([[1.0],[0.0]])
        N = np.array([[0.0],[1.0]])

        jet = extract_sparse_kn_jet(action, theta_star, K, N, order=4, h=1e-4)
        assert jet.certified_order == 2, (
            f"Sparse jet should have certified_order=2, got {jet.certified_order}"
        )
        assert jet.metadata is not None
        assert jet.metadata.get("certified") is False
        assert "sparse" in jet.metadata.get("extraction", "")

    def test_sparse_reducer_does_not_overcertify(self):
        """When sparse jet (certified_order=2) has cubic coupling, the
        reducer must NOT compute the quartic correction, because the
        jet cannot certify that no higher normal nonlinearities exist."""
        # Psi = u^4 + 3*u^2*y + 2*y^2 — has cubic coupling
        def action(theta):
            return theta[0]**4 + 3*theta[0]**2*theta[1] + 2*theta[1]**2

        theta_star = np.zeros(2)
        K = np.array([[1.0],[0.0]])
        N = np.array([[0.0],[1.0]])

        jet = extract_sparse_kn_jet(action, theta_star, K, N, order=4, h=1e-4)
        assert jet.certified_order == 2  # honest

        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=1,
            kernel_basis=K,
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            positive_normal_dim=1,
            positive_normal_basis=N,
            positive_normal_eigenvalues=np.array([4.0]),
            log_normal_prefactor=0.5*np.log(2*np.pi) - 0.5*np.log(4.0),
        )
        reduced = reduce_kernel_action_kn(seed, jet)

        # The reducer should NOT apply the quartic correction because
        # certified_order < 4.  So the output u^4 coefficient should be
        # the raw pure-kernel value (1.0), NOT the corrected value (-0.125).
        u4 = reduced.kernel_action.terms.get((4,), 0.0)
        assert abs(u4 - 1.0) < 1e-3, (
            f"Sparse jet should yield uncorrected u^4=1.0, got {u4}"
        )
        # And the output certified_order should also be 2.
        assert reduced.kernel_action.certified_order == 2

    def test_sparse_plus_kn_reduction_matches_analytic(self):
        """Sparse extraction → kn reduction must give the analytically
        known reduced germ."""
        # A1 singularity: Psi = u^4 + y^2
        # Reduced germ: Phi = u^4 (no coupling, so no correction)
        def action(theta):
            return theta[0]**4 + theta[1]**2

        theta_star = np.zeros(2)
        K = np.array([[1.0],[0.0]])
        N = np.array([[0.0],[1.0]])

        jet = extract_sparse_kn_jet(action, theta_star, K, N, order=4, h=1e-4)

        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=1,
            kernel_basis=K,
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            positive_normal_dim=1,
            positive_normal_basis=N,
            positive_normal_eigenvalues=np.array([2.0]),
            log_normal_prefactor=0.5*np.log(2*np.pi) - 0.5*np.log(2.0),
        )

        reduced = reduce_kernel_action_kn(seed, jet)
        assert abs(reduced.kernel_action.terms.get((4,), 0.0) - 1.0) < 1e-3


# ═══════════════════════════════════════════════════════════════════════
# Family D: Coefficient conventions
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyD:
    """Known analytic singularities must give the correct reduced germ
    coefficients under the 1/α! normalisation convention."""

    def test_pure_quartic_no_coupling(self):
        """Psi = a*u^4 + b*y^2 → Phi = a*u^4 (no cubic coupling)."""
        a, b = 2.5, 3.0
        jet = PolynomialJet(dim=2, terms={(4,0): a, (0,2): b}, certified_order=4)
        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=1,
            kernel_basis=np.array([[1.0],[0.0]]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            positive_normal_dim=1,
            positive_normal_basis=np.array([[0.0],[1.0]]),
            positive_normal_eigenvalues=np.array([2*b]),
            log_normal_prefactor=0.5*np.log(2*np.pi) - 0.5*np.log(2*b),
        )
        reduced = reduce_kernel_action_kn(seed, jet)
        assert abs(reduced.kernel_action.terms.get((4,), 0.0) - a) < 1e-10

    def test_quartic_with_cubic_coupling(self):
        """Psi = a*u^4 + c*u^2*y + b*y^2
        → Phi = a*u^4 - (1/2) * c^2 / (2b) * u^4
             = (a - c^2/(4b)) * u^4"""
        a, b, c = 1.0, 2.0, 3.0
        jet = PolynomialJet(
            dim=2,
            terms={(4,0): a, (0,2): b, (2,1): c},
            certified_order=4,
        )
        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=1,
            kernel_basis=np.array([[1.0],[0.0]]),
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=1),
            positive_normal_dim=1,
            positive_normal_basis=np.array([[0.0],[1.0]]),
            positive_normal_eigenvalues=np.array([2*b]),
            log_normal_prefactor=0.5*np.log(2*np.pi) - 0.5*np.log(2*b),
        )
        reduced = reduce_kernel_action_kn(seed, jet)
        expected_u4 = a - c**2 / (4*b)
        got_u4 = reduced.kernel_action.terms.get((4,), 0.0)
        assert abs(got_u4 - expected_u4) < 1e-10, f"Expected {expected_u4}, got {got_u4}"

    def test_2d_kernel_pure_quartic(self):
        """Psi = u1^4 + u2^4 + y^2 in R^3, r=2, p=1.
        No coupling → Phi = u1^4 + u2^4."""
        jet = PolynomialJet(
            dim=3,
            terms={(4,0,0): 1.0, (0,4,0): 1.0, (0,0,2): 2.0},
            certified_order=4,
        )
        from nomogeo.regime_types import QuadraticKernelSeed
        seed = QuadraticKernelSeed(
            kernel_dim=2,
            kernel_basis=np.eye(3)[:, :2],
            kernel_cone=ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=2),
            positive_normal_dim=1,
            positive_normal_basis=np.eye(3)[:, 2:],
            positive_normal_eigenvalues=np.array([4.0]),
            log_normal_prefactor=0.5*np.log(2*np.pi) - 0.5*np.log(4.0),
        )
        reduced = reduce_kernel_action_kn(seed, jet)
        assert abs(reduced.kernel_action.terms.get((4,0), 0.0) - 1.0) < 1e-10
        assert abs(reduced.kernel_action.terms.get((0,4), 0.0) - 1.0) < 1e-10


# ═══════════════════════════════════════════════════════════════════════
# Family E: Singular classifier path
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyE:
    """End-to-end: synthetic action → Layer 3 → Layer 4 → Layer 5."""

    def test_a1_singularity_full_path(self):
        """A1 singularity u^4 + y^2: Layer 3 detects kernel,
        Layer 4 reduces, Layer 5 classifies."""
        from nomogeo.singular import classify_singular_kernel

        # Build Hessian: kernel in direction 1, normal eigenvalue 2
        H = np.diag([0.0, 2.0])
        result = classify_from_hessian(H)
        assert isinstance(result, UnresolvedKernel)
        seed = result.kernel

        # Build jet in kn-coords (jet is already in eigenbasis = standard basis)
        jet = PolynomialJet(
            dim=2,
            terms={(4,0): 1.0, (0,2): 1.0},
            certified_order=4,
        )

        # Layer 4: reduce
        reduced = reduce_kernel_action_kn(seed, jet)
        assert reduced.kernel_dim == 1
        assert abs(reduced.kernel_action.terms.get((4,), 0.0) - 1.0) < 1e-10

        # Layer 5: classify
        singular = classify_singular_kernel(reduced)
        # The classifier should recognise this as a weighted homogeneous germ
        # or return a refusal if the coefficient is not certified.
        # Either way, it should not crash.
        assert singular is not None

    def test_degenerate_quartic_with_coupling(self):
        """u^4 + 3*u^2*y + 2*y^2: Layer 4 should produce a quartic
        reduced germ with the correct coupling correction."""
        H = np.diag([0.0, 4.0])  # H_NN = 4.0 (since y^2 coeff = 2 → H_NN_jj = 2*2 = 4)
        result = classify_from_hessian(H)
        assert isinstance(result, UnresolvedKernel)
        seed = result.kernel

        jet = PolynomialJet(
            dim=2,
            terms={(4,0): 1.0, (2,1): 3.0, (0,2): 2.0},
            certified_order=4,
        )
        reduced = reduce_kernel_action_kn(seed, jet)

        # Expected: u^4 coeff = 1 - 3^2/(4*2) = 1 - 9/8 = -1/8
        expected = 1.0 - 9.0/8.0
        got = reduced.kernel_action.terms.get((4,), 0.0)
        assert abs(got - expected) < 1e-10, f"Expected {expected}, got {got}"


# ═══════════════════════════════════════════════════════════════════════
# Family F: Tolerance stress
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyF:
    """Sweep one eigenvalue through +ε → 0 → −ε under different overall
    Hessian scales, and verify the three-way split is stable.

    The spectral tolerance is:
        ε = max(atol, rtol · max(1, ||H||_op))

    These tests verify:
      - Classification is invariant under uniform positive scaling of H
        (the split should track the relative spectrum, not absolute values)
      - The borderline region (|λ| ≈ ε) is narrow and predictable
      - Large-scale Hessians don't shift the boundary in unexpected ways
    """

    def _classify_min_eig(self, bulk_eigenvalue, test_eigenvalue, d=4):
        """Build a d×d Hessian with (d-1) copies of bulk_eigenvalue
        and one copy of test_eigenvalue, classify it."""
        from nomogeo.validation import Tolerances
        eigs = np.array([bulk_eigenvalue] * (d - 1) + [test_eigenvalue])
        # Random rotation to avoid axis-alignment artefacts
        rng = np.random.RandomState(42)
        Q, _ = np.linalg.qr(rng.randn(d, d))
        H = Q @ np.diag(eigs) @ Q.T
        H = 0.5 * (H + H.T)
        tol = Tolerances(atol=1e-8, rtol=1e-6)
        return classify_regime(_datum(H), tolerances=tol)

    def test_clearly_positive_is_regular(self):
        """All eigenvalues well above ε → RegularInterior."""
        result = self._classify_min_eig(5.0, 0.1)
        assert isinstance(result, RegularInterior)

    def test_clearly_negative_is_indefinite(self):
        """One eigenvalue well below -ε → IndefiniteStationaryPoint."""
        result = self._classify_min_eig(5.0, -0.1)
        assert isinstance(result, IndefiniteStationaryPoint)
        assert result.negative_dim == 1

    def test_exact_zero_is_kernel(self):
        """One eigenvalue exactly zero → UnresolvedKernel."""
        result = self._classify_min_eig(5.0, 0.0)
        assert isinstance(result, UnresolvedKernel)
        assert result.kernel.kernel_dim == 1

    def test_scaling_invariance_regular(self):
        """Scaling H by α > 0 must not change RegularInterior classification.
        Test with α = 1, 100, 1e-4."""
        for scale in [1.0, 100.0, 1e-4]:
            result = self._classify_min_eig(5.0 * scale, 0.1 * scale)
            assert isinstance(result, RegularInterior), (
                f"scale={scale}: expected RegularInterior, got {type(result).__name__}"
            )

    def test_scaling_invariance_indefinite(self):
        """Scaling H by α > 0 must not change IndefiniteStationaryPoint."""
        for scale in [1.0, 100.0, 1e-4]:
            result = self._classify_min_eig(5.0 * scale, -0.1 * scale)
            assert isinstance(result, IndefiniteStationaryPoint), (
                f"scale={scale}: expected Indefinite, got {type(result).__name__}"
            )

    def test_scaling_invariance_kernel(self):
        """Scaling H by α > 0 must not change UnresolvedKernel."""
        for scale in [1.0, 100.0, 1e-4]:
            result = self._classify_min_eig(5.0 * scale, 0.0)
            assert isinstance(result, UnresolvedKernel), (
                f"scale={scale}: expected UnresolvedKernel, got {type(result).__name__}"
            )

    def test_borderline_positive_is_regular_or_kernel(self):
        """Eigenvalue at +2ε should be regular (positive).
        Eigenvalue at +0.5ε should be kernel (within tolerance)."""
        from nomogeo.validation import Tolerances
        bulk = 5.0
        tol = Tolerances(atol=1e-8, rtol=1e-6)
        eps = max(tol.atol, tol.rtol * max(1.0, bulk))  # ≈ 5e-6

        # Well above ε → regular
        result = self._classify_min_eig(bulk, 10 * eps)
        assert isinstance(result, RegularInterior)

        # Within ε → kernel
        result = self._classify_min_eig(bulk, 0.5 * eps)
        assert isinstance(result, UnresolvedKernel)

    def test_borderline_negative_is_indefinite(self):
        """Eigenvalue at -2ε should be indefinite."""
        from nomogeo.validation import Tolerances
        bulk = 5.0
        tol = Tolerances(atol=1e-8, rtol=1e-6)
        eps = max(tol.atol, tol.rtol * max(1.0, bulk))

        result = self._classify_min_eig(bulk, -2 * eps)
        assert isinstance(result, IndefiniteStationaryPoint)

    def test_large_hessian_scale_threshold_tracks(self):
        """When the Hessian scale is 1e6, the tolerance scales to
        ~1.0 (rtol * 1e6).  A test eigenvalue of 0.1 should be
        within tolerance (kernel), not regular."""
        from nomogeo.validation import Tolerances
        bulk = 1e6
        tol = Tolerances(atol=1e-8, rtol=1e-6)
        eps = max(tol.atol, tol.rtol * max(1.0, bulk))  # ≈ 1.0

        # 0.1 < eps ≈ 1.0 → kernel
        result = self._classify_min_eig(bulk, 0.1)
        assert isinstance(result, UnresolvedKernel), (
            f"Expected kernel at scale=1e6 with test_eig=0.1, "
            f"got {type(result).__name__}"
        )

        # 5.0 > eps ≈ 1.0 → regular
        result = self._classify_min_eig(bulk, 5.0)
        assert isinstance(result, RegularInterior)

        # -5.0 < -eps → indefinite
        result = self._classify_min_eig(bulk, -5.0)
        assert isinstance(result, IndefiniteStationaryPoint)


# ═══════════════════════════════════════════════════════════════════════
# Family G: Stage A — Exact regime benchmark ladder
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyG:
    """Exact germs where the correct singular channel is known analytically.

    These are the Stage A benchmarks from the regime benchmark ladder.
    Each test constructs a ReducedKernelActionDatum with a known kernel
    action jet and verifies that classify_singular_kernel returns the
    expected channel (or refusal).

    Cases:
      G1: Clean weighted-homogeneous quartic (r=1) → WH channel
      G2: Critical quadratic branch-restiffening (r=2) → CriticalBranchTemplate
      G3: Under-certified germ (certified_order=2) → refusal (not WH)
      G4: Multi-dimensional separable WH (r=2) → WH channel
    """

    def _make_reduced_datum(self, r, p, jet, *, cone_kind=ConeKind.FULL_SPACE):
        """Build a ReducedKernelActionDatum from a kernel jet."""
        eigs = np.array([2.0] * p)
        log_prefactor = 0.5 * p * np.log(2 * np.pi) - 0.5 * np.sum(np.log(eigs))
        return ReducedKernelActionDatum(
            kernel_dim=r,
            kernel_action=jet,
            kernel_cone=ConeSpec(kind=cone_kind, ambient_dim=r),
            positive_normal_dim=p,
            log_normal_prefactor=log_prefactor,
            total_dim=r + p,
        )

    # ── G1: Clean WH quartic (r=1) ─────────────────────────────────

    def test_g1_clean_wh_quartic(self):
        """Φ(u) = u⁴ with certified_order=4 → WeightedHomogeneousTemplate.

        Expected: weight = 1/4, total_weight = 0.25, degree = (4,).
        The principal kernel integral factors as
        ∫_{-∞}^{∞} e^{-|s|⁴} ds = 2·Γ(1/4)/4.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate

        jet = PolynomialJet(dim=1, terms={(4,): 1.0}, certified_order=4)
        datum = self._make_reduced_datum(r=1, p=2, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, WeightedHomogeneousTemplate), (
            f"Expected WeightedHomogeneousTemplate, got {type(result).__name__}"
        )
        assert result.weights == (0.25,)
        assert result.total_weight == pytest.approx(0.25)
        assert result.leading_monomial_degrees == (4,)
        assert result.kernel_dim == 1
        assert result.positive_normal_dim == 2

        # Verify the principal integral: 2·Γ(1/4)/4
        from math import lgamma
        expected_log = np.log(2.0) + lgamma(0.25) - np.log(4.0)
        assert result.log_principal_integral == pytest.approx(expected_log, rel=1e-10)

    def test_g1_wh_with_subordinate_terms(self):
        """Φ(u) = u⁴ + 0.5·u⁶ with certified_order=6 → still WH.

        u⁶ has weighted degree 6·(1/4) = 1.5 > 1, so it is strictly
        subordinate to the u⁴ principal part.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate

        jet = PolynomialJet(
            dim=1, terms={(4,): 1.0, (6,): 0.5}, certified_order=6
        )
        datum = self._make_reduced_datum(r=1, p=1, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.weights == (0.25,)
        assert result.leading_monomial_degrees == (4,)

    # ── G2: Critical quadratic branch-restiffening (r=2) ───────────

    def test_g2_critical_branch(self):
        """Φ(u,v) = u⁴ + u²v² + v⁶ with certified_order=6, branch_straightened=True.

        WH analysis:
          Pure powers: u⁴ (ω₁=1/4) and v⁶ (ω₂=1/6).
          Mixed term u²v² has weighted degree = 2·(1/4) + 2·(1/6) = 5/6 < 1.
          Subordination fails → WH refuses.

        Branch analysis (c=2):
          u² coefficient as function of v: λ(v) = v².
          λ(0) = 0 → stiffness vanishes.  stiffness_ord = 2, H₀ = 1.0.
          Branch action: Φ(0,v) = v⁶ → d = 6, B₀ = 1.0.
          Higher transverse: Φ(u,0) = u⁴ → a = 4, C₀ = 1.0.
          κ = (6·(4-2) - 2·4) / (2·4) = (12 - 8)/8 = 0.5.

        Expected: CriticalBranchTemplate with κ = 0.5.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import CriticalBranchTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (2, 2): 1.0, (0, 6): 1.0},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)

        assert isinstance(result, CriticalBranchTemplate), (
            f"Expected CriticalBranchTemplate, got {type(result).__name__}: "
            f"{getattr(result, 'reason', '')}"
        )
        assert result.kappa == pytest.approx(0.5)
        assert result.branch_datum.branch_order == 2
        assert result.branch_datum.stiffness_vanishing_order == 2
        assert result.branch_datum.branch_action_degree == 6
        assert result.branch_datum.higher_transverse_degree == 4

    def test_g2_branch_without_straightened_flag_refuses(self):
        """Same germ as g2_critical_branch, but branch_straightened=False.

        WH fails (subordination), branch not attempted → refusal.
        This validates the safety guard: we don't enter the branch
        channel unless the caller certifies straightened coordinates.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import UnresolvedKernelJetRefusal

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (2, 2): 1.0, (0, 6): 1.0},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=False)

        assert isinstance(result, UnresolvedKernelJetRefusal)
        assert "branch" in result.reason.lower()

    # ── G3: Under-certified germ → refusal ─────────────────────────

    def test_g3_undercertified_quartic_refuses(self):
        """Φ(u) = u⁴ with certified_order=2 → must REFUSE, not WH.

        The u⁴ term is present but the jet is only certified to degree 2.
        The quartic coefficient may be incomplete (missing perturbative
        corrections from kernel reduction).  The classifier must check
        certified_order and refuse when the principal part lives above it.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import UnresolvedKernelJetRefusal

        jet = PolynomialJet(dim=1, terms={(4,): 1.0}, certified_order=2)
        datum = self._make_reduced_datum(r=1, p=1, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, UnresolvedKernelJetRefusal), (
            f"Under-certified jet (certified_order=2 with quartic principal "
            f"part) must refuse, but got {type(result).__name__}"
        )

    def test_g3_undercertified_branch_refuses(self):
        """2D germ with certified_order=2 → branch channel must also refuse.

        Even with branch_straightened=True, the branch classifier needs
        to read structure at degree ≥ 4, so certified_order=2 is too low.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import UnresolvedKernelJetRefusal

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (2, 2): 1.0, (0, 6): 1.0},
            certified_order=2,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)

        assert isinstance(result, UnresolvedKernelJetRefusal), (
            f"Under-certified 2D jet must refuse, got {type(result).__name__}"
        )

    def test_g3_certified_at_exact_boundary_passes(self):
        """Φ(u) = u⁴ with certified_order=4 → WH succeeds.

        This is the boundary case: the principal degree equals
        certified_order, so the coefficient IS trusted.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate

        jet = PolynomialJet(dim=1, terms={(4,): 1.0}, certified_order=4)
        datum = self._make_reduced_datum(r=1, p=1, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, WeightedHomogeneousTemplate)

    # ── G4: Multi-dimensional separable WH (r=2) ──────────────────

    def test_g4_2d_separable_wh(self):
        """Φ(u₁,u₂) = u₁⁴ + 2·u₂⁴ with certified_order=4 → WH.

        Separable product form.  Weights = (1/4, 1/4), total = 1/2.
        Principal integral = Π_i ∫ e^{-cᵢ|sᵢ|^4} dsᵢ.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (0, 4): 2.0},
            certified_order=4,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, WeightedHomogeneousTemplate), (
            f"Expected WeightedHomogeneousTemplate, got {type(result).__name__}"
        )
        assert result.weights == (0.25, 0.25)
        assert result.total_weight == pytest.approx(0.5)
        assert result.leading_monomial_degrees == (4, 4)
        assert result.kernel_dim == 2

        # Verify principal integral factors correctly
        from math import lgamma
        log_int_1 = np.log(2) + lgamma(0.25) - np.log(4) - 0.25 * np.log(1.0)
        log_int_2 = np.log(2) + lgamma(0.25) - np.log(4) - 0.25 * np.log(2.0)
        expected_log = log_int_1 + log_int_2
        assert result.log_principal_integral == pytest.approx(expected_log, rel=1e-10)

    def test_g4_2d_with_nonseparable_mixed_refuses_wh(self):
        """Φ(u₁,u₂) = u₁⁴ + u₂⁴ + u₁²u₂² with certified_order=4.

        The mixed term u₁²u₂² has weighted degree = 2·(1/4) + 2·(1/4) = 1.
        This sits ON the Newton boundary, making the principal germ
        non-separable.  The separable WH classifier must refuse.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import (
            UnresolvedKernelJetRefusal,
            WeightedHomogeneousTemplate,
        )

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (0, 4): 1.0, (2, 2): 1.0},
            certified_order=4,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum)

        # WH should refuse due to non-separable principal germ.
        # Without branch_straightened=True, branch is not attempted.
        # Expected: refusal.
        assert not isinstance(result, WeightedHomogeneousTemplate), (
            "Non-separable principal germ u₁⁴ + u₂⁴ + u₁²u₂² must not "
            "be classified as separable WH"
        )

    def test_g4_mixed_degree_wh(self):
        """Φ(u₁,u₂) = u₁⁴ + u₂⁶ with certified_order=6 → WH.

        Different degrees: ω₁=1/4, ω₂=1/6, total = 5/12.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(4, 0): 1.0, (0, 6): 1.0},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=2, jet=jet)
        result = classify_singular_kernel(datum)

        assert isinstance(result, WeightedHomogeneousTemplate)
        assert result.weights == pytest.approx((0.25, 1.0 / 6.0))
        assert result.total_weight == pytest.approx(5.0 / 12.0)


# ═══════════════════════════════════════════════════════════════════════
# Family H: Comparator API tests
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyH:
    """Certified local asymptotic comparator API.

    Tests that learning_exponent, multiplicity, and log_local_evidence
    return the correct values for each template class.
    """

    # ── Regular templates ──────────────────────────────────────────────

    def test_h1_regular_interior_comparator(self):
        """Regular interior: λ = m/2, mult = 1, full formula check."""
        from nomogeo.evidence import local_evidence_from_hessian
        from nomogeo.regime_types import RegularEvidenceTemplate

        H = np.diag([2.0, 3.0, 5.0])
        result = local_evidence_from_hessian(H)
        assert isinstance(result, RegularEvidenceTemplate)

        # λ = 3/2
        assert result.local_learning_exponent == pytest.approx(1.5)
        assert result.local_multiplicity == 1

        # Full formula at n=100, S_vis=1.5:
        # log Z = -100*1.5 - 1.5*log(100) + (3/2)*log(2π) - (1/2)*log(30)
        import math
        n, S_vis = 100, 1.5
        expected = (
            -n * S_vis
            - 1.5 * math.log(n)
            + 1.5 * math.log(2 * math.pi)
            - 0.5 * math.log(2.0 * 3.0 * 5.0)
        )
        assert result.log_local_evidence(n, S_vis) == pytest.approx(expected)

    def test_h2_regular_quotient_comparator(self):
        """Regular quotient: λ = m_⊥/2, orbit volume additive."""
        from nomogeo.evidence import local_evidence_from_hessian
        from nomogeo.regime_types import (
            ConeKind, ConeSpec, OrbitSpec, JacobianConvention,
            RegularEvidenceTemplate,
        )

        H_perp = np.diag([4.0, 9.0])
        orbit = OrbitSpec(
            orbit_dim=1,
            jacobian_convention=JacobianConvention.SLICE_LEBESGUE,
            log_orbit_volume=np.log(2 * np.pi),
        )
        result = local_evidence_from_hessian(H_perp, orbit=orbit)
        assert isinstance(result, RegularEvidenceTemplate)

        # λ = 2/2 = 1 (transverse dim is 2)
        assert result.local_learning_exponent == pytest.approx(1.0)
        assert result.local_multiplicity == 1

        import math
        n, S_vis = 50, 0.8
        expected = (
            -n * S_vis
            - 1.0 * math.log(n)
            + 1.0 * math.log(2 * math.pi)
            - 0.5 * math.log(36.0)
            + math.log(2 * math.pi)  # orbit volume
        )
        assert result.log_local_evidence(n, S_vis) == pytest.approx(expected)

    # ── Singular WH template ──────────────────────────────────────────

    def _make_reduced_datum(self, r, p, jet, *, cone_kind=ConeKind.FULL_SPACE):
        eigs = np.array([2.0] * p)
        log_prefactor = 0.5 * p * np.log(2 * np.pi) - 0.5 * np.sum(np.log(eigs))
        return ReducedKernelActionDatum(
            kernel_dim=r, kernel_action=jet,
            kernel_cone=ConeSpec(kind=cone_kind, ambient_dim=r),
            positive_normal_dim=p, log_normal_prefactor=log_prefactor, total_dim=r + p,
        )

    def test_h3_wh_quartic_comparator(self):
        """WH u⁴: λ = p/2 + 1/4, mult = 1, full formula."""
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import WeightedHomogeneousTemplate
        from nomogeo.evidence import dispatch_singular_evidence
        from nomogeo.regime_types import SingularEvidenceTemplate, RegimeKind

        jet = PolynomialJet(dim=1, terms={(4,): 1.0}, certified_order=4)
        datum = self._make_reduced_datum(r=1, p=2, jet=jet)
        sc = classify_singular_kernel(datum)
        assert isinstance(sc, WeightedHomogeneousTemplate)

        ev = dispatch_singular_evidence(sc)
        assert isinstance(ev, SingularEvidenceTemplate)
        assert ev.regime == RegimeKind.SINGULAR_WH

        # λ = 2/2 + 1/4 = 1.25
        assert ev.local_learning_exponent == pytest.approx(1.25)
        assert ev.local_multiplicity == 1

        import math
        n, S_vis = 200, 0.5
        # log_normal_prefactor = (2/2)*log(2π) - (1/2)*log(4) = log(2π) - log(2)
        # log_principal_integral = log(2*Γ(1/4)/4) = log(2) + lgamma(0.25) - log(4)
        log_np = math.log(2 * math.pi) - math.log(2.0)
        log_pi = math.log(2) + math.lgamma(0.25) - math.log(4)

        expected = (
            -n * S_vis
            - 1.25 * math.log(n)
            + log_np + log_pi
        )
        assert ev.log_local_evidence(n, S_vis) == pytest.approx(expected)

    def test_h4_critical_branch_comparator(self):
        """Critical branch u⁴ + u²v² + v⁶: λ = p/2 + 1/2, mult = 2."""
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import CriticalBranchTemplate
        from nomogeo.evidence import dispatch_singular_evidence
        from nomogeo.regime_types import SingularEvidenceTemplate, RegimeKind

        jet = PolynomialJet(
            dim=2,
            terms={(0, 6): 1.0, (2, 2): 1.0, (4, 0): 1.0},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=3, jet=jet)
        sc = classify_singular_kernel(datum, branch_straightened=True)
        assert isinstance(sc, CriticalBranchTemplate)

        ev = dispatch_singular_evidence(sc)
        assert isinstance(ev, SingularEvidenceTemplate)
        assert ev.regime == RegimeKind.SINGULAR_CRITICAL_BRANCH

        # λ = 3/2 + 1/2 = 2.0
        assert ev.local_learning_exponent == pytest.approx(2.0)
        assert ev.local_multiplicity == 2

        import math
        n, S_vis = 500, 0.3
        log_np = 0.5 * 3 * math.log(2 * math.pi) - 0.5 * 3 * math.log(2.0)
        # κ = (6*(4-2) - 2*4)/(2*4) = (12-8)/8 = 0.5
        # log_coefficient = log(√π/6 · 1.0^{-1/2} · 0.5)
        log_coeff = 0.5 * math.log(math.pi) - math.log(6.0) - 0.5 * math.log(1.0) + math.log(0.5)
        expected = (
            -n * S_vis
            - 2.0 * math.log(n)
            + math.log(math.log(n))  # (m-1)*log(log n) with m=2
            + log_np + log_coeff
        )
        assert ev.log_local_evidence(n, S_vis) == pytest.approx(expected)


# ═══════════════════════════════════════════════════════════════════════
# Family I: Hard-scope condition tests
# ═══════════════════════════════════════════════════════════════════════

class TestFamilyI:
    """Hard scope condition for critical quadratic branch channel."""

    def _make_reduced_datum(self, r, p, jet, *, cone_kind=ConeKind.FULL_SPACE):
        eigs = np.array([2.0] * p)
        log_prefactor = 0.5 * p * np.log(2 * np.pi) - 0.5 * np.sum(np.log(eigs))
        return ReducedKernelActionDatum(
            kernel_dim=r, kernel_action=jet,
            kernel_cone=ConeSpec(kind=cone_kind, ambient_dim=r),
            positive_normal_dim=p, log_normal_prefactor=log_prefactor, total_dim=r + p,
        )

    def test_i1_clean_critical_branch_passes_scope(self):
        """u⁴ + u²v² + v⁶ (only principal terms) → passes hard scope."""
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import CriticalBranchTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(0, 6): 1.0, (2, 2): 1.0, (4, 0): 1.0},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)
        assert isinstance(result, CriticalBranchTemplate)

    def test_i2_stray_monomial_violates_scope(self):
        """u⁴ + u²v² + v⁶ + u²v⁴ violates hard scope → refusal.

        Stray term u²v⁴: weighted order = 2/4 + 4/6 = 0.5 + 0.667 = 1.167.
        κ = 0.5, threshold = 1.5.  Since 1.167 ≤ 1.5, hard scope fails.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import (
            CriticalBranchTemplate,
            UnresolvedKernelJetRefusal,
        )

        jet = PolynomialJet(
            dim=2,
            terms={(0, 6): 1.0, (2, 2): 1.0, (4, 0): 1.0, (2, 4): 0.5},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)
        assert not isinstance(result, CriticalBranchTemplate), (
            "Stray u²v⁴ at weighted order 1.167 ≤ 1+κ=1.5 must trigger "
            "hard scope refusal"
        )

    def test_i3_high_order_stray_passes_scope(self):
        """u⁴ + u²v² + v⁶ + u⁴v⁴ passes hard scope.

        Stray term u⁴v⁴: weighted order = 4/4 + 4/6 = 1.667.
        κ = 0.5, threshold = 1.5.  Since 1.667 > 1.5, hard scope passes.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import CriticalBranchTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(0, 6): 1.0, (2, 2): 1.0, (4, 0): 1.0, (4, 4): 0.5},
            certified_order=8,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)
        assert isinstance(result, CriticalBranchTemplate)

    def test_i4_boundary_stray_violates_scope(self):
        """u⁴ + u²v² + v⁶ + u³v³ violates hard scope.

        Stray term u³v³: weighted order = 3/4 + 3/6 = 0.75 + 0.5 = 1.25.
        κ = 0.5, threshold = 1.5.  Since 1.25 ≤ 1.5, hard scope fails.
        """
        from nomogeo.singular import classify_singular_kernel
        from nomogeo.singular_types import CriticalBranchTemplate

        jet = PolynomialJet(
            dim=2,
            terms={(0, 6): 1.0, (2, 2): 1.0, (4, 0): 1.0, (3, 3): 0.3},
            certified_order=6,
        )
        datum = self._make_reduced_datum(r=2, p=1, jet=jet)
        result = classify_singular_kernel(datum, branch_straightened=True)
        assert not isinstance(result, CriticalBranchTemplate)
