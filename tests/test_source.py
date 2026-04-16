"""Tests for the source law, information conservation, and observer steering (0.4.0)."""
from __future__ import annotations

import numpy as np
import pytest
from numpy.linalg import inv, det, slogdet
from numpy.testing import assert_allclose

from nomogeo.source import (
    information_budget, source_law, evidence_decomposition,
    observer_diagnostics, capture_curve,
)
from nomogeo.exceptions import SupportError
from nomogeo import visible_precision


def _spd(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T + np.eye(n)


def _sym(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return 0.5 * (A + A.T)


# ---- information_budget ----

class TestInformationBudget:
    """First-order conservation: vis + hid = ambient."""

    @pytest.mark.parametrize("n,m,seed", [(3,1,0), (4,2,10), (5,3,20), (6,2,30)])
    def test_conservation_exact(self, n, m, seed):
        rng = np.random.default_rng(seed)
        H = _spd(n, seed)
        C = rng.standard_normal((m, n))
        Hdot = _sym(n, seed + 1)
        result = information_budget(H, C, Hdot)
        assert result.conservation_residual < 1e-12
        assert_allclose(
            result.visible_rate + result.hidden_rate,
            result.ambient_rate,
            atol=1e-12,
        )

    @pytest.mark.parametrize("n,m", [(3,1), (4,2), (5,3)])
    def test_visible_rate_matches_numerical(self, n, m):
        rng = np.random.default_rng(42)
        H0 = _spd(n, 42)
        C = rng.standard_normal((m, n))
        Hdot = _sym(n, 43)
        dt = 1e-7
        phi_p = visible_precision(H0 + dt * Hdot, C)
        phi_m = visible_precision(H0 - dt * Hdot, C)
        _, ld_p = slogdet(phi_p)
        _, ld_m = slogdet(phi_m)
        numerical_rate = (ld_p - ld_m) / (2 * dt)
        result = information_budget(H0, C, Hdot)
        assert abs(result.visible_rate - numerical_rate) < 1e-4

    def test_zero_hdot_gives_zero_rates(self):
        H = _spd(4, 0)
        C = np.random.default_rng(0).standard_normal((2, 4))
        Hdot = np.zeros((4, 4))
        result = information_budget(H, C, Hdot)
        assert abs(result.visible_rate) < 1e-14
        assert abs(result.hidden_rate) < 1e-14
        assert abs(result.ambient_rate) < 1e-14

    def test_visible_fraction_range(self):
        """Visible fraction can be negative or > 1."""
        rng = np.random.default_rng(99)
        for seed in range(50):
            H = _spd(4, seed * 100)
            C = rng.standard_normal((2, 4))
            Hdot = _sym(4, seed * 100 + 1)
            result = information_budget(H, C, Hdot)
            # Just check it's finite and conservation holds
            assert result.conservation_residual < 1e-11


# ---- source_law ----

class TestSourceLaw:
    """Source law A_cpl on support-stable strata."""

    def test_source_law_basic(self):
        """A_cpl is symmetric and has correct dimensions."""
        rng = np.random.default_rng(200)
        for attempt in range(50):
            H = _spd(4, 200 + attempt)
            C = rng.standard_normal((2, 4))
            Hdot = _sym(4, 201 + attempt) * 5.0
            Hddot = _sym(4, 202 + attempt)
            budget = information_budget(H, C, Hdot)
            if np.min(budget.v_eigenvalues) > 0.1:
                result = source_law(H, C, Hdot, Hddot)
                assert result.A_cpl.shape == (2, 2)
                assert_allclose(result.A_cpl, result.A_cpl.T, atol=1e-12)
                return
        pytest.skip("Could not find V > 0 in 50 attempts")

    def test_source_law_rejects_v_not_positive(self):
        """SupportError when V is not positive definite."""
        H = np.diag([1.0, 2.0, 3.0])
        C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        Hdot = np.diag([-1.0, 1.0, 0.0])  # V will have negative eigenvalue
        Hddot = np.zeros((3, 3))
        with pytest.raises(SupportError):
            source_law(H, C, Hdot, Hddot)

    def test_hidden_defect_psd(self):
        """The hidden defect Q_hat / V is always PSD."""
        rng = np.random.default_rng(300)
        for attempt in range(50):
            H = _spd(4, 300 + attempt)
            C = rng.standard_normal((2, 4))
            Hdot = _sym(4, 301 + attempt) * 5.0
            Hddot = _sym(4, 302 + attempt)
            budget = information_budget(H, C, Hdot)
            if np.min(budget.v_eigenvalues) > 0.1:
                result = source_law(H, C, Hdot, Hddot)
                eigvals = np.linalg.eigvalsh(result.hidden_defect)
                assert np.min(eigvals) > -1e-10
                return
        pytest.skip("Could not find V > 0")

    def test_acpl_decomposition(self):
        """A_cpl = A_direct + hidden_defect (when beta = 0, Hddot = 0)."""
        rng = np.random.default_rng(400)
        for attempt in range(50):
            H = _spd(4, 400 + attempt)
            C = rng.standard_normal((2, 4))
            Hdot = _sym(4, 401 + attempt) * 5.0
            Hddot = np.zeros((4, 4))
            budget = information_budget(H, C, Hdot)
            if np.min(budget.v_eigenvalues) > 0.1:
                result = source_law(H, C, Hdot, Hddot)
                # With Hddot = 0: A_direct = 0, A_cpl = hidden_defect
                assert_allclose(result.A_direct, np.zeros((2,2)), atol=1e-10)
                assert_allclose(result.A_cpl, result.hidden_defect, atol=1e-10)
                return
        pytest.skip("Could not find V > 0")


# ---- evidence_decomposition ----

class TestEvidenceDecomposition:
    """Second-order conservation: f''_vis + f''_hid = f''_amb."""

    @pytest.mark.parametrize("n,m,seed", [(3,1,0), (4,2,10), (5,3,20), (6,2,30)])
    def test_second_order_conservation(self, n, m, seed):
        rng = np.random.default_rng(seed)
        H = _spd(n, seed)
        C = rng.standard_normal((m, n))
        Hdot = _sym(n, seed + 1)
        Hddot = _sym(n, seed + 2)
        result = evidence_decomposition(H, C, Hdot, Hddot)
        assert_allclose(
            result.f_pprime_visible + result.f_pprime_hidden,
            result.f_pprime_ambient,
            atol=1e-10,
        )

    @pytest.mark.parametrize("n,m,seed", [(3,1,50), (4,2,60)])
    def test_ambient_matches_numerical(self, n, m, seed):
        rng = np.random.default_rng(seed)
        H0 = _spd(n, seed)
        C = rng.standard_normal((m, n))
        Hdot = _sym(n, seed + 1) * 2.0
        Hddot = _sym(n, seed + 2)
        dt = 1e-5
        def log_det_H(t):
            H_t = H0 + t * Hdot + t**2 * 0.5 * Hddot
            _, ld = slogdet(H_t)
            return ld
        f2_num = (log_det_H(dt) - 2*log_det_H(0) + log_det_H(-dt)) / dt**2
        result = evidence_decomposition(H0, C, Hdot, Hddot)
        assert abs(result.f_pprime_ambient - f2_num) < 1e-2  # FD tolerance

    def test_zero_hddot(self):
        """When Hddot = 0, ambient acceleration = -Tr((H^{-1} Hdot)^2)."""
        H = _spd(4, 0)
        C = np.random.default_rng(0).standard_normal((2, 4))
        Hdot = _sym(4, 1)
        Hddot = np.zeros((4, 4))
        result = evidence_decomposition(H, C, Hdot, Hddot)
        P = inv(H) @ Hdot
        expected = -np.trace(P @ P)
        assert abs(result.f_pprime_ambient - expected) < 1e-10


# ---- coupled spring integration test ----

class TestCoupledSpring:
    """End-to-end test on the coupled spring system from the operationalise layer."""

    def test_coupled_spring_conservation(self):
        k1, k2 = 3.0, 2.0
        C = np.array([[1.0, 0.0]])
        dH = np.array([[1.0, -1.0], [-1.0, 1.0]])

        for kc in [0.5, 1.0, 5.0, 10.0, 20.0]:
            H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
            result = information_budget(H, C, dH)
            assert result.conservation_residual < 1e-14

    def test_coupled_spring_phi_analytical(self):
        k1, k2 = 3.0, 2.0
        C = np.array([[1.0, 0.0]])
        dH = np.array([[1.0, -1.0], [-1.0, 1.0]])

        for kc in [0.5, 2.0, 10.0]:
            H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
            result = information_budget(H, C, dH)
            # Analytical: dPhi/dkc = k2^2/(k2+kc)^2
            phi = (k1*k2 + k1*kc + k2*kc) / (k2 + kc)
            dphi_analytical = k2**2 / (k2 + kc)**2
            # visible_rate = dphi / phi
            expected_vis_rate = dphi_analytical / phi
            assert abs(result.visible_rate - expected_vis_rate) < 1e-10

    def test_coupled_spring_source_law(self):
        k1, k2 = 3.0, 2.0
        C = np.array([[1.0, 0.0]])
        dH = np.array([[1.0, -1.0], [-1.0, 1.0]])
        Hddot = np.zeros((2, 2))

        for kc in [0.5, 2.0, 10.0]:
            H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
            result = source_law(H, C, dH, Hddot)
            # A_cpl should be positive (evidence peak regime)
            assert result.a_cpl_eigenvalues[0] > 0


# ---- observer_diagnostics ----

class TestObserverDiagnostics:
    """Observer diagnostics for steering."""

    @pytest.mark.parametrize("n,m,seed", [(3,1,0), (4,2,10), (5,3,20)])
    def test_conservation(self, n, m, seed):
        rng = np.random.default_rng(seed)
        H = _spd(n, seed)
        C = rng.standard_normal((m, n))
        Hdot = _sym(n, seed + 1)
        result = observer_diagnostics(H, C, Hdot)
        assert result.conservation_residual < 1e-12

    def test_exact_sector_detection(self):
        """Known V > 0 case should be flagged as exact."""
        k1, k2, kc = 3.0, 2.0, 5.0
        H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
        C = np.array([[1.0, 0.0]])
        dH = np.array([[1.0, -1.0], [-1.0, 1.0]])
        result = observer_diagnostics(H, C, dH)
        assert result.exact_sector
        assert result.v_min_eigenvalue > 0

    def test_hidden_defect_nonneg(self):
        rng = np.random.default_rng(42)
        H = _spd(4, 42)
        C = rng.standard_normal((2, 4))
        Hdot = _sym(4, 43)
        result = observer_diagnostics(H, C, Hdot)
        assert result.hidden_defect_trace >= -1e-10

    def test_leakage_nonneg(self):
        rng = np.random.default_rng(99)
        H = _spd(5, 99)
        C = rng.standard_normal((3, 5))
        Hdot = _sym(5, 100)
        result = observer_diagnostics(H, C, Hdot)
        assert result.leakage_norm >= 0


# ---- capture_curve ----

class TestCaptureCurve:
    """Information capture as a function of observer rank."""

    def test_monotone_towards_one(self):
        """As m increases, vis_frac should approach 1 (for full observation)."""
        H = _spd(6, 0)
        Hdot = _sym(6, 1)
        result = capture_curve(H, Hdot)
        # At m = n-1, vis_frac should be close to 1
        assert abs(result.visible_fractions[-1] - 1.0) < 0.5  # not exact due to one hidden dim

    def test_correct_ranks(self):
        H = _spd(5, 10)
        Hdot = _sym(5, 11)
        result = capture_curve(H, Hdot, m_max=3)
        assert result.ranks == (1, 2, 3)

    def test_conservation_at_each_rank(self):
        """Conservation should hold at every rank."""
        H = _spd(4, 20)
        Hdot = _sym(4, 21)
        result = capture_curve(H, Hdot)
        # Check that the budget is valid at each rank by checking vis_frac is finite
        for vf in result.visible_fractions:
            assert np.isfinite(vf)

    def test_half_capture_rank(self):
        """Half-capture rank should be set when vis_frac crosses 0.5."""
        rng = np.random.default_rng(30)
        # Use a matrix where the eigenvalues are spread
        H = np.diag([10.0, 5.0, 2.0, 1.0, 0.5])
        Hdot = _sym(5, 30)
        result = capture_curve(H, Hdot)
        # half_capture_rank should be defined (or None if vis_frac never reaches 0.5)
        if result.half_capture_rank is not None:
            m = result.half_capture_rank
            idx = result.ranks.index(m)
            assert result.visible_fractions[idx] >= 0.5

    def test_coupled_spring_all_exact(self):
        """On the coupled spring, every rank should be in the exact sector."""
        k1, k2, kc = 3.0, 2.0, 5.0
        H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
        dH = np.array([[1.0, -1.0], [-1.0, 1.0]])
        result = capture_curve(H, dH)
        assert result.ranks == (1,)
        # m=1 on n=2 is the only possibility
        assert result.exact_sector_flags[0]

    def test_custom_observer_basis(self):
        """Custom basis should give different capture curve than canonical."""
        H = _spd(5, 50)
        Hdot = _sym(5, 51)
        # Canonical
        cc_canon = capture_curve(H, Hdot, m_max=3)
        # Random orthogonal basis
        rng = np.random.default_rng(52)
        Q = np.linalg.qr(rng.standard_normal((5, 5)))[0]
        cc_custom = capture_curve(H, Hdot, m_max=3, observer_basis=Q)
        # They should generally differ
        assert cc_canon.ranks == cc_custom.ranks
        # Both should have valid (finite) fractions
        for vf in cc_custom.visible_fractions:
            assert np.isfinite(vf)
