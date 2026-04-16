"""
Regression test for Strike 1B: regular-Gaussian benchmark.

Asserts the deterministic sign pattern (seed=1, n=60) and verifies
that the comparator API on RegularEvidenceTemplate reproduces the
exact local evidence score from the benchmark.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from nomogeo.regime_types import RegularEvidenceTemplate, RegimeKind

# Import benchmark machinery
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "benchmarks"))
from strike1b_regular_gaussian import run_benchmark, run_sweep


class TestStrike1BDeterministic:
    """Deterministic flagship case: seed=1, n=60."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.m1, self.m2 = run_benchmark(n=60, seed=1, rho=0.98)
        self.n = 60

    def test_sign_pattern_aic_prefers_m1(self):
        """AIC (log L - k) should prefer the smaller model M1.

        AIC = -2 log L + 2k.  Maximisation score = -AIC/2 = log L - k.
        With Δ log L ≈ +0.05 and Δk = 1, the AIC gap is ≈ -0.95,
        so AIC prefers M1.
        """
        assert self.m1.aic_score > self.m2.aic_score

    def test_sign_pattern_bic_prefers_m1(self):
        """BIC should prefer the smaller model M1."""
        assert self.m1.bic_score > self.m2.bic_score

    def test_sign_pattern_exact_prefers_m2(self):
        """Exact local Laplace evidence (with determinant correction)
        should prefer M2 — the true generating model."""
        assert self.m2.exact_score > self.m1.exact_score

    def test_both_surrogates_reversed(self):
        """Both AIC and BIC pick M1, but exact picks M2.

        This is the core demonstration: the omitted Laplace determinant
        correction reverses both count-penalty surrogates.
        """
        aic_prefers_m1 = self.m1.aic_score > self.m2.aic_score
        bic_prefers_m1 = self.m1.bic_score > self.m2.bic_score
        exact_prefers_m2 = self.m2.exact_score > self.m1.exact_score
        assert aic_prefers_m1 and bic_prefers_m1 and exact_prefers_m2

    def test_determinant_correction_positive(self):
        """The determinant correction (favouring M2) should be positive
        and large enough to flip BIC's ranking."""
        det_corr = (
            -0.5 * self.m2.log_det_J + 0.5 * self.m1.log_det_J
            + (self.m2.k / 2.0) * np.log(2 * np.pi)
            - (self.m1.k / 2.0) * np.log(2 * np.pi)
        )
        bic_gap = self.m2.bic_score - self.m1.bic_score  # negative (BIC prefers M1)
        assert det_corr > 0, "Determinant correction should favour M2"
        assert det_corr > abs(bic_gap), "Correction must be large enough to flip BIC"

    def test_comparator_api_matches_benchmark(self):
        """The comparator API on RegularEvidenceTemplate reproduces
        the exact local evidence score from the benchmark."""
        for m in (self.m1, self.m2):
            # Build a RegularEvidenceTemplate from the benchmark data.
            # log_quadratic_det_term = (k/2)*log(2π) - (1/2)*log det(J)
            log_quad = (m.k / 2.0) * math.log(2 * math.pi) - 0.5 * m.log_det_J
            tmpl = RegularEvidenceTemplate(
                regime=RegimeKind.REGULAR_INTERIOR,
                active_dim=m.k,
                log_quadratic_det_term=log_quad,
            )
            assert tmpl.local_learning_exponent == pytest.approx(m.k / 2.0)
            assert tmpl.local_multiplicity == 1

            computed = tmpl.log_local_evidence(self.n, m.S_vis)
            assert computed == pytest.approx(m.exact_score, abs=1e-8)


class TestStrike1BSweep:
    """Sweep over 200 seeds: structural properties of the reversal."""

    @pytest.fixture(autouse=True, scope="class")
    def setup(self, request):
        request.cls.results = run_sweep(n_seeds=200, n=60, rho=0.98)

    def test_aic_reversal_rate_above_70_percent(self):
        """AIC→Exact reversal should occur in the majority of seeds."""
        reversals = sum(
            1 for m1, m2 in self.results
            if m1.aic_score > m2.aic_score and m2.exact_score > m1.exact_score
        )
        rate = reversals / len(self.results)
        assert rate > 0.70, f"AIC→Exact reversal rate {rate:.1%} is too low"

    def test_bic_reversal_rate_above_80_percent(self):
        """BIC→Exact reversal should occur in the large majority of seeds."""
        reversals = sum(
            1 for m1, m2 in self.results
            if m1.bic_score > m2.bic_score and m2.exact_score > m1.exact_score
        )
        rate = reversals / len(self.results)
        assert rate > 0.80, f"BIC→Exact reversal rate {rate:.1%} is too low"

    def test_exact_prefers_m2_majority(self):
        """Exact local evidence should prefer M2 in the large majority."""
        exact_m2 = sum(1 for m1, m2 in self.results if m2.exact_score > m1.exact_score)
        rate = exact_m2 / len(self.results)
        assert rate > 0.90, f"Exact prefers M2 only {rate:.1%}"

    def test_det_correction_consistently_positive(self):
        """The determinant correction should be positive in every seed."""
        for m1, m2 in self.results:
            det_corr = (
                -0.5 * m2.log_det_J + 0.5 * m1.log_det_J
                + (m2.k / 2.0) * np.log(2 * np.pi)
                - (m1.k / 2.0) * np.log(2 * np.pi)
            )
            assert det_corr > 0
