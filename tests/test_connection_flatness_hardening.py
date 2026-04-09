from __future__ import annotations

from tools.connection_flatness_hardening import run_connection_flatness_hardening


def test_connection_flatness_hardening_smoke() -> None:
    summary = run_connection_flatness_hardening(seed=20260409, synthetic_trials=12)

    assert summary.max_factorisation_residual < 1e-10
    assert summary.max_phi_recovery_residual < 1e-9
    assert summary.max_metric_split_residual < 1e-9
    assert summary.max_q_from_dk_residual < 1e-9
    assert summary.max_transition_residual < 1e-10
    assert summary.max_hidden_gauge_residual < 1e-8
    assert summary.max_visible_gauge_residual < 1e-10
    assert summary.max_current_variation < 1e-7
    assert summary.max_q_from_current_residual < 1e-7
    assert summary.max_geodesic_deviation_t3_scaled < 5e-2
    assert summary.max_horizontal_projection_residual < 1e-10
    assert summary.max_fixed_k_free_sector_q < 1e-12
    assert summary.min_fixed_k_nonhorizontal_cross_norm > 1e-3

    assert summary.conditioning_acceptance["cond_1e4"] > 0
    assert summary.conditioning_acceptance["cond_1e6"] > 0
    assert summary.conditioning_acceptance["cond_1e8"] > 0

    assert all(value < 1e-9 for value in summary.empirical_transition_residuals.values())
    assert all(value < 1e-6 for value in summary.empirical_current_variation.values())
    assert all(value < 1e-6 for value in summary.empirical_q_from_current_residual.values())
