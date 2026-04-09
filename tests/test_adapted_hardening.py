from __future__ import annotations

from tools.adapted_hardening import run_adapted_hardening


def test_adapted_hardening_smoke() -> None:
    summary = run_adapted_hardening(seed=20260409, commuting_trials=8, bridge_trials=6)

    assert summary.max_commuting_eta < 1e-10
    assert summary.max_commuting_visible_score_gap < 1e-8
    assert summary.max_bridge_eta < 1e-10
    assert summary.max_bridge_scaled_quartic_residual < 1e-8

    assert all(summary.real_bundle_noncommuting_rejections.values())
    assert all(value > 1e-3 for value in summary.real_bundle_commutator_norms.values())
    assert all(value < 1e-10 for value in summary.real_bundle_exact_success_eta.values())
    assert all(value > 0.0 for value in summary.real_bundle_exact_success_score.values())

    assert summary.rotated_condition_acceptance["cond_1e4"] == 12
    assert summary.rotated_condition_acceptance["cond_1e6"] == 12
    assert summary.rotated_condition_acceptance["cond_1e9"] == 0
