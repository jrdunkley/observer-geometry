from __future__ import annotations

from tools.validation_sweep import run_validation_sweep


def test_validation_sweep_smoke() -> None:
    summary = run_validation_sweep(seed=123, dims=[2, 3], draws_per_dim=1)

    assert summary.max_lift_projector_residual < 1e-10
    assert summary.max_h_projector_symmetry_residual < 1e-10
    assert summary.max_tower_law_residual < 1e-10
    assert summary.max_hidden_rank_mismatch == 0
    assert summary.max_transport_two_step_residual < 1e-10
    assert summary.max_divergence_contraction_violation == 0.0
