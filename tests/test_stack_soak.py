from __future__ import annotations

from tools.stack_soak import run_stack_soak


def test_stack_soak_smoke() -> None:
    summary = run_stack_soak(
        seed=1234,
        kernel_count=8,
        descent_count=6,
        evidence_count=6,
        example_count=6,
        max_workers=2,
    )

    assert summary.backend_agreement["serial_vs_thread"]["all_match"] is True
    assert summary.backend_agreement["serial_vs_process"]["all_match"] is True
    assert summary.failure_propagation["serial"]["index"] == 1
    assert summary.failure_propagation["thread"]["index"] == 1
    assert summary.failure_propagation["process"]["index"] == 1
    assert summary.forced_process_fallback_matches_serial is True
    assert summary.validation_sweep["max_divergence_contraction_violation"] == 0.0
