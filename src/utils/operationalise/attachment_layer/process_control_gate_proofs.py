from __future__ import annotations

import json
from pathlib import Path
from statistics import mean
from typing import Any


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "process_control"


def _z_scores(values: list[float], baseline_mean: float, sigma: float) -> list[float]:
    if sigma <= 0.0:
        raise ValueError("sigma must be positive for a process-control gate")
    return [(value - baseline_mean) / sigma for value in values]


def _gate_readout(values: list[float], baseline_mean: float, sigma: float, action_limit_z: float) -> dict[str, Any]:
    z_scores = _z_scores(values, baseline_mean, sigma)
    max_abs_z = max(abs(item) for item in z_scores)
    return {
        "observation_count": len(values),
        "observed_mean": mean(values),
        "baseline_mean": baseline_mean,
        "repeatability_sigma": sigma,
        "z_scores": z_scores,
        "max_abs_z": max_abs_z,
        "action_limit_z": action_limit_z,
        "passed_action_limit": max_abs_z <= action_limit_z,
    }


def check_standard_process_gate() -> dict[str, Any]:
    baseline_mean = 10.000
    sigma = 0.010
    action_limit_z = 3.0
    observations = [10.004, 9.994, 10.008, 10.001, 9.998, 10.011, 10.000, 9.993]
    readout = _gate_readout(observations, baseline_mean, sigma, action_limit_z)
    if not readout["passed_action_limit"]:
        raise AssertionError("internal process-control proof should pass")
    return {
        "id": "check_standard_process_gate",
        "route": "process_control_gate",
        "status": "gate_passed_no_module_call",
        "check_standard": {
            "quantity": "declared check-standard indication",
            "nominal_value": 10.0,
            "unit": "declared_indication_unit",
        },
        "declarations": {
            "historical_baseline_mean": baseline_mean,
            "repeatability_sigma": sigma,
            "action_limit_z": action_limit_z,
            "control_rule": "refuse the measurement route if any check-standard observation exceeds the three-sigma action limit",
        },
        "observations": observations,
        "gate_readout": readout,
        "emitted_object": "measurement-process validity decision only",
        "allowed_route": "gate_or_refusal_only",
        "theorem_local_calls_licensed": False,
        "refused_routes": ["visible_precision", "hidden_load", "static storage Hessian", "Fisher precision", "theorem-local kernel call requested"],
        "boundary": "This gate says the measurement process is not obviously out of control under the declared check-standard rule. It does not create a module Hessian, Fisher object, or hidden-load object.",
    }


def negative_controls() -> list[dict[str, Any]]:
    baseline_mean = 10.000
    sigma = 0.010
    action_limit_z = 3.0
    drift_values = [10.012, 10.018, 10.026, 10.035, 10.041]
    drift_readout = _gate_readout(drift_values, baseline_mean, sigma, action_limit_z)
    if drift_readout["passed_action_limit"]:
        raise AssertionError("internal drift control should refuse")
    return [
        {
            "id": "check_standard_drift_refusal",
            "route": "process_control_gate",
            "status": "refused",
            "reason": "At least one check-standard observation exceeds the declared three-sigma action limit.",
            "observations": drift_values,
            "gate_readout": drift_readout,
            "refused_route": "gate_or_refusal_only",
            "theorem_local_calls_licensed": False,
        },
        {
            "id": "check_standard_missing_baseline_refusal",
            "route": "process_control_gate",
            "status": "refused",
            "reason": "A process-control gate cannot be evaluated without a stable check standard, historical baseline, and variability declaration.",
            "missing_declarations": ["historical_baseline_mean", "repeatability_sigma", "bias_or_drift_rule"],
            "refused_route": "gate_or_refusal_only",
            "theorem_local_calls_licensed": False,
        },
    ]


def _markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Process-Control Gate Proofs v0",
        "",
        "Status: generated proof packet for measurement-process gates.",
        "",
        "These are not attachment cards. They test whether a measurement route may proceed at all. A passing gate emits a validity decision only; it does not license theorem-local kernel calls.",
        "",
    ]
    for proof in payload["proofs"]:
        readout = proof["gate_readout"]
        lines.extend(
            [
                f"## {proof['id']}",
                "",
                f"- route: `{proof['route']}`",
                f"- status: `{proof['status']}`",
                f"- emitted object: {proof['emitted_object']}",
                f"- observed mean: `{readout['observed_mean']:.12g}`",
                f"- max absolute z: `{readout['max_abs_z']:.12g}`",
                f"- action limit z: `{readout['action_limit_z']:.12g}`",
                f"- theorem-local calls licensed: `{str(proof['theorem_local_calls_licensed']).lower()}`",
                f"- boundary: {proof['boundary']}",
                "",
            ]
        )
    lines.append("## Negative Controls")
    lines.append("")
    for control in payload["negative_controls"]:
        lines.extend(
            [
                f"### {control['id']}",
                "",
                f"- status: `{control['status']}`",
                f"- reason: {control['reason']}",
                f"- refused route: `{control['refused_route']}`",
                f"- theorem-local calls licensed: `{str(control['theorem_local_calls_licensed']).lower()}`",
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    payload = {
        "schema_version": "process_control_gate_proofs.v0",
        "purpose": "Proof packet for process-control gates before empirical measurement objects are trusted as live attachment candidates.",
        "proofs": [check_standard_process_gate()],
        "negative_controls": negative_controls(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "process_control_gate_proofs.json"
    md_path = OUT_DIR / "process_control_gate_proofs.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "proofs": len(payload["proofs"]),
                "negative_controls": len(payload["negative_controls"]),
                "outputs": [str(json_path), str(md_path)],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
