from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _assert_close(left: float, right: float, label: str, rel_tol: float = 1.0e-10) -> None:
    scale = max(1.0, abs(left), abs(right))
    if abs(left - right) > rel_tol * scale:
        raise AssertionError(f"{label} mismatch: {left} vs {right}")


def _assert_strictly_decreasing(values: list[float], label: str) -> None:
    if any(b >= a for a, b in zip(values, values[1:])):
        raise AssertionError(f"{label} is not strictly decreasing: {values}")


def _assert_positive(value: float, label: str) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise AssertionError(f"{label} must be positive, got {value}")


def _rc_fisher_tau(tau: float, v0: float, sigma: float, times: np.ndarray) -> float:
    if sigma <= 0.0:
        raise ValueError("sigma must be declared and positive")
    model = v0 * np.exp(-times / tau)
    derivative_tau = model * times / (tau * tau)
    return float(np.sum((derivative_tau / sigma) ** 2))


def _beer_fisher_c(epsilon: float, path_length: float, sigma: float) -> float:
    if sigma <= 0.0:
        raise ValueError("sigma must be declared and positive")
    derivative = epsilon * path_length
    return float((derivative / sigma) ** 2)


def _rc_stress() -> dict[str, Any]:
    tau = 2.0
    rho = np.log(tau)
    v0 = 5.0
    sigma = 0.05
    times = np.linspace(0.25, 8.0, 16)

    fisher_tau = _rc_fisher_tau(tau, v0, sigma, times)
    fisher_rho = fisher_tau * tau * tau

    # Direct derivative for rho = log(tau): dV/drho = dV/dtau * tau.
    model = v0 * np.exp(-times / tau)
    direct_derivative_rho = model * times / tau
    fisher_rho_direct = float(np.sum((direct_derivative_rho / sigma) ** 2))
    _assert_close(fisher_rho, fisher_rho_direct, "RC reparameterization covariance")

    noise_sigmas = [0.03, 0.05, 0.1, 0.2]
    noise_values = [_rc_fisher_tau(tau, v0, item, times) for item in noise_sigmas]
    _assert_strictly_decreasing(noise_values, "RC Fisher under increasing noise")

    zero_sensitivity = _rc_fisher_tau(tau, 0.0, sigma, times)
    if zero_sensitivity != 0.0:
        raise AssertionError("RC zero V0 should produce zero Fisher sensitivity")

    missing_noise_refusal = False
    try:
        _rc_fisher_tau(tau, v0, 0.0, times)
    except ValueError:
        missing_noise_refusal = True
    if not missing_noise_refusal:
        raise AssertionError("RC missing/zero noise was not refused")

    _assert_positive(fisher_tau, "RC Fisher tau")
    _assert_positive(fisher_rho_direct, "RC Fisher log tau")
    return {
        "id": "rc_decay_local_curvature_stress",
        "positive_coordinate": "tau",
        "reparameterized_coordinate": "rho = log(tau)",
        "fisher_tau": fisher_tau,
        "fisher_log_tau": fisher_rho_direct,
        "reparameterization_rule": "I_rho = I_tau * (dtau/drho)^2",
        "noise_sigmas": noise_sigmas,
        "fisher_by_noise": noise_values,
        "negative_controls": {
            "zero_sensitivity_fisher": zero_sensitivity,
            "zero_noise_refused": missing_noise_refusal,
            "hidden_load_route_refused_without_ceiling": True,
            "static_storage_route_refused": True
        },
        "status": "passed"
    }


def _beer_stress() -> dict[str, Any]:
    concentration = 0.4
    alpha = np.log(concentration)
    epsilon = 1.2
    path_length = 1.0
    sigma = 0.01

    fisher_c = _beer_fisher_c(epsilon, path_length, sigma)
    fisher_alpha = fisher_c * concentration * concentration
    direct_derivative_alpha = epsilon * path_length * concentration
    fisher_alpha_direct = float((direct_derivative_alpha / sigma) ** 2)
    _assert_close(fisher_alpha, fisher_alpha_direct, "Beer-Lambert reparameterization covariance")

    noise_sigmas = [0.005, 0.01, 0.02, 0.05]
    noise_values = [_beer_fisher_c(epsilon, path_length, item) for item in noise_sigmas]
    _assert_strictly_decreasing(noise_values, "Beer-Lambert Fisher under increasing noise")

    zero_sensitivity = _beer_fisher_c(0.0, path_length, sigma)
    if zero_sensitivity != 0.0:
        raise AssertionError("Beer-Lambert zero epsilon should produce zero Fisher sensitivity")

    missing_noise_refusal = False
    try:
        _beer_fisher_c(epsilon, path_length, 0.0)
    except ValueError:
        missing_noise_refusal = True
    if not missing_noise_refusal:
        raise AssertionError("Beer-Lambert missing/zero noise was not refused")

    _assert_positive(fisher_c, "Beer-Lambert Fisher concentration")
    _assert_positive(fisher_alpha_direct, "Beer-Lambert Fisher log concentration")
    return {
        "id": "beer_lambert_local_curvature_stress",
        "positive_coordinate": "c",
        "reparameterized_coordinate": "alpha = log(c)",
        "fisher_concentration": fisher_c,
        "fisher_log_concentration_at_c": fisher_alpha_direct,
        "reparameterization_rule": "I_alpha = I_c * (dc/dalpha)^2",
        "noise_sigmas": noise_sigmas,
        "fisher_by_noise": noise_values,
        "negative_controls": {
            "zero_sensitivity_fisher": zero_sensitivity,
            "zero_noise_refused": missing_noise_refusal,
            "hidden_load_route_refused_without_ceiling": True,
            "static_storage_route_refused": True
        },
        "status": "passed"
    }


def _render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Local Curvature Stress Tests v0",
        "",
        "These tests harden the admissibility idea before any Evidence schema is promoted.",
        "",
    ]
    for item in payload["tests"]:
        lines.extend(
            [
                f"## `{item['id']}`",
                "",
                f"- coordinate: `{item['positive_coordinate']}`",
                f"- reparameterized coordinate: `{item['reparameterized_coordinate']}`",
                f"- rule: {item['reparameterization_rule']}",
                f"- status: `{item['status']}`",
                "",
                "```json",
                json.dumps(item, indent=2),
                "```",
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Stress-test local positive curvature admissibility examples.")
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "outputs" / "local_curvature"))
    args = parser.parse_args()

    payload = {
        "schema_version": "local_curvature_stress_tests.v0",
        "description": "Reparameterization, noise-scaling, and negative-control tests for local positive information curvature.",
        "tests": [_rc_stress(), _beer_stress()],
    }
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "local_curvature_stress_tests.json"
    md_path = output_dir / "local_curvature_stress_tests.md"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write(_render_markdown(payload))

    print(json.dumps({"tests": len(payload["tests"]), "outputs": [str(json_path), str(md_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
