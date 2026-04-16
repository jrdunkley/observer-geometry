from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "observation_equation"


def _positive(value: float, label: str) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise AssertionError(f"{label} must be finite and positive, got {value}")


def _spd(matrix: np.ndarray, label: str) -> dict[str, Any]:
    eigvals = np.linalg.eigvalsh(matrix)
    if np.any(eigvals <= 0.0):
        raise AssertionError(f"{label} must be SPD, eigenvalues={eigvals}")
    return {
        "eigenvalues": [float(item) for item in eigvals],
        "condition_number": float(eigvals[-1] / eigvals[0]),
    }


def _residual_gate(residuals: np.ndarray, sigma: float, *, parameter_count: int) -> dict[str, Any]:
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    standardized = residuals / sigma
    dof = int(residuals.size - parameter_count)
    return {
        "max_abs_standardized_residual": float(np.max(np.abs(standardized))),
        "sum_squared_standardized_residuals": float(np.sum(standardized**2)),
        "degrees_of_freedom": dof,
        "passed_simple_gate": bool(np.max(np.abs(standardized)) <= 2.5),
        "gate_note": "Simple deterministic gate for this proof packet; not a substitute for a full residual diagnostic workflow.",
    }


def rc_decay_tau_observation_equation() -> dict[str, Any]:
    tau = 2.0
    v0 = 5.0
    sigma = 0.05
    times = np.linspace(0.25, 8.0, 16)
    model = v0 * np.exp(-times / tau)
    residuals = sigma * np.array([
        0.18,
        -0.22,
        0.12,
        0.07,
        -0.15,
        0.19,
        -0.09,
        0.04,
        -0.11,
        0.14,
        -0.06,
        0.02,
        0.10,
        -0.08,
        0.05,
        -0.03,
    ])
    observations = model + residuals
    derivative_tau = model * times / (tau * tau)
    fisher_tau = float(np.sum((derivative_tau / sigma) ** 2))
    _positive(fisher_tau, "RC decay Fisher precision")
    return {
        "id": "rc_decay_tau_observation_equation",
        "route": "observation_equation_fisher",
        "status": "admissible_local_precision",
        "observation_equation": "V_i = V0 exp(-t_i / tau) + E_i, E_i ~ N(0, sigma^2)",
        "measurand_or_parameter": "tau",
        "controls": {
            "V0": v0,
            "sigma": sigma,
            "times": [float(item) for item in times],
        },
        "observations": [float(item) for item in observations],
        "fisher_object": {
            "kind": "scalar Fisher precision",
            "coordinate": "tau",
            "value": fisher_tau,
            "local_variance": 1.0 / fisher_tau,
        },
        "residual_gate": _residual_gate(residuals, sigma, parameter_count=1),
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["static storage Hessian", "hidden_load without a ceiling"],
        "boundary": "This compiles a sampled voltage trace into local information about tau; it is not capacitor storage or a hidden-load object.",
    }


def repeated_current_mean_observation_equation() -> dict[str, Any]:
    observations = np.array([19.663, 19.639, 19.640, 19.685, 19.678])
    estimate = float(np.mean(observations))
    sigma_hat = float(np.std(observations, ddof=1))
    fisher_mean = float(observations.size / (sigma_hat * sigma_hat))
    _positive(fisher_mean, "current mean Fisher precision")
    residuals = observations - estimate
    return {
        "id": "repeated_current_mean_observation_equation",
        "route": "observation_equation_fisher",
        "status": "admissible_local_precision_with_estimated_sigma",
        "observation_equation": "I_j = I0 + E_j, E_j ~ N(0, sigma^2)",
        "measurand_or_parameter": "I0",
        "observations": [float(item) for item in observations],
        "estimate": estimate,
        "estimated_sigma": sigma_hat,
        "degrees_of_freedom_for_sigma": int(observations.size - 1),
        "fisher_object": {
            "kind": "scalar Fisher precision",
            "coordinate": "I0",
            "value": fisher_mean,
            "local_variance": 1.0 / fisher_mean,
        },
        "residual_gate": _residual_gate(residuals, sigma_hat, parameter_count=1),
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["static storage Hessian", "hidden_load without a ceiling"],
        "boundary": "This is a statistical observation model for a repeated measurement, not a physical storage law.",
    }


def load_cell_linear_calibration_observation_equation() -> dict[str, Any]:
    forces = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])
    alpha = 0.2
    beta = 0.05
    sigma = 0.02
    residuals = sigma * np.array([0.12, -0.18, 0.06, 0.10, -0.15, 0.04])
    indications = alpha + beta * forces + residuals

    design = np.column_stack([np.ones_like(forces), forces])
    fisher = (design.T @ design) / (sigma * sigma)
    spd = _spd(fisher, "load-cell calibration Fisher matrix")
    covariance = np.linalg.inv(fisher)
    beta_hat = np.linalg.lstsq(design, indications, rcond=None)[0]
    fitted = design @ beta_hat
    fit_residuals = indications - fitted
    return {
        "id": "load_cell_linear_calibration_observation_equation",
        "route": "observation_equation_fisher",
        "status": "admissible_local_precision_matrix",
        "observation_equation": "D_i = alpha + beta F_i + E_i, E_i ~ N(0, sigma^2)",
        "measurand_or_parameter": ["alpha", "beta"],
        "controls": {
            "force_values": [float(item) for item in forces],
            "sigma": sigma,
        },
        "observations": [float(item) for item in indications],
        "fit_estimate": {
            "alpha_hat": float(beta_hat[0]),
            "beta_hat": float(beta_hat[1]),
        },
        "fisher_object": {
            "kind": "2x2 Fisher precision matrix",
            "coordinates": ["alpha", "beta"],
            "matrix": fisher.tolist(),
            "covariance": covariance.tolist(),
            "spd_check": spd,
        },
        "residual_gate": _residual_gate(fit_residuals, sigma, parameter_count=2),
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["static storage Hessian", "hidden_load without a ceiling", "analysis-function inversion without validity range"],
        "boundary": "This proves calibration-parameter precision for a declared linear response; using it as an inverse force-measurement function is a later analysis-function step.",
    }


def negative_controls() -> list[dict[str, Any]]:
    sigma = 0.05
    tau = 2.0
    v0 = 5.0
    zero_times = np.zeros(6)
    zero_model = v0 * np.exp(-zero_times / tau)
    zero_derivative = zero_model * zero_times / (tau * tau)
    zero_fisher = float(np.sum((zero_derivative / sigma) ** 2))

    repeated_force = np.array([10.0, 10.0, 10.0, 10.0])
    rank_deficient_design = np.column_stack([np.ones_like(repeated_force), repeated_force])
    rank_deficient_fisher = (rank_deficient_design.T @ rank_deficient_design) / (sigma * sigma)
    rank_deficient_eigvals = np.linalg.eigvalsh(rank_deficient_fisher)

    return [
        {
            "id": "rc_decay_zero_sensitivity_grid_refusal",
            "route": "observation_equation_fisher",
            "status": "refused",
            "reason": "All samples are at t = 0, so dV/dtau is zero and no local precision for tau is exposed.",
            "computed_fisher": zero_fisher,
            "refused_route": "local_positive_curvature",
        },
        {
            "id": "load_cell_rank_deficient_design_refusal",
            "route": "observation_equation_fisher",
            "status": "refused",
            "reason": "All force standards have the same force value, so intercept and slope are not jointly identifiable.",
            "fisher_eigenvalues": [float(item) for item in rank_deficient_eigvals],
            "refused_route": "local_positive_curvature_matrix",
        },
    ]


def _markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Observation Equation Fisher Proofs v0",
        "",
        "Status: generated proof packet for observation-equation local Fisher precision.",
        "",
        "These are not attachment cards. They test when declared observations plus a noise model expose local positive curvature.",
        "",
    ]
    for proof in payload["proofs"]:
        fisher = proof["fisher_object"]
        lines.extend(
            [
                f"## {proof['id']}",
                "",
                f"- route: `{proof['route']}`",
                f"- status: `{proof['status']}`",
                f"- equation: `{proof['observation_equation']}`",
                f"- object: `{fisher['kind']}`",
            ]
        )
        if "value" in fisher:
            lines.append(f"- local precision: `{fisher['value']:.12g}`")
        else:
            lines.append(f"- SPD eigenvalues: `{fisher['spd_check']['eigenvalues']}`")
        lines.extend(
            [
                f"- residual gate passed: `{proof['residual_gate']['passed_simple_gate']}`",
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
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    proofs = [
        rc_decay_tau_observation_equation(),
        repeated_current_mean_observation_equation(),
        load_cell_linear_calibration_observation_equation(),
    ]
    payload = {
        "schema_version": "observation_equation_fisher_proofs.v0",
        "purpose": "Proof packet for observation-equation Fisher curvature before any generic evidence schema is promoted.",
        "proofs": proofs,
        "negative_controls": negative_controls(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "observation_equation_fisher_proofs.json"
    md_path = OUT_DIR / "observation_equation_fisher_proofs.md"
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
