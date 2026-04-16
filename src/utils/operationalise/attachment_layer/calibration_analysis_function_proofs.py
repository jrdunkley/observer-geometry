from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "calibration_analysis"


def _residual_gate(residuals: np.ndarray, sigma: float, *, parameter_count: int) -> dict[str, Any]:
    standardized = residuals / sigma
    return {
        "max_abs_standardized_residual": float(np.max(np.abs(standardized))),
        "sum_squared_standardized_residuals": float(np.sum(standardized**2)),
        "degrees_of_freedom": int(residuals.size - parameter_count),
        "passed_simple_gate": bool(np.max(np.abs(standardized)) <= 2.5),
        "gate_note": "Simple deterministic proof gate; a real calibration workflow should use residual plots and model checks.",
    }


def load_cell_linear_analysis_function() -> dict[str, Any]:
    forces = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])
    true_alpha = 0.2
    true_beta = 0.05
    sigma = 0.02
    residuals = sigma * np.array([0.12, -0.18, 0.06, 0.10, -0.15, 0.04])
    indications = true_alpha + true_beta * forces + residuals

    design = np.column_stack([np.ones_like(forces), forces])
    fit = np.linalg.lstsq(design, indications, rcond=None)[0]
    alpha_hat = float(fit[0])
    beta_hat = float(fit[1])
    fitted = design @ fit
    fit_residuals = indications - fitted
    fisher = (design.T @ design) / (sigma * sigma)
    covariance = np.linalg.inv(fisher)

    indication_new = 1.45
    validity_range = [float(np.min(indications)), float(np.max(indications))]
    if not (validity_range[0] <= indication_new <= validity_range[1]):
        raise AssertionError("internal proof indication should be inside validity range")
    force_estimate = (indication_new - alpha_hat) / beta_hat
    gradient = np.array([-1.0 / beta_hat, -(indication_new - alpha_hat) / (beta_hat * beta_hat)])
    parameter_variance = float(gradient @ covariance @ gradient)
    indication_variance = float((1.0 / beta_hat) ** 2 * sigma * sigma)
    total_variance = parameter_variance + indication_variance

    return {
        "id": "load_cell_linear_analysis_function",
        "route": "calibration_analysis_function",
        "status": "admissible_analysis_function_record",
        "calibration_function": "D = alpha + beta F",
        "analysis_function": "F = (D - alpha) / beta",
        "standards": {
            "force_values": [float(item) for item in forces],
            "indications": [float(item) for item in indications],
            "sigma_indication": sigma,
        },
        "fit": {
            "alpha_hat": alpha_hat,
            "beta_hat": beta_hat,
            "parameter_covariance": covariance.tolist(),
            "residual_gate": _residual_gate(fit_residuals, sigma, parameter_count=2),
        },
        "validity_range": {
            "indication_min": validity_range[0],
            "indication_max": validity_range[1],
            "rule": "analysis function use is refused outside the calibrated indication span unless a separate theoretical extrapolation argument is declared",
        },
        "analysis_readout": {
            "new_indication": indication_new,
            "force_estimate": float(force_estimate),
            "parameter_variance": parameter_variance,
            "indication_variance": indication_variance,
            "total_local_variance": total_variance,
            "local_standard_uncertainty": float(np.sqrt(total_variance)),
            "local_precision": float(1.0 / total_variance),
        },
        "allowed_route": "analysis_function_record",
        "refused_routes": ["static storage Hessian", "hidden_load without a ceiling"],
        "boundary": "This proves inversion of a declared linear calibration inside its calibration range. It does not make the load cell a storage Hessian or hidden-load object.",
    }


def negative_controls() -> list[dict[str, Any]]:
    sigma = 0.02
    repeated_force = np.array([10.0, 10.0, 10.0, 10.0])
    rank_deficient_design = np.column_stack([np.ones_like(repeated_force), repeated_force])
    rank_deficient_fisher = (rank_deficient_design.T @ rank_deficient_design) / (sigma * sigma)
    eigvals = np.linalg.eigvalsh(rank_deficient_fisher)
    return [
        {
            "id": "load_cell_analysis_extrapolation_refusal",
            "route": "calibration_analysis_function",
            "status": "refused",
            "reason": "A new indication outside the calibrated indication span cannot use the empirical analysis function without a separate extrapolation argument.",
            "calibrated_indication_span": [0.2024, 2.7008],
            "requested_indication": 3.2,
            "refused_route": "analysis_function_record",
        },
        {
            "id": "load_cell_analysis_rank_deficient_refusal",
            "route": "calibration_analysis_function",
            "status": "refused",
            "reason": "Repeated use of one force standard cannot identify both intercept and slope for the analysis function.",
            "fisher_eigenvalues": [float(item) for item in eigvals],
            "refused_route": "analysis_function_record",
        },
    ]


def _markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Calibration Analysis Function Proofs v0",
        "",
        "Status: generated proof packet for calibration-analysis functions.",
        "",
        "These are not attachment cards. They test when a calibrated sensor can expose an analysis function with local uncertainty.",
        "",
    ]
    for proof in payload["proofs"]:
        readout = proof["analysis_readout"]
        lines.extend(
            [
                f"## {proof['id']}",
                "",
                f"- route: `{proof['route']}`",
                f"- status: `{proof['status']}`",
                f"- calibration function: `{proof['calibration_function']}`",
                f"- analysis function: `{proof['analysis_function']}`",
                f"- force estimate: `{readout['force_estimate']:.12g}`",
                f"- local standard uncertainty: `{readout['local_standard_uncertainty']:.12g}`",
                f"- local precision: `{readout['local_precision']:.12g}`",
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
    payload = {
        "schema_version": "calibration_analysis_function_proofs.v0",
        "purpose": "Proof packet for calibration-analysis function admissibility before a measurement attachment schema is promoted.",
        "proofs": [load_cell_linear_analysis_function()],
        "negative_controls": negative_controls(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "calibration_analysis_function_proofs.json"
    md_path = OUT_DIR / "calibration_analysis_function_proofs.md"
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
