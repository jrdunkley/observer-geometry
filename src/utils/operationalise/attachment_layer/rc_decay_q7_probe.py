from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "datasets" / "Q7_CapacitorVoltageDecay.m"
METADATA_PATH = ROOT / "datasets" / "Q7_CapacitorVoltageDecay_metadata.json"
OUTPUT_DIR = ROOT / "outputs" / "empirical_datasets"


def _extract_array(text: str, name: str) -> list[float]:
    match = re.search(rf"{name}\s*=\s*\[([^\]]+)\]", text)
    if not match:
        raise ValueError(f"could not find MATLAB array {name}")
    values = [item.strip() for item in match.group(1).replace(",", " ").split()]
    return [float(item) for item in values]


def _fit_loglinear(time: list[float], voltage: list[float]) -> dict[str, float]:
    logs = [math.log(value) for value in voltage]
    n = len(time)
    mean_t = sum(time) / n
    mean_l = sum(logs) / n
    s_tt = sum((t - mean_t) ** 2 for t in time)
    s_tl = sum((t - mean_t) * (l - mean_l) for t, l in zip(time, logs))
    slope = s_tl / s_tt
    intercept = mean_l - slope * mean_t
    residuals = [l - (intercept + slope * t) for t, l in zip(time, logs)]
    sigma_log = math.sqrt(sum(r * r for r in residuals) / (n - 2))
    return {
        "V0_hat_V": math.exp(intercept),
        "decay_rate_per_s": slope,
        "tau_hat_s": -1.0 / slope,
        "sigma_log_residual": sigma_log,
    }


def _voltage_sse_fixed_v0(time: list[float], voltage: list[float], v0: float, tau: float) -> float:
    return sum((v - v0 * math.exp(-t / tau)) ** 2 for t, v in zip(time, voltage))


def _fit_tau_fixed_v0(time: list[float], voltage: list[float], v0: float) -> dict[str, Any]:
    lo = 1.0e-6
    hi = max(10.0, 10.0 * max(time))
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - golden * (hi - lo)
    d = lo + golden * (hi - lo)
    for _ in range(200):
        if _voltage_sse_fixed_v0(time, voltage, v0, c) < _voltage_sse_fixed_v0(time, voltage, v0, d):
            hi = d
            d = c
            c = hi - golden * (hi - lo)
        else:
            lo = c
            c = d
            d = lo + golden * (hi - lo)
    tau = 0.5 * (lo + hi)
    fitted = [v0 * math.exp(-t / tau) for t in time]
    residuals = [v - f for v, f in zip(voltage, fitted)]
    sse = sum(r * r for r in residuals)
    dof = len(time) - 1
    sigma_v = math.sqrt(sse / dof)
    sensitivities = [v0 * math.exp(-t / tau) * t / (tau * tau) for t in time]
    fisher_tau = sum((s * s) / (sigma_v * sigma_v) for s in sensitivities)
    return {
        "tau_hat_s": tau,
        "sse_V2": sse,
        "sigma_voltage_residual_V": sigma_v,
        "fitted_voltage_V": fitted,
        "residual_voltage_V": residuals,
        "fisher_tau": fisher_tau,
        "local_variance_tau_s2": 1.0 / fisher_tau,
        "max_abs_standardized_residual": max(abs(r) / sigma_v for r in residuals),
    }


def _nonmonotone_points(time: list[float], voltage: list[float]) -> list[dict[str, float]]:
    warnings: list[dict[str, float]] = []
    for previous, current in zip(range(len(voltage) - 1), range(1, len(voltage))):
        if voltage[current] > voltage[previous]:
            warnings.append(
                {
                    "previous_time_s": time[previous],
                    "previous_voltage_V": voltage[previous],
                    "current_time_s": time[current],
                    "current_voltage_V": voltage[current],
                }
            )
    return warnings


def _write_markdown(summary: dict[str, Any], path: Path) -> None:
    fit = summary["fixed_V0_voltage_fit"]
    log_fit = summary["loglinear_comparison"]
    lines = [
        "# Q7 RC Decay Probe v0",
        "",
        f"- dataset: `{summary['dataset']}`",
        f"- route: `{summary['route']}`",
        f"- status: `{summary['status']}`",
        f"- observations: {summary['observation_count']}",
        "",
        "## Admission",
        "",
        summary["admission"],
        "",
        "Refused routes:",
        "",
    ]
    for item in summary["mandatory_refusals"]:
        lines.append(f"- `{item}`")
    lines.extend(
        [
            "",
            "## Readout",
            "",
            f"- fixed V0: `{summary['declared_controls']['V0_V']} V`",
            f"- tau_hat: `{fit['tau_hat_s']:.8g} s`",
            f"- residual sigma: `{fit['sigma_voltage_residual_V']:.8g} V`",
            f"- Fisher precision for tau: `{fit['fisher_tau']:.8g}`",
            f"- local variance for tau: `{fit['local_variance_tau_s2']:.8g} s^2`",
            f"- max abs standardized residual: `{fit['max_abs_standardized_residual']:.8g}`",
            "",
            "## Log-Linear Comparison",
            "",
            f"- V0_hat: `{log_fit['V0_hat_V']:.8g} V`",
            f"- tau_hat: `{log_fit['tau_hat_s']:.8g} s`",
            f"- decay rate: `{log_fit['decay_rate_per_s']:.8g} s^-1`",
            "",
            "## Warnings",
            "",
        ]
    )
    if summary["data_warnings"]:
        for warning in summary["data_warnings"]:
            lines.append(
                "- nonmonotone point: "
                f"{warning['previous_time_s']} s, {warning['previous_voltage_V']} V -> "
                f"{warning['current_time_s']} s, {warning['current_voltage_V']} V"
            )
    else:
        lines.append("None.")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            summary["boundary"],
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    text = DATA_PATH.read_text(encoding="utf-8")
    time = _extract_array(text, "time")
    voltage = _extract_array(text, "voltage")
    if len(time) != len(voltage) or len(time) < 3:
        raise ValueError("time and voltage arrays must have matching length >= 3")
    if any(v <= 0.0 for v in voltage):
        raise ValueError("all voltages must be positive for this probe")
    v0 = float(metadata["measurement_model"]["V0_V"])
    fixed_fit = _fit_tau_fixed_v0(time, voltage, v0)
    log_fit = _fit_loglinear(time, voltage)
    summary = {
        "schema_version": "rc_decay_q7_probe.v0",
        "dataset": "datasets/Q7_CapacitorVoltageDecay.m",
        "route": "observation_equation_fisher",
        "status": "admissible_empirical_probe",
        "source": metadata["source"],
        "admission": "The MATLAB problem data compile into a small RC exponential observation-equation probe after declaring V0 = 100 V, positive tau, and residual-estimated Gaussian voltage noise.",
        "observation_count": len(time),
        "declared_controls": {
            "V0_V": v0,
            "time_unit": "s",
            "voltage_unit": "V",
        },
        "observations": [{"time_s": t, "voltage_V": v} for t, v in zip(time, voltage)],
        "fixed_V0_voltage_fit": fixed_fit,
        "loglinear_comparison": log_fit,
        "data_warnings": _nonmonotone_points(time, voltage),
        "mandatory_refusals": metadata["mandatory_refusals"],
        "boundary": "This is a small textbook regression dataset. It supports a local Fisher precision for tau only under a residual-estimated voltage-noise convention. It does not identify capacitance without a declared resistance, does not expose static storage Hessian geometry, and does not license hidden-load calls.",
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "rc_decay_q7_probe.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_markdown(summary, OUTPUT_DIR / "rc_decay_q7_probe.md")
    print(
        json.dumps(
            {
                "status": summary["status"],
                "tau_hat_s": fixed_fit["tau_hat_s"],
                "fisher_tau": fixed_fit["fisher_tau"],
                "data_warnings": len(summary["data_warnings"]),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
