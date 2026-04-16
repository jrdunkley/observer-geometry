from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DATASET_DIR = ROOT / "datasets" / "edinburgh_data_driven_chemistry_unit09_section3"
OUTPUT_DIR = ROOT / "outputs" / "empirical_datasets"
KNOWN_SAMPLES = {"A", "D", "E", "F"}
UNKNOWN_SAMPLE = "C"


def _read_concentrations(path: Path) -> dict[str, float | None]:
    rows = path.read_text(encoding="utf-8").splitlines()
    header = rows[1].lstrip("#").split(",")
    values = rows[2].split(",")
    concentrations: dict[str, float | None] = {}
    for name, value in zip(header, values):
        stripped = value.strip()
        concentrations[name.strip()] = None if stripped.lower() == "nan" else float(stripped)
    return concentrations


def _read_spectrum(path: Path) -> dict[float, float]:
    spectrum: dict[float, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        header = next(reader, None)
        if header is None or len(header) < 2 or "Wavelength" not in header[0] or "Abs" not in header[1]:
            raise ValueError(f"{path.name}: missing wavelength/absorbance header")
        for row in reader:
            if len(row) < 2:
                continue
            try:
                spectrum[float(row[0])] = float(row[1])
            except ValueError:
                continue
    return spectrum


def _load_records() -> tuple[dict[str, float | None], list[dict[str, Any]]]:
    concentrations = _read_concentrations(DATASET_DIR / "concentrations.csv")
    records: list[dict[str, Any]] = []
    for path in sorted(DATASET_DIR.glob("Rhodamine_6G_in_methanol_*.csv")):
        stem = path.stem
        sample_rep = stem.rsplit("_", 1)[-1]
        sample = sample_rep[0]
        replicate = int(sample_rep[1:])
        records.append(
            {
                "file": path.name,
                "sample": sample,
                "replicate": replicate,
                "concentration_mM": concentrations.get(sample),
                "spectrum": _read_spectrum(path),
            }
        )
    return concentrations, records


def _select_lambda_star(records: list[dict[str, Any]]) -> float:
    sums: dict[float, float] = {}
    counts: dict[float, int] = {}
    for record in records:
        if record["sample"] not in KNOWN_SAMPLES:
            continue
        for wavelength, absorbance in record["spectrum"].items():
            if 450.0 <= wavelength <= 650.0:
                sums[wavelength] = sums.get(wavelength, 0.0) + absorbance
                counts[wavelength] = counts.get(wavelength, 0) + 1
    if not sums:
        raise ValueError("no visible-region spectra available for lambda_star selection")
    return max(sums, key=lambda wavelength: sums[wavelength] / counts[wavelength])


def _linear_fit(xs: list[float], ys: list[float]) -> dict[str, Any]:
    n = len(xs)
    if n < 3:
        raise ValueError("at least three calibration observations are needed")
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    s_xx = sum((x - x_mean) ** 2 for x in xs)
    s_xy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    if s_xx <= 0.0:
        raise ValueError("calibration concentrations are rank deficient")
    beta = s_xy / s_xx
    alpha = y_mean - beta * x_mean
    residuals = [y - (alpha + beta * x) for x, y in zip(xs, ys)]
    sse = sum(r * r for r in residuals)
    dof = n - 2
    sigma = math.sqrt(sse / dof) if dof > 0 else float("nan")
    total = sum((y - y_mean) ** 2 for y in ys)
    r2 = 1.0 - sse / total if total > 0.0 else float("nan")
    fisher = [
        [n / (sigma * sigma), sum(xs) / (sigma * sigma)],
        [sum(xs) / (sigma * sigma), sum(x * x for x in xs) / (sigma * sigma)],
    ]
    det = fisher[0][0] * fisher[1][1] - fisher[0][1] * fisher[1][0]
    covariance = [
        [fisher[1][1] / det, -fisher[0][1] / det],
        [-fisher[1][0] / det, fisher[0][0] / det],
    ]
    return {
        "alpha": alpha,
        "beta_absorbance_per_mM": beta,
        "sigma_absorbance_residual": sigma,
        "r_squared": r2,
        "fisher_alpha_beta": fisher,
        "covariance_alpha_beta": covariance,
        "residuals": residuals,
    }


def _unknown_estimates(records: list[dict[str, Any]], lambda_star: float, fit: dict[str, Any]) -> list[dict[str, float]]:
    estimates: list[dict[str, float]] = []
    beta = fit["beta_absorbance_per_mM"]
    alpha = fit["alpha"]
    covariance = fit["covariance_alpha_beta"]
    sigma2 = fit["sigma_absorbance_residual"] ** 2
    for record in records:
        if record["sample"] != UNKNOWN_SAMPLE:
            continue
        absorbance = record["spectrum"][lambda_star]
        c_hat = (absorbance - alpha) / beta
        # Delta-method variance for c = (A - alpha) / beta.
        d_alpha = -1.0 / beta
        d_beta = -(absorbance - alpha) / (beta * beta)
        parameter_variance = (
            d_alpha * d_alpha * covariance[0][0]
            + 2.0 * d_alpha * d_beta * covariance[0][1]
            + d_beta * d_beta * covariance[1][1]
        )
        observation_variance = sigma2 / (beta * beta)
        total_variance = parameter_variance + observation_variance
        estimates.append(
            {
                "replicate": record["replicate"],
                "absorbance": absorbance,
                "estimated_concentration_mM": c_hat,
                "local_variance_mM2": total_variance,
                "local_precision_inverse_mM2": 1.0 / total_variance if total_variance > 0.0 else float("inf"),
            }
        )
    return estimates


def _write_markdown(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# Beer-Lambert UV-Vis Probe v0",
        "",
        f"- dataset: `{summary['dataset']}`",
        f"- source: {summary['source_repository']}",
        f"- route: `{summary['route']}`",
        f"- status: `{summary['status']}`",
        f"- lambda_star_nm: `{summary['lambda_star_nm']}`",
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
            "## Calibration",
            "",
            f"- known calibration observations: {summary['calibration']['known_observation_count']}",
            f"- alpha: {summary['calibration']['alpha']:.8g}",
            f"- beta_absorbance_per_mM: {summary['calibration']['beta_absorbance_per_mM']:.8g}",
            f"- residual sigma absorbance: {summary['calibration']['sigma_absorbance_residual']:.8g}",
            f"- r_squared: {summary['calibration']['r_squared']:.8g}",
            "",
            "## Unknown Sample C",
            "",
            "| replicate | absorbance | estimated concentration mM | local precision inverse mM2 |",
            "| ---: | ---: | ---: | ---: |",
        ]
    )
    for estimate in summary["unknown_sample_estimates"]:
        lines.append(
            "| {replicate} | {absorbance:.8g} | {estimated_concentration_mM:.8g} | {local_precision_inverse_mM2:.8g} |".format(
                **estimate
            )
        )
    lines.extend(
        [
            "",
            f"Mean estimated concentration: `{summary['unknown_sample_summary']['mean_estimated_concentration_mM']:.8g} mM`.",
            "",
            "## Boundary",
            "",
            summary["boundary"],
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    metadata = json.loads((DATASET_DIR / "metadata.json").read_text(encoding="utf-8"))
    concentrations, records = _load_records()
    lambda_star = _select_lambda_star(records)
    calibration_records = [
        record for record in records if record["sample"] in KNOWN_SAMPLES and record["concentration_mM"] is not None
    ]
    xs = [float(record["concentration_mM"]) for record in calibration_records]
    ys = [record["spectrum"][lambda_star] for record in calibration_records]
    fit = _linear_fit(xs, ys)
    unknown_estimates = _unknown_estimates(records, lambda_star, fit)
    known_span = [min(xs), max(xs)]
    unknown_values = [item["estimated_concentration_mM"] for item in unknown_estimates]
    inside_span = all(known_span[0] <= value <= known_span[1] for value in unknown_values)
    summary = {
        "schema_version": "beer_lambert_uvvis_probe.v0",
        "dataset": str(DATASET_DIR.relative_to(ROOT)),
        "source_repository": metadata["source"]["repository"],
        "source_directory": metadata["source"]["source_directory"],
        "route": "observation_equation_fisher",
        "status": "admissible_empirical_probe" if inside_span else "refused_unknown_outside_calibration_span",
        "admission": "Spectra and concentrations compile into a Beer-Lambert calibration/analysis-function probe after declaring a common wavelength, known standards, residual noise model, and unknown-sample inversion.",
        "lambda_star_nm": lambda_star,
        "known_concentration_span_mM": known_span,
        "calibration": {
            "known_observation_count": len(xs),
            "alpha": fit["alpha"],
            "beta_absorbance_per_mM": fit["beta_absorbance_per_mM"],
            "sigma_absorbance_residual": fit["sigma_absorbance_residual"],
            "r_squared": fit["r_squared"],
            "fisher_alpha_beta": fit["fisher_alpha_beta"],
            "covariance_alpha_beta": fit["covariance_alpha_beta"],
        },
        "unknown_sample_estimates": unknown_estimates,
        "unknown_sample_summary": {
            "mean_estimated_concentration_mM": sum(unknown_values) / len(unknown_values),
            "min_estimated_concentration_mM": min(unknown_values),
            "max_estimated_concentration_mM": max(unknown_values),
            "inside_known_calibration_span": inside_span,
        },
        "mandatory_refusals": metadata["mandatory_refusals"],
        "boundary": "This is a sensor-law information object and analysis-function probe. It is not a thermodynamic concentration model, not a storage Hessian, and not a hidden-load object. Molar absorptivity is not claimed because path length is not declared in the downloaded data.",
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "beer_lambert_uvvis_probe.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_markdown(summary, OUTPUT_DIR / "beer_lambert_uvvis_probe.md")
    print(json.dumps({k: summary[k] for k in ("status", "lambda_star_nm", "known_concentration_span_mM")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
