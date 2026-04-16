from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "measurement_equation"


def _central_sensitivity(
    fn: Callable[[dict[str, float]], float],
    values: dict[str, float],
    name: str,
    standard_uncertainty: float,
) -> float:
    if standard_uncertainty <= 0.0:
        raise ValueError(f"{name}: standard_uncertainty must be positive")
    hi = dict(values)
    lo = dict(values)
    hi[name] += standard_uncertainty
    lo[name] -= standard_uncertainty
    return (fn(hi) - fn(lo)) / (2.0 * standard_uncertainty)


def _linearised_case(
    *,
    case_id: str,
    measurement_equation: str,
    values: dict[str, float],
    uncertainties: dict[str, float],
    units: dict[str, str],
    fn: Callable[[dict[str, float]], float],
    analytic_sensitivities: dict[str, float],
    positivity_checks: dict[str, str],
    boundary: str,
) -> dict[str, Any]:
    estimate = fn(values)
    numeric_sensitivities = {
        name: _central_sensitivity(fn, values, name, uncertainty)
        for name, uncertainty in uncertainties.items()
    }
    variance_terms = {
        name: (analytic_sensitivities[name] * uncertainty) ** 2
        for name, uncertainty in uncertainties.items()
    }
    variance = sum(variance_terms.values())
    max_sensitivity_error = max(
        abs(analytic_sensitivities[name] - numeric_sensitivities[name])
        for name in uncertainties
    )
    return {
        "id": case_id,
        "route": "measurement_equation_local_linearisation",
        "status": "admissible_local_precision",
        "measurement_equation": measurement_equation,
        "values": values,
        "standard_uncertainties": uncertainties,
        "units": units,
        "estimate": estimate,
        "sensitivities": analytic_sensitivities,
        "numeric_sensitivity_check": numeric_sensitivities,
        "max_sensitivity_error": max_sensitivity_error,
        "variance_terms": variance_terms,
        "local_variance": variance,
        "local_standard_uncertainty": math.sqrt(variance),
        "local_precision": 1.0 / variance,
        "positivity_checks": positivity_checks,
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["static storage Hessian", "hidden_load without a ceiling"],
        "boundary": boundary,
    }


def beer_lambert_case() -> dict[str, Any]:
    values = {"A_abs": 0.72, "epsilon": 15000.0, "path_length": 1.0}
    uncertainties = {"A_abs": 0.01, "epsilon": 300.0, "path_length": 0.005}

    def concentration(v: dict[str, float]) -> float:
        return v["A_abs"] / (v["epsilon"] * v["path_length"])

    estimate = concentration(values)
    sensitivities = {
        "A_abs": 1.0 / (values["epsilon"] * values["path_length"]),
        "epsilon": -estimate / values["epsilon"],
        "path_length": -estimate / values["path_length"],
    }
    return _linearised_case(
        case_id="beer_lambert_concentration_measurement_equation",
        measurement_equation="c = A_abs / (epsilon path_length)",
        values=values,
        uncertainties=uncertainties,
        units={
            "A_abs": "dimensionless absorbance",
            "epsilon": "L mol^-1 cm^-1",
            "path_length": "cm",
            "c": "mol L^-1",
        },
        fn=concentration,
        analytic_sensitivities=sensitivities,
        positivity_checks={
            "epsilon": "positive denominator; relative standard uncertainty is 2 percent",
            "path_length": "positive denominator; relative standard uncertainty is 0.5 percent",
        },
        boundary="This is a local concentration precision, not a storage Hessian; it must remain inside the calibrated optical range.",
    )


def cadmium_standard_case() -> dict[str, Any]:
    values = {"mass_cd_mg": 100.28, "purity": 0.9999, "volume_ml": 100.0}
    uncertainties = {"mass_cd_mg": 0.05, "purity": 0.000058, "volume_ml": 0.07}

    def concentration(v: dict[str, float]) -> float:
        return 1000.0 * v["mass_cd_mg"] * v["purity"] / v["volume_ml"]

    estimate = concentration(values)
    sensitivities = {
        "mass_cd_mg": estimate / values["mass_cd_mg"],
        "purity": estimate / values["purity"],
        "volume_ml": -estimate / values["volume_ml"],
    }
    return _linearised_case(
        case_id="cadmium_calibration_standard_measurement_equation",
        measurement_equation="c_Cd = 1000 mass_cd_mg purity / volume_ml",
        values=values,
        uncertainties=uncertainties,
        units={
            "mass_cd_mg": "mg",
            "purity": "dimensionless",
            "volume_ml": "mL",
            "c_Cd": "mg L^-1",
        },
        fn=concentration,
        analytic_sensitivities=sensitivities,
        positivity_checks={
            "volume_ml": "positive denominator; relative standard uncertainty is 0.07 percent",
            "purity": "bounded positive factor near one",
        },
        boundary="This is an uncertainty-propagation object for a calibration solution, not a physical hidden-load object.",
    )


def bromine_negative_control() -> dict[str, Any]:
    x = 0.5
    u = 0.005
    estimate_first_order = 2.0 * x * (1.0 - x)
    derivative = 2.0 - 4.0 * x
    first_order_variance = (derivative * u) ** 2
    second_order_mean = 0.5 - 2.0 * u**2
    second_order_standard_uncertainty = math.sqrt(8.0 * u**4)
    return {
        "id": "bromine_abundance_linearisation_negative_control",
        "route": "measurement_equation_local_linearisation",
        "status": "first_order_refused",
        "measurement_equation": "x160 = 2 x79 (1 - x79)",
        "values": {"x79": x},
        "standard_uncertainties": {"x79": u},
        "estimate_first_order": estimate_first_order,
        "sensitivity": {"x79": derivative},
        "first_order_variance": first_order_variance,
        "second_order_mean_for_normal_local_model": second_order_mean,
        "second_order_standard_uncertainty_for_normal_local_model": second_order_standard_uncertainty,
        "allowed_route": "higher_order_or_monte_carlo_only",
        "refused_routes": [
            "first_order local precision",
            "static storage Hessian",
            "hidden_load without a ceiling",
        ],
        "lesson": "A zero first derivative can make first-order uncertainty propagation report zero variance even when the nonlinear model has nonzero uncertainty.",
    }


def _markdown(proofs: list[dict[str, Any]]) -> str:
    lines = [
        "# Measurement Equation Linearisation Proofs v0",
        "",
        "Status: generated proof packet for measurement-equation local linearisation.",
        "",
        "These are not attachment cards. They test when a measurement equation can emit a local positive curvature object.",
        "",
    ]
    for proof in proofs:
        lines.extend(
            [
                f"## {proof['id']}",
                "",
                f"- route: `{proof['route']}`",
                f"- status: `{proof['status']}`",
                f"- equation: `{proof['measurement_equation']}`",
            ]
        )
        if "estimate" in proof:
            lines.extend(
                [
                    f"- estimate: `{proof['estimate']:.12g}`",
                    f"- local standard uncertainty: `{proof['local_standard_uncertainty']:.12g}`",
                    f"- local precision: `{proof['local_precision']:.12g}`",
                    f"- max sensitivity error: `{proof['max_sensitivity_error']:.12g}`",
                    f"- boundary: {proof['boundary']}",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    f"- first-order variance: `{proof['first_order_variance']:.12g}`",
                    f"- second-order standard uncertainty: `{proof['second_order_standard_uncertainty_for_normal_local_model']:.12g}`",
                    f"- allowed route: `{proof['allowed_route']}`",
                    f"- lesson: {proof['lesson']}",
                    "",
                ]
            )
    return "\n".join(lines)


def main() -> int:
    proofs = [
        beer_lambert_case(),
        cadmium_standard_case(),
        bromine_negative_control(),
    ]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "measurement_equation_linearisation_proofs.json"
    md_path = OUT_DIR / "measurement_equation_linearisation_proofs.md"
    payload = {
        "schema_version": "measurement_equation_linearisation_proofs.v0",
        "purpose": "Proof packet for measurement-equation local linearisation before any evidence schema is promoted.",
        "proofs": proofs,
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(_markdown(proofs), encoding="utf-8")
    print(json.dumps({"proofs": len(proofs), "outputs": [str(json_path), str(md_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
