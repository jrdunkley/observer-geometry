from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _scalar_matrix_value(payload: dict[str, Any], field: str) -> float:
    return float(payload["readout"][field][0][0])


def build_demonstration(root: Path) -> dict[str, Any]:
    spring = _load_json(root / "outputs" / "mass_spring_coupled_observe_one.json")
    lc = _load_json(root / "outputs" / "coupled_lc_resonators_observe_one_node.json")
    string_fixed = _load_json(root / "outputs" / "string_fixed_end_three_mode.json")
    string_point = _load_json(root / "outputs" / "string_point_sensor_modal_observer.json")
    string_hidden = _load_json(root / "outputs" / "string_modal_hidden_load_with_declared_ceiling.json")
    rc = _load_json(root / "outputs" / "empirical_datasets" / "rc_decay_q7_probe.json")
    beer = _load_json(root / "outputs" / "empirical_datasets" / "beer_lambert_uvvis_probe.json")

    rc_fit = rc["fixed_V0_voltage_fit"]
    beer_cal = beer["calibration"]
    beer_summary = beer["unknown_sample_summary"]

    return {
        "schema_version": "gold_quartet_demonstration.v0",
        "status": "external_facing_demonstration",
        "purpose": "Compact evidence bundle showing that declared local positive contact is a route discipline rather than a list of examples.",
        "demonstrations": [
            {
                "id": "cross_substrate_hidden_renormalisation",
                "route": "reference_ceiling_hidden_load",
                "source_outputs": [
                    "outputs/mass_spring_coupled_observe_one.json",
                    "outputs/coupled_lc_resonators_observe_one_node.json",
                ],
                "key_readout": {
                    "spring_visible_precision": _scalar_matrix_value(spring, "visible_precision"),
                    "spring_ceiling_minus_visible_precision": _scalar_matrix_value(spring, "ceiling_minus_visible_precision"),
                    "spring_hidden_load": _scalar_matrix_value(spring, "hidden_load"),
                    "lc_visible_precision": _scalar_matrix_value(lc, "visible_precision"),
                    "lc_ceiling_minus_visible_precision": _scalar_matrix_value(lc, "ceiling_minus_visible_precision"),
                    "lc_hidden_load": _scalar_matrix_value(lc, "hidden_load"),
                },
                "lesson": "Different substrates produce the same observer-facing pattern only after partial observation and a declared visible reference ceiling are present.",
                "mandatory_refusal": "hidden-load language without a ceiling/reference",
            },
            {
                "id": "finite_modal_sensor_reference",
                "route": "finite_modal_curvature_to_reference_ceiling_hidden_load",
                "source_outputs": [
                    "outputs/string_fixed_end_three_mode.json",
                    "outputs/string_point_sensor_modal_observer.json",
                    "outputs/string_modal_hidden_load_with_declared_ceiling.json",
                ],
                "key_readout": {
                    "three_mode_stiffness_diagonal": [
                        float(string_fixed["readout"]["visible_precision"][index][index])
                        for index in range(3)
                    ],
                    "point_sensor_visible_precision": _scalar_matrix_value(string_point, "visible_precision"),
                    "reference_ceiling_minus_visible_precision": _scalar_matrix_value(string_hidden, "ceiling_minus_visible_precision"),
                    "reference_hidden_load": _scalar_matrix_value(string_hidden, "hidden_load"),
                },
                "lesson": "Continuum wave language becomes module-facing only after boundary, basis, truncation, observer row, and reference ceiling are declared.",
                "mandatory_refusal": "raw wave formula, full modal inference from one point sensor, or hidden load without reference",
            },
            {
                "id": "rc_decay_information_curvature",
                "route": "observation_equation_fisher",
                "source_outputs": ["outputs/empirical_datasets/rc_decay_q7_probe.json"],
                "key_readout": {
                    "tau_hat_s": float(rc_fit["tau_hat_s"]),
                    "fisher_tau": float(rc_fit["fisher_tau"]),
                    "local_variance_tau_s2": float(rc_fit["local_variance_tau_s2"]),
                    "data_warnings": len(rc["data_warnings"]),
                },
                "lesson": "A voltage decay trace can become local Fisher precision for tau after declaring V0, positive tau, sample times, and a residual-estimated noise convention.",
                "mandatory_refusal": "storage Hessian, hidden load, capacitance without separately declared resistance, and instrument-noise claim without instrument metadata",
            },
            {
                "id": "beer_lambert_sensor_information_curvature",
                "route": "observation_equation_fisher_or_analysis_function",
                "source_outputs": ["outputs/empirical_datasets/beer_lambert_uvvis_probe.json"],
                "key_readout": {
                    "lambda_star_nm": float(beer["lambda_star_nm"]),
                    "known_observation_count": int(beer_cal["known_observation_count"]),
                    "beta_absorbance_per_mM": float(beer_cal["beta_absorbance_per_mM"]),
                    "sigma_absorbance_residual": float(beer_cal["sigma_absorbance_residual"]),
                    "mean_estimated_concentration_mM": float(beer_summary["mean_estimated_concentration_mM"]),
                    "inside_known_calibration_span": bool(beer_summary["inside_known_calibration_span"]),
                },
                "lesson": "Spectra and standards become a local concentration/calibration object after wavelength, standards, residual noise, and calibration span are declared.",
                "mandatory_refusal": "storage Hessian, hidden load, thermodynamic concentration model, and molar absorptivity without path length",
            },
        ],
        "scientific_conclusion": "The same local-positive interface is reached through energy storage, finite modal declaration, and information curvature, but each route licenses a different kind of claim and carries a neighbouring refusal.",
    }


def render_markdown(demo: dict[str, Any]) -> str:
    by_id = {item["id"]: item for item in demo["demonstrations"]}
    cross = by_id["cross_substrate_hidden_renormalisation"]["key_readout"]
    modal = by_id["finite_modal_sensor_reference"]["key_readout"]
    rc = by_id["rc_decay_information_curvature"]["key_readout"]
    beer = by_id["beer_lambert_sensor_information_curvature"]["key_readout"]

    return f"""# Gold Quartet Demonstration v0

Status: compact external-facing demonstration.

This file records the smallest current evidence bundle for the attachment layer. It is not an expanded theory digest. It shows the same contact discipline across four routes:

```text
declared physical surface
-> finite positive object
-> observer or sensor route
-> licensed readout
-> neighbouring refusal
```

The source of truth for the numeric readouts is `gold_quartet_demonstration_v0.json`, generated by `gold_quartet_demonstration.py`, and checked against the generated card and empirical-probe outputs by `validate_operationalise_metadata.py`.

## 1. Cross-Substrate Hidden Renormalisation

Source outputs:

- `outputs/mass_spring_coupled_observe_one.json`
- `outputs/coupled_lc_resonators_observe_one_node.json`

Readout:

| Substrate | Visible Precision | Ceiling Minus Visible | Hidden Load |
| --- | ---: | ---: | ---: |
| Coupled springs | {cross["spring_visible_precision"]} | {cross["spring_ceiling_minus_visible_precision"]} | {cross["spring_hidden_load"]} |
| Coupled LC resonators | {cross["lc_visible_precision"]} | {cross["lc_ceiling_minus_visible_precision"]} | {cross["lc_hidden_load"]} |

Lesson: different substrates produce the same observer-facing pattern only after partial observation and a declared visible reference ceiling are present.

Refusal: hidden-load language without a ceiling/reference.

## 2. Finite Modal Sensor Reference

Source outputs:

- `outputs/string_fixed_end_three_mode.json`
- `outputs/string_point_sensor_modal_observer.json`
- `outputs/string_modal_hidden_load_with_declared_ceiling.json`

Readout:

| Quantity | Value |
| --- | ---: |
| Three-mode stiffness diagonal 1 | {modal["three_mode_stiffness_diagonal"][0]} |
| Three-mode stiffness diagonal 2 | {modal["three_mode_stiffness_diagonal"][1]} |
| Three-mode stiffness diagonal 3 | {modal["three_mode_stiffness_diagonal"][2]} |
| Point-sensor visible precision | {modal["point_sensor_visible_precision"]} |
| Reference ceiling minus visible precision | {modal["reference_ceiling_minus_visible_precision"]} |
| Reference hidden load | {modal["reference_hidden_load"]} |

Lesson: continuum wave language becomes module-facing only after boundary, basis, truncation, observer row, and reference ceiling are declared.

Refusal: raw wave formula, full modal inference from one point sensor, or hidden load without reference.

## 3. RC Decay Information Curvature

Source output:

- `outputs/empirical_datasets/rc_decay_q7_probe.json`

Readout:

| Quantity | Value |
| --- | ---: |
| `tau_hat_s` | {rc["tau_hat_s"]} |
| `fisher_tau` | {rc["fisher_tau"]} |
| `local_variance_tau_s2` | {rc["local_variance_tau_s2"]} |
| data warnings | {rc["data_warnings"]} |

Lesson: a voltage decay trace can become local Fisher precision for `tau` after declaring `V0`, positive `tau`, sample times, and a residual-estimated noise convention.

Refusal: storage Hessian, hidden load, capacitance without separately declared resistance, and instrument-noise claim without instrument metadata.

## 4. Beer-Lambert Sensor Information Curvature

Source output:

- `outputs/empirical_datasets/beer_lambert_uvvis_probe.json`

Readout:

| Quantity | Value |
| --- | ---: |
| `lambda_star_nm` | {beer["lambda_star_nm"]} |
| known observation count | {beer["known_observation_count"]} |
| `beta_absorbance_per_mM` | {beer["beta_absorbance_per_mM"]} |
| `sigma_absorbance_residual` | {beer["sigma_absorbance_residual"]} |
| mean estimated concentration mM | {beer["mean_estimated_concentration_mM"]} |
| inside known calibration span | {str(beer["inside_known_calibration_span"]).lower()} |

Lesson: spectra and standards become a local concentration/calibration object after wavelength, standards, residual noise, and calibration span are declared.

Refusal: storage Hessian, hidden load, thermodynamic concentration model, and molar absorptivity without path length.

## Conclusion

The important result is not that four cases work. It is that three different provenance routes now converge on finite positive contact while preserving different licensed claims:

- storage/reference curvature licenses visible precision and hidden-load readout only with a ceiling;
- finite modal curvature licenses finite modal and observer-map readouts only after boundary, basis, and truncation;
- empirical information curvature licenses local parameter or analysis-function precision, not storage geometry.

The refusals are part of the result. They show that the interface is beginning to act like a contact discipline rather than a permissive encoding language.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the gold quartet demonstration from validated output artifacts.")
    parser.add_argument("--root", default=str(ROOT))
    args = parser.parse_args()
    root = Path(args.root)

    demo = build_demonstration(root)
    (root / "gold_quartet_demonstration_v0.json").write_text(json.dumps(demo, indent=2) + "\n", encoding="utf-8")
    (root / "gold_quartet_demonstration_v0.md").write_text(render_markdown(demo), encoding="utf-8")
    print(json.dumps({"generated": "gold_quartet_demonstration_v0", "demonstrations": len(demo["demonstrations"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
