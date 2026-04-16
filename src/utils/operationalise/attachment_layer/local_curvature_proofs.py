from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _positive_scalar(value: float, label: str) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise AssertionError(f"{label} must be finite and positive, got {value}")


def _rc_decay_proof() -> dict[str, Any]:
    tau = 2.0
    v0 = 5.0
    sigma = 0.05
    times = np.linspace(0.25, 8.0, 16)

    model = v0 * np.exp(-times / tau)
    derivative_tau = model * times / (tau * tau)
    fisher_tau = float(np.sum((derivative_tau / sigma) ** 2))
    _positive_scalar(fisher_tau, "RC decay Fisher curvature")

    return {
        "id": "rc_decay_local_curvature",
        "source_hook": "V(t; tau) = V0 exp(-t/tau)",
        "bare_formula_refusal": "The decay formula alone is a time-evolution surface, not a module object.",
        "declarations": {
            "coordinate": "tau",
            "coordinate_meaning": "RC time constant",
            "locality": "finite declared sample times",
            "curvature_source": "Gaussian voltage likelihood",
            "sensor_surface": "voltage samples V(t_i)",
            "noise_model": "independent Gaussian voltage noise with standard deviation sigma",
            "known_controls": {
                "V0": v0,
                "sample_times": times.tolist(),
                "sigma": sigma
            }
        },
        "derivation": [
            "For independent Gaussian observations with fixed sigma, the local Fisher curvature is sum_i (1/sigma^2) (dV(t_i;tau)/dtau)^2.",
            "For V(t;tau) = V0 exp(-t/tau), dV/dtau = V(t;tau) * t / tau^2.",
            "Substituting the declared sample grid gives a positive scalar Fisher object."
        ],
        "positive_object": {
            "kind": "scalar Fisher precision",
            "coordinate": "tau",
            "value": fisher_tau,
            "units_note": "inverse tau^2 units under the declared voltage-noise scale"
        },
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["hidden_load without a ceiling", "static storage Hessian"],
        "boundary": "This proof does not infer tau from data and does not model resistor/capacitor storage; it only shows admissibility of declared local information curvature.",
    }


def _beer_lambert_proof() -> dict[str, Any]:
    epsilon = 1.2
    path_length = 1.0
    sigma = 0.01
    derivative_c = epsilon * path_length
    fisher_c = float((derivative_c / sigma) ** 2)
    _positive_scalar(fisher_c, "Beer-Lambert Fisher curvature")

    return {
        "id": "beer_lambert_local_curvature",
        "source_hook": "A(c) = epsilon b c",
        "bare_formula_refusal": "The absorbance formula alone is a sensor law, not a module object.",
        "declarations": {
            "coordinate": "c",
            "coordinate_meaning": "solute concentration",
            "locality": "single declared absorbance measurement protocol; linear law is global in this coordinate under the stated regime",
            "curvature_source": "Gaussian absorbance likelihood",
            "sensor_surface": "absorbance A",
            "noise_model": "Gaussian absorbance noise with standard deviation sigma_A",
            "known_controls": {
                "molar_absorptivity": epsilon,
                "path_length": path_length,
                "sigma_absorbance": sigma
            }
        },
        "derivation": [
            "For one Gaussian absorbance observation with fixed sigma_A, the Fisher curvature is (1/sigma_A^2) (dA/dc)^2.",
            "For A(c) = epsilon b c, dA/dc = epsilon b.",
            "Substituting the declared sensor constants gives a positive scalar Fisher object."
        ],
        "positive_object": {
            "kind": "scalar Fisher precision",
            "coordinate": "c",
            "value": fisher_c,
            "units_note": "inverse concentration^2 units under the declared absorbance-noise scale"
        },
        "allowed_route": "local_positive_curvature",
        "refused_routes": ["hidden_load without a ceiling", "static storage Hessian"],
        "boundary": "This proof does not model chemistry thermodynamics or concentration dynamics; it only shows admissibility of declared local information curvature from a sensor law.",
    }


def _render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Local Positive Curvature Proofs v0",
        "",
        "These examples test the admissibility idea before promoting any evidence schema.",
        "",
        "A bare formula is refused. A declared local positive curvature object is admitted.",
        "",
    ]
    for proof in payload["proofs"]:
        lines.extend(
            [
                f"## `{proof['id']}`",
                "",
                f"- source hook: `{proof['source_hook']}`",
                f"- bare formula refusal: {proof['bare_formula_refusal']}",
                f"- allowed route: `{proof['allowed_route']}`",
                f"- refused routes: {', '.join(proof['refused_routes'])}",
                f"- positive object: `{proof['positive_object']['kind']}` on `{proof['positive_object']['coordinate']}` = {proof['positive_object']['value']:.12g}",
                f"- boundary: {proof['boundary']}",
                "",
                "### Declarations",
                "",
                "```json",
                json.dumps(proof["declarations"], indent=2),
                "```",
                "",
                "### Derivation",
                "",
            ]
        )
        for line in proof["derivation"]:
            lines.append(f"- {line}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate local positive curvature admissibility proofs.")
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "outputs" / "local_curvature"))
    args = parser.parse_args()

    payload = {
        "schema_version": "local_positive_curvature_proofs.v0",
        "description": "Proof-of-admissibility examples for declared local positive information curvature.",
        "proofs": [_rc_decay_proof(), _beer_lambert_proof()],
    }
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "local_positive_curvature_proofs.json"
    md_path = output_dir / "local_positive_curvature_proofs.md"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write(_render_markdown(payload))

    print(json.dumps({"proofs": len(payload["proofs"]), "outputs": [str(json_path), str(md_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
