from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def _positive(value: float, name: str) -> None:
    if not value > 0.0:
        raise AssertionError(f"{name} should be positive, got {value}")


def _rc_decay_probe() -> dict[str, Any]:
    tau = 2.0
    v0 = 5.0
    sigma = 0.05
    times = np.linspace(0.25, 8.0, 16)
    values = v0 * np.exp(-times / tau)
    derivative = values * times / (tau * tau)
    fisher_tau = float(np.sum((derivative / sigma) ** 2))
    _positive(fisher_tau, "rc fisher_tau")
    return {
        "id": "rc_charge_decay",
        "route_class": "dynamic_evidence",
        "surface_formula": "V(t) = V0 exp(-t/tau)",
        "refusal_without_extra_declarations": "The exponential law is not a static Hessian.",
        "minimal_extra_declarations": ["sample times", "Gaussian voltage noise sigma", "fit parameter tau"],
        "declared_object_after_extra_declarations": "scalar Fisher precision for tau",
        "numeric_readout": {
            "tau": tau,
            "v0": v0,
            "sigma": sigma,
            "sample_count": int(times.size),
            "fisher_tau": fisher_tau,
        },
        "lesson": "The failure points to an evidence/Fisher attachment schema, not to a static hidden-load card.",
    }


def _beer_lambert_probe() -> dict[str, Any]:
    epsilon = 1.2
    path_length = 1.0
    sigma_a = 0.01
    fisher_c = float((epsilon * path_length / sigma_a) ** 2)
    _positive(fisher_c, "beer-lambert fisher_c")
    return {
        "id": "absorbance_beer_lambert",
        "route_class": "sensor_or_observer_map",
        "surface_formula": "A = epsilon * b * c",
        "refusal_without_extra_declarations": "The absorbance equation is a sensor law, not a storage Hessian.",
        "minimal_extra_declarations": ["molar absorptivity", "path length", "Gaussian absorbance noise", "concentration parameter"],
        "declared_object_after_extra_declarations": "scalar measurement Fisher precision for concentration",
        "numeric_readout": {
            "molar_absorptivity": epsilon,
            "path_length": path_length,
            "sigma_absorbance": sigma_a,
            "fisher_concentration": fisher_c,
        },
        "lesson": "Sensor laws can become module-facing evidence only after a noise model and estimated parameter are declared.",
    }


def _nuclear_decay_probe() -> dict[str, Any]:
    decay_rate = 0.18
    initial_rate = 200.0
    times = np.arange(1.0, 9.0)
    expected_counts = initial_rate * np.exp(-decay_rate * times)
    derivative = -times * expected_counts
    fisher_lambda = float(np.sum((derivative**2) / expected_counts))
    _positive(fisher_lambda, "nuclear fisher_lambda")
    return {
        "id": "nuclear_decay_poisson",
        "route_class": "dynamic_evidence",
        "surface_formula": "N(t) = N0 exp(-lambda t)",
        "refusal_without_extra_declarations": "The decay law is a stochastic counting process, not a static quadratic storage object.",
        "minimal_extra_declarations": ["counting windows", "Poisson likelihood", "initial rate", "decay parameter lambda"],
        "declared_object_after_extra_declarations": "scalar Poisson Fisher precision for lambda",
        "numeric_readout": {
            "lambda": decay_rate,
            "initial_rate": initial_rate,
            "sample_count": int(times.size),
            "fisher_lambda": fisher_lambda,
        },
        "lesson": "This is a strong evidence hook; the missing object is likelihood structure, not more formula digestion.",
    }


def _standing_wave_probe() -> dict[str, Any]:
    tension = 20.0
    length = 1.0
    linear_density = 0.05
    modes = np.array([1.0, 2.0, 3.0])
    stiffness = tension * (modes * math.pi / length) ** 2
    mass_metric = np.full_like(stiffness, linear_density * length / 2.0)
    omega_squared = stiffness / mass_metric
    eig_min = float(np.min(stiffness))
    _positive(eig_min, "standing-wave stiffness min")
    return {
        "id": "standing_wave_string_modes",
        "route_class": "finite_weighted_family",
        "surface_formula": "string harmonic frequencies and wave speed",
        "refusal_without_extra_declarations": "The wave formula is not finite until boundary conditions, basis, and truncation are declared.",
        "minimal_extra_declarations": ["fixed-end boundary", "sine mode basis", "three-mode truncation", "tension", "linear density"],
        "declared_object_after_extra_declarations": "finite diagonal modal stiffness family with mass metric",
        "numeric_readout": {
            "modes": modes.astype(int).tolist(),
            "stiffness_diagonal": stiffness.tolist(),
            "mass_metric_diagonal": mass_metric.tolist(),
            "omega_squared": omega_squared.tolist(),
            "min_stiffness": eig_min,
        },
        "lesson": "The obstruction is infinite/basis-free description; a finite mode declaration converts it into a clean quadratic family.",
    }


def _gravity_local_probe() -> dict[str, Any]:
    mu = 1.0
    r0 = 2.0
    raw_second_derivative = -2.0 * mu / (r0**3)
    angular_momentum_squared = mu * r0
    effective_second_derivative = 3.0 * angular_momentum_squared / (r0**4) - 2.0 * mu / (r0**3)
    if raw_second_derivative >= 0.0:
        raise AssertionError("raw gravitational second derivative should be negative in this probe")
    _positive(effective_second_derivative, "effective gravitational Hessian")
    return {
        "id": "gravity_inverse_square_local",
        "route_class": "local_quadratic",
        "surface_formula": "U(r) = -mu/r",
        "refusal_without_extra_declarations": "The raw inverse-square potential is not globally positive quadratic; local radial curvature has the wrong sign in this simple coordinate.",
        "minimal_extra_declarations": ["operating radius", "effective potential", "angular momentum/constraint convention", "local stability sector"],
        "declared_object_after_extra_declarations": "positive local radial Hessian of an effective potential at a stable circular orbit",
        "numeric_readout": {
            "mu": mu,
            "r0": r0,
            "raw_potential_second_derivative": raw_second_derivative,
            "effective_potential_second_derivative": effective_second_derivative,
        },
        "lesson": "The failure is informative: local field attachment depends on stability and effective-potential declarations.",
    }


def _ideal_gas_probe() -> dict[str, Any]:
    n_moles = 1.0
    gas_constant = 8.314462618
    temperature = 300.0
    volume = 0.024
    helmholtz_volume_hessian = n_moles * gas_constant * temperature / (volume**2)
    _positive(helmholtz_volume_hessian, "ideal gas Helmholtz Hessian")
    return {
        "id": "ideal_gas_law_boundary",
        "route_class": "thermodynamic_fluctuation",
        "surface_formula": "p V = n R T",
        "refusal_without_extra_declarations": "The equation of state alone is not a module object.",
        "minimal_extra_declarations": ["fixed n,T ensemble", "Helmholtz free-energy model", "volume coordinate", "operating point"],
        "declared_object_after_extra_declarations": "local scalar Hessian d2F/dV2 in declared isothermal Helmholtz model",
        "numeric_readout": {
            "n_moles": n_moles,
            "temperature_K": temperature,
            "volume_m3": volume,
            "helmholtz_volume_hessian": helmholtz_volume_hessian,
        },
        "lesson": "Thermodynamic formulas attach only after ensemble and potential are declared; the refusal remains correct.",
    }


def _electrochemistry_probe() -> dict[str, Any]:
    gas_constant = 8.314462618
    temperature = 298.15
    faraday = 96485.33212
    electrons = 1.0
    sigma_v = 0.002
    derivative_wrt_log_q = -gas_constant * temperature / (electrons * faraday)
    fisher_log_q = float((derivative_wrt_log_q / sigma_v) ** 2)
    _positive(fisher_log_q, "electrochemistry fisher_log_q")
    return {
        "id": "electrochemistry_nernst_voltage",
        "route_class": "thermodynamic_fluctuation",
        "surface_formula": "E = E0 - (RT/nF) ln Q",
        "refusal_without_extra_declarations": "Voltage is measurable, but reaction coordinate and thermodynamic convention are missing.",
        "minimal_extra_declarations": ["reaction quotient coordinate log Q", "electron count", "temperature", "voltage noise"],
        "declared_object_after_extra_declarations": "scalar voltage-measurement Fisher precision for log Q",
        "numeric_readout": {
            "temperature_K": temperature,
            "electrons": electrons,
            "sigma_voltage": sigma_v,
            "dE_dlogQ": derivative_wrt_log_q,
            "fisher_logQ": fisher_log_q,
        },
        "lesson": "Electrochemistry is a strong bridge, but it should enter through a declared thermodynamic/evidence schema.",
    }


def _relativity_metric_probe() -> dict[str, Any]:
    metric = np.diag([-1.0, 1.0, 1.0, 1.0])
    eigenvalues = np.linalg.eigvalsh(metric)
    if float(np.min(eigenvalues)) >= 0.0:
        raise AssertionError("Lorentz metric probe should be indefinite")
    return {
        "id": "relativity_lorentz_metric",
        "route_class": "exact_special_sector",
        "surface_formula": "Minkowski metric diag(-1,1,1,1)",
        "refusal_without_extra_declarations": "The current static kernel calls require SPD objects; this metric is indefinite.",
        "minimal_extra_declarations": ["positive subproblem or special sector", "signature convention", "observable map"],
        "declared_object_after_extra_declarations": "none in v0 static attachment runner",
        "numeric_readout": {
            "metric_diagonal": np.diag(metric).tolist(),
            "eigenvalues": eigenvalues.tolist(),
            "min_eigenvalue": float(np.min(eigenvalues)),
        },
        "lesson": "This is a real boundary for the current SPD route, not merely missing units.",
    }


PROBES = [
    _rc_decay_probe,
    _beer_lambert_probe,
    _nuclear_decay_probe,
    _standing_wave_probe,
    _gravity_local_probe,
    _ideal_gas_probe,
    _electrochemistry_probe,
    _relativity_metric_probe,
]


def _render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Failure Probe Suite v0",
        "",
        "These are not attachment cards. They interrogate representative failures to determine what extra declaration would be needed before a module-facing object exists.",
        "",
        "| probe | route | declared object after extra declarations | lesson |",
        "| --- | --- | --- | --- |",
    ]
    for probe in payload["probes"]:
        lines.append(
            "| `{id}` | `{route_class}` | {declared_object_after_extra_declarations} | {lesson} |".format(**probe)
        )
    lines.extend(["", "## Readouts", ""])
    for probe in payload["probes"]:
        lines.extend(
            [
                f"### `{probe['id']}`",
                "",
                f"- surface formula: `{probe['surface_formula']}`",
                f"- refusal: {probe['refusal_without_extra_declarations']}",
                f"- minimal extra declarations: {', '.join(probe['minimal_extra_declarations'])}",
                "",
                "```json",
                json.dumps(probe["numeric_readout"], indent=2),
                "```",
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Interrogate representative non-integrating physical hooks.")
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "outputs" / "failure_probes"))
    args = parser.parse_args()

    probes = [probe() for probe in PROBES]
    payload = {
        "schema_version": "failure_probe_suite.v0",
        "description": "Representative probes for non-integrating AP/A-level/formulary hooks.",
        "probes": probes,
    }

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "failure_probe_suite.json"
    md_path = output_dir / "failure_probe_suite.md"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write(_render_markdown(payload))

    print(json.dumps({"probes": len(probes), "outputs": [str(json_path), str(md_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
