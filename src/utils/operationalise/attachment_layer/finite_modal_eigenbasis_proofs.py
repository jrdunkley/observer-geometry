from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "finite_modal"


def _spd(matrix: np.ndarray, label: str) -> dict[str, Any]:
    eigvals = np.linalg.eigvalsh(matrix)
    if np.any(eigvals <= 0.0):
        raise AssertionError(f"{label} must be SPD, eigenvalues={eigvals}")
    return {
        "eigenvalues": [float(item) for item in eigvals],
        "condition_number": float(eigvals[-1] / eigvals[0]),
    }


def string_fixed_end_three_mode() -> dict[str, Any]:
    length = 1.0
    tension = 100.0
    linear_density = 0.01
    mode_count = 3
    modes = np.arange(1, mode_count + 1, dtype=float)

    # Basis: phi_n(x) = sqrt(2/L) sin(n pi x / L), so int_0^L phi_i phi_j dx = delta_ij.
    mass_metric = linear_density * np.eye(mode_count)
    stiffness = np.diag(tension * (modes * math.pi / length) ** 2)
    angular_frequencies = np.sqrt(np.diag(stiffness) / np.diag(mass_metric))
    frequencies = angular_frequencies / (2.0 * math.pi)

    return {
        "id": "string_fixed_end_three_mode",
        "route": "finite_modal_eigenbasis",
        "status": "admissible_finite_modal_family",
        "source_hook_refusal": "The raw standing-wave formula is infinite and boundary-dependent; it is not a finite module object.",
        "boundary_conditions": "fixed displacement at x = 0 and x = L",
        "basis": "orthonormal sine modes phi_n(x) = sqrt(2/L) sin(n pi x / L)",
        "truncation": {
            "mode_count": mode_count,
            "included_modes": [int(item) for item in modes],
            "boundary": "Modes above n = 3 are explicitly outside this finite proof object.",
        },
        "controls": {
            "length_m": length,
            "tension_N": tension,
            "linear_density_kg_per_m": linear_density,
        },
        "module_objects": {
            "mass_metric": mass_metric.tolist(),
            "stiffness_hessian": stiffness.tolist(),
            "mass_metric_spd": _spd(mass_metric, "string modal mass metric"),
            "stiffness_spd": _spd(stiffness, "string modal stiffness Hessian"),
        },
        "modal_readout": {
            "angular_frequencies_rad_per_s": [float(item) for item in angular_frequencies],
            "frequencies_hz": [float(item) for item in frequencies],
            "frequency_rule": "f_n = n/(2L) sqrt(T/mu)",
        },
        "sensor_surface": {
            "ideal": "modal amplitudes q_1, q_2, q_3",
            "practical": "point displacement, pickup, or microphone channel requires an additional observer map",
        },
        "allowed_route": "finite_modal_positive_curvature",
        "refused_routes": ["raw infinite wave formula", "hidden_load without a ceiling"],
        "boundary": "This proof admits only the declared three-mode truncation. A physical sensor still needs an observer map before visible precision or hidden-load claims.",
    }


def two_mode_projection_example() -> dict[str, Any]:
    length = 1.0
    sensor_x = 0.25
    modes = np.array([1.0, 2.0])
    observer = math.sqrt(2.0 / length) * np.sin(modes * math.pi * sensor_x / length)
    rank = int(np.linalg.matrix_rank(observer.reshape(1, -1)))
    return {
        "id": "string_point_sensor_two_mode_projection",
        "route": "finite_modal_eigenbasis",
        "status": "observer_map_recorded_not_hidden_load",
        "basis": "first two fixed-end sine modes",
        "sensor_surface": "point displacement y = phi_1(x0) q_1 + phi_2(x0) q_2 at x0 = L/4",
        "observer_map": [float(item) for item in observer],
        "rank": rank,
        "lesson": "A point sensor exposes a projection of modal coordinates. It is an observer map, but it is not enough for hidden-load without a declared ceiling/reference.",
        "allowed_route": "observer_map_record_only",
        "refused_routes": ["hidden_load without a ceiling"],
    }


def acoustic_tube_three_mode_pressure_observer() -> dict[str, Any]:
    length = 1.0
    area = 0.01
    density = 1.2
    sound_speed = 343.0
    sensor_x = length / 3.0
    mode_count = 3
    modes = np.arange(1, mode_count + 1, dtype=float)

    mass_metric = density * area * np.eye(mode_count)
    stiffness = np.diag(density * sound_speed * sound_speed * area * (modes * math.pi / length) ** 2)
    observer = -density * sound_speed * sound_speed * np.sqrt(2.0 / length) * (modes * math.pi / length) * np.cos(modes * math.pi * sensor_x / length)
    visible_precision = np.linalg.inv(observer.reshape(1, -1) @ np.linalg.inv(stiffness) @ observer.reshape(-1, 1))
    angular_frequencies = np.sqrt(np.diag(stiffness) / np.diag(mass_metric))
    frequencies = angular_frequencies / (2.0 * math.pi)

    return {
        "id": "acoustic_tube_three_mode_pressure_observer",
        "route": "finite_modal_eigenbasis",
        "status": "admissible_finite_modal_pressure_observer",
        "source_hook_refusal": "Raw acoustic impedance, intensity, attenuation, and ultrasound reflection formulas are not admitted as module objects without a finite acoustic field model and sensor convention.",
        "boundary_conditions": "rigid closed ends with zero longitudinal displacement at x = 0 and x = L",
        "basis": "orthonormal longitudinal displacement sine modes phi_n(x) = sqrt(2/L) sin(n pi x / L)",
        "truncation": {
            "mode_count": mode_count,
            "included_modes": [int(item) for item in modes],
            "boundary": "Modes above n = 3, damping, radiation losses, and drive response are explicitly outside this finite proof object.",
        },
        "controls": {
            "length_m": length,
            "area_m2": area,
            "density_kg_per_m3": density,
            "sound_speed_m_per_s": sound_speed,
            "sensor_x_m": sensor_x,
        },
        "module_objects": {
            "mass_metric": mass_metric.tolist(),
            "stiffness_hessian": stiffness.tolist(),
            "pressure_observer": [float(item) for item in observer],
            "visible_pressure_precision": visible_precision.tolist(),
            "mass_metric_spd": _spd(mass_metric, "acoustic modal mass metric"),
            "stiffness_spd": _spd(stiffness, "acoustic modal stiffness Hessian"),
        },
        "modal_readout": {
            "angular_frequencies_rad_per_s": [float(item) for item in angular_frequencies],
            "frequencies_hz": [float(item) for item in frequencies],
            "frequency_rule": "f_n = n c / (2L) for the declared rigid-ended tube displacement modes",
        },
        "sensor_surface": {
            "ideal": "point pressure p(x0) = -rho c^2 d xi / dx at x0 = L/3",
            "practical": "calibrated small microphone approximated as a point pressure sensor",
        },
        "allowed_route": "finite_modal_pressure_observer",
        "refused_routes": ["raw acoustic impedance", "ultrasound reflection", "hidden_load without a ceiling"],
        "boundary": "This proof admits only the declared three-mode closed-tube pressure observer. It does not admit impedance, reflection, attenuation, damping, or full field reconstruction.",
    }


def negative_controls() -> list[dict[str, Any]]:
    zero_tension_stiffness = np.diag([0.0, 0.0, 0.0])
    no_truncation_control = {
        "id": "raw_standing_wave_no_truncation_refusal",
        "route": "finite_modal_eigenbasis",
        "status": "refused",
        "reason": "A standing-wave formula without boundary conditions, basis, and finite truncation is not a finite matrix object.",
        "refused_route": "finite_modal_positive_curvature",
    }
    zero_tension_control = {
        "id": "string_zero_tension_stiffness_refusal",
        "route": "finite_modal_eigenbasis",
        "status": "refused",
        "reason": "Zero tension gives zero stiffness eigenvalues, so the declared stiffness Hessian is not positive.",
        "stiffness_eigenvalues": [float(item) for item in np.linalg.eigvalsh(zero_tension_stiffness)],
        "refused_route": "finite_modal_positive_curvature",
    }
    return [no_truncation_control, zero_tension_control]


def _markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Finite Modal Eigenbasis Proofs v0",
        "",
        "Status: generated proof packet for finite modal attachment.",
        "",
        "These are not attachment cards. They test when a wave surface becomes a finite positive matrix object.",
        "",
    ]
    for proof in payload["proofs"]:
        lines.extend(
            [
                f"## {proof['id']}",
                "",
                f"- route: `{proof['route']}`",
                f"- status: `{proof['status']}`",
            ]
        )
        if "module_objects" in proof:
            lines.extend(
                [
                    f"- stiffness eigenvalues: `{proof['module_objects']['stiffness_spd']['eigenvalues']}`",
                    f"- frequencies Hz: `{proof['modal_readout']['frequencies_hz']}`",
                    f"- boundary: {proof['boundary']}",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    f"- observer map: `{proof['observer_map']}`",
                    f"- lesson: {proof['lesson']}",
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
        "schema_version": "finite_modal_eigenbasis_proofs.v0",
        "purpose": "Proof packet for finite modal positive curvature across the string and acoustic wave cards.",
        "proofs": [string_fixed_end_three_mode(), two_mode_projection_example(), acoustic_tube_three_mode_pressure_observer()],
        "negative_controls": negative_controls(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "finite_modal_eigenbasis_proofs.json"
    md_path = OUT_DIR / "finite_modal_eigenbasis_proofs.md"
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
