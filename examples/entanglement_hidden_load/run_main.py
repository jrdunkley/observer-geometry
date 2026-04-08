from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import scipy.linalg as la

from nomogeo import (
    clock,
    hidden_contraction,
    hidden_load,
    load_from_hidden_contraction,
    visible_precision,
)

from examples.common import line_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def selector_a() -> np.ndarray:
    return np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])


def two_mode_squeezer(r: float) -> np.ndarray:
    c = np.cosh(r)
    s = np.sinh(r)
    return np.array(
        [
            [c, 0.0, s, 0.0],
            [0.0, c, 0.0, -s],
            [s, 0.0, c, 0.0],
            [0.0, -s, 0.0, c],
        ]
    )


def local_squeezer(sa: float, sb: float) -> np.ndarray:
    return np.diag([np.exp(sa), np.exp(-sa), np.exp(sb), np.exp(-sb)])


def local_rotation(theta_a: float, theta_b: float) -> np.ndarray:
    rot_a = np.array([[np.cos(theta_a), -np.sin(theta_a)], [np.sin(theta_a), np.cos(theta_a)]])
    rot_b = np.array([[np.cos(theta_b), -np.sin(theta_b)], [np.sin(theta_b), np.cos(theta_b)]])
    return la.block_diag(rot_a, rot_b)


def two_mode_squeezed_thermal_covariance(r: float, n_a: float = 0.0, n_b: float = 0.0) -> np.ndarray:
    v_a = 2.0 * n_a + 1.0
    v_b = 2.0 * n_b + 1.0
    sigma0 = np.diag([v_a, v_a, v_b, v_b])
    s = two_mode_squeezer(r)
    sigma = s @ sigma0 @ s.T
    return 0.5 * (sigma + sigma.T)


def gaussian_mutual_information(sigma: np.ndarray) -> float:
    a = sigma[:2, :2]
    b = sigma[2:, 2:]
    sign_a, logdet_a = np.linalg.slogdet(a)
    sign_b, logdet_b = np.linalg.slogdet(b)
    sign_s, logdet_s = np.linalg.slogdet(sigma)
    if min(sign_a, sign_b, sign_s) <= 0:
        raise RuntimeError("covariance must be positive definite")
    return 0.5 * (logdet_a + logdet_b - logdet_s)


def entanglement_hidden_load(sigma: np.ndarray) -> dict[str, np.ndarray | float]:
    H = np.linalg.inv(sigma)
    phi = visible_precision(H, selector_a())
    ceiling = H[:2, :2]
    load = hidden_load(ceiling, phi, support_mode="ambient")
    return {
        "H": H,
        "phi": phi,
        "ceiling": ceiling,
        "lambda": load.reduced_lambda,
        "tau": load.clock,
        "mutual_information": gaussian_mutual_information(sigma),
    }


def _write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rs = np.linspace(0.0, 1.6, 33)
    vacuum_rows: list[dict[str, float]] = []
    thermal_rows: list[dict[str, float]] = []

    for r in rs:
        vacuum = entanglement_hidden_load(two_mode_squeezed_thermal_covariance(r))
        vacuum_rows.append(
            {
                "r": float(r),
                "tau": float(vacuum["tau"]),
                "two_mutual_information": float(2.0 * vacuum["mutual_information"]),
                "residual_abs": abs(float(vacuum["tau"] - 2.0 * vacuum["mutual_information"])),
            }
        )
        thermal = entanglement_hidden_load(two_mode_squeezed_thermal_covariance(r, n_a=0.4, n_b=0.9))
        thermal_rows.append(
            {
                "r": float(r),
                "tau": float(thermal["tau"]),
                "two_mutual_information": float(2.0 * thermal["mutual_information"]),
                "residual_abs": abs(float(thermal["tau"] - 2.0 * thermal["mutual_information"])),
            }
        )

    _write_csv(OUT / "tmsv_identity.csv", vacuum_rows)
    _write_csv(OUT / "thermal_identity.csv", thermal_rows)

    line_chart_svg(
        OUT / "tmsv_tau_vs_2mi.svg",
        rs,
        [
            ("tau", np.array([row["tau"] for row in vacuum_rows]), "#004488"),
            ("2 MI", np.array([row["two_mutual_information"] for row in vacuum_rows]), "#bb5566"),
        ],
        title="Two-Mode Squeezed Vacuum: tau and 2 MI",
        x_label="squeezing r",
        y_label="value",
    )
    line_chart_svg(
        OUT / "tmsv_residual.svg",
        rs,
        [("|tau - 2 MI|", np.array([row["residual_abs"] for row in vacuum_rows]), "#228833")],
        title="Two-Mode Squeezed Vacuum Residual",
        x_label="squeezing r",
        y_label="absolute residual",
    )

    sigma_1 = two_mode_squeezed_thermal_covariance(0.65)
    sigma_2 = two_mode_squeezed_thermal_covariance(0.85, n_a=0.3, n_b=0.4)
    load_1 = entanglement_hidden_load(sigma_1)
    load_2 = entanglement_hidden_load(sigma_2)
    contraction_pair = la.block_diag(hidden_contraction(load_1["lambda"]), hidden_contraction(load_2["lambda"]))
    load_composed = load_from_hidden_contraction(contraction_pair)
    additive_summary = {
        "tau_pair_1": float(load_1["tau"]),
        "tau_pair_2": float(load_2["tau"]),
        "tau_composed": float(clock(load_composed)),
        "composition_residual_abs": abs(float(clock(load_composed) - load_1["tau"] - load_2["tau"])),
    }
    headline_rows = [
        {
            "check": "tmsv_identity",
            "parameter_range": "r in [0.0, 1.6]",
            "value": max(row["residual_abs"] for row in vacuum_rows),
        },
        {
            "check": "thermal_identity",
            "parameter_range": "r in [0.0, 1.6], n_a = 0.4, n_b = 0.9",
            "value": max(row["residual_abs"] for row in thermal_rows),
        },
        {
            "check": "r_zero_tau",
            "parameter_range": "r = 0",
            "value": vacuum_rows[0]["tau"],
        },
        {
            "check": "r_zero_two_mutual_information",
            "parameter_range": "r = 0",
            "value": vacuum_rows[0]["two_mutual_information"],
        },
        {
            "check": "small_r_identity_residual",
            "parameter_range": "r = 0.05",
            "value": vacuum_rows[1]["residual_abs"],
        },
        {
            "check": "composition_residual",
            "parameter_range": "independent pair composition",
            "value": additive_summary["composition_residual_abs"],
        },
    ]
    _write_csv(OUT / "headline_summary.csv", headline_rows)

    summary = {
        "claim": "tau = 2 I(A:B) for the explicit Gaussian families in this example",
        "parameter_ranges": {
            "vacuum_squeezing_r": [float(rs[0]), float(rs[-1])],
            "thermal_squeezing_r": [float(rs[0]), float(rs[-1])],
        },
        "tmsv_max_residual_abs": max(row["residual_abs"] for row in vacuum_rows),
        "thermal_max_residual_abs": max(row["residual_abs"] for row in thermal_rows),
        "tmsv_monotone_tau": bool(
            np.all(np.diff(np.array([row["tau"] for row in vacuum_rows], dtype=float)) >= -1e-12)
        ),
        "additive_composition": additive_summary,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
