from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.linalg as la

from nomogeo import clock, hidden_contraction, hidden_load, load_from_hidden_contraction, visible_precision

from examples.entanglement_hidden_load.run_main import (
    entanglement_hidden_load,
    gaussian_mutual_information,
    local_rotation,
    local_squeezer,
    two_mode_squeezed_thermal_covariance,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rs = np.linspace(0.0, 1.5, 16)

    max_residual = 0.0
    max_local_symplectic_residual = 0.0
    for r in rs:
        sigma = two_mode_squeezed_thermal_covariance(r, n_a=0.2, n_b=0.6)
        base = entanglement_hidden_load(sigma)
        max_residual = max(max_residual, abs(float(base["tau"] - 2.0 * base["mutual_information"])))

        local = local_rotation(0.23, -0.41) @ local_squeezer(0.37, -0.18)
        transformed_sigma = local @ sigma @ local.T
        transformed = entanglement_hidden_load(transformed_sigma)
        max_local_symplectic_residual = max(
            max_local_symplectic_residual,
            abs(float(base["tau"] - transformed["tau"])),
            abs(float(gaussian_mutual_information(sigma) - gaussian_mutual_information(transformed_sigma))),
        )

    sigma_1 = two_mode_squeezed_thermal_covariance(0.55)
    sigma_2 = two_mode_squeezed_thermal_covariance(0.95, n_a=0.1, n_b=0.7)
    s_global = la.block_diag(sigma_1, sigma_2)
    h_global = np.linalg.inv(s_global)
    c_global = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        ]
    )
    phi_global = visible_precision(h_global, c_global)
    ceiling_global = h_global[np.ix_([0, 1, 4, 5], [0, 1, 4, 5])]
    global_load = hidden_load(ceiling_global, phi_global, support_mode="ambient")

    pair_1 = entanglement_hidden_load(sigma_1)
    pair_2 = entanglement_hidden_load(sigma_2)
    k_global = la.block_diag(hidden_contraction(pair_1["lambda"]), hidden_contraction(pair_2["lambda"]))
    lambda_from_factors = load_from_hidden_contraction(k_global)

    report = {
        "max_abs_residual_tau_minus_2mi": max_residual,
        "max_local_symplectic_invariance_residual": max_local_symplectic_residual,
        "global_clock": float(global_load.clock),
        "sum_pair_clocks": float(pair_1["tau"] + pair_2["tau"]),
        "composition_clock_residual": abs(float(global_load.clock - pair_1["tau"] - pair_2["tau"])),
        "contraction_recovery_residual": float(np.linalg.norm(global_load.reduced_lambda - lambda_from_factors, ord="fro")),
    }
    (OUT / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    thresholds = {
        "identity_residual_abs_max": 1e-10,
        "local_symplectic_invariance_residual_max": 1e-10,
        "composition_clock_residual_max": 1e-10,
        "contraction_recovery_residual_max": 1e-10,
    }
    audit = {
        "example": "entanglement_hidden_load",
        "claim": "For the stated Gaussian families, the nomogeo hidden-load clock equals twice the standard Gaussian mutual information.",
        "parameter_ranges": {
            "validation_squeezing_r": [float(rs[0]), float(rs[-1])],
            "thermal_occupancies": {"n_a": 0.2, "n_b": 0.6},
        },
        "headline_quantities": report,
        "thresholds": thresholds,
        "passes": {
            "identity": bool(report["max_abs_residual_tau_minus_2mi"] <= thresholds["identity_residual_abs_max"]),
            "local_symplectic_invariance": bool(
                report["max_local_symplectic_invariance_residual"] <= thresholds["local_symplectic_invariance_residual_max"]
            ),
            "composition": bool(report["composition_clock_residual"] <= thresholds["composition_clock_residual_max"]),
            "contraction_recovery": bool(
                report["contraction_recovery_residual"] <= thresholds["contraction_recovery_residual_max"]
            ),
        },
        "scripts": {
            "main": "python -m examples.entanglement_hidden_load.run_main",
            "validate": "python -m examples.entanglement_hidden_load.validate",
            "audit": "python -m examples.entanglement_hidden_load.run_all",
        },
    }
    audit["passes"]["all"] = bool(all(audit["passes"].values()))
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
