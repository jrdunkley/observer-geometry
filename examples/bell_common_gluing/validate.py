from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from examples.bell_common_gluing.run_main import (
    CLASS_LABELS,
    arcsine_margin,
    best_psd_completion,
    classify_family,
    family_covariances,
    linear_completion_residual,
    normalized_correlators,
    pair_precision_residuals,
    variance_gap,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    compatible = family_covariances(0.0, 0.6)
    variance_only = family_covariances(0.22, 0.6)
    correlator_only = family_covariances(0.0, 0.82)
    both = family_covariances(0.22, 0.82)

    compatible_sigma, _alpha, _beta, compatible_gap = best_psd_completion(0.6)
    _, _, _, correlator_gap = best_psd_completion(0.82)

    neighborhood = []
    for delta in np.linspace(0.18, 0.26, 5):
        for rho in np.linspace(0.56, 0.66, 5):
            neighborhood.append(classify_family(float(delta), float(rho)))

    report = {
        "compatible_class": CLASS_LABELS[classify_family(0.0, 0.6)],
        "variance_only_class": CLASS_LABELS[classify_family(0.22, 0.6)],
        "correlator_only_class": CLASS_LABELS[classify_family(0.0, 0.82)],
        "both_class": CLASS_LABELS[classify_family(0.22, 0.82)],
        "correlator_matrix_match_residual": float(
            np.linalg.norm(normalized_correlators(compatible) - normalized_correlators(variance_only), ord=np.inf)
        ),
        "compatible_variance_gap": variance_gap(compatible),
        "variance_only_variance_gap": variance_gap(variance_only),
        "compatible_arcsine_margin": arcsine_margin(compatible),
        "variance_only_arcsine_margin": arcsine_margin(variance_only),
        "compatible_linear_completion_residual": linear_completion_residual(compatible),
        "variance_only_linear_completion_residual": linear_completion_residual(variance_only),
        "correlator_only_linear_completion_residual": linear_completion_residual(correlator_only),
        "compatible_best_psd_gap": compatible_gap,
        "correlator_only_best_psd_gap": correlator_gap,
        "compatible_pair_precision_residual_max": max(pair_precision_residuals(compatible_sigma).values()),
        "variance_open_region_confirmed": bool(all(cls == 1 for cls in neighborhood)),
        "both_obstructions_present": bool(
            variance_gap(both) > 1e-12 and arcsine_margin(both) < -1e-12
        ),
        "open_region_grid": {
            "delta_range": [0.18, 0.26],
            "rho_range": [0.56, 0.66],
            "grid_shape": [5, 5],
        },
    }
    (OUT / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    thresholds = {
        "correlator_matrix_match_residual_max": 1e-12,
        "compatible_pair_precision_residual_max": 1e-10,
        "compatible_linear_completion_residual_max": 1e-10,
        "variance_only_linear_completion_residual_min": 1e-3,
        "correlator_only_best_psd_gap_max": -1e-3,
    }
    audit = {
        "example": "bell_common_gluing",
        "claim": (
            "There exists an explicit Gaussian family where the normalized correlator summary is unchanged, "
            "the arcsine-CHSH obstruction is absent, and common Gaussian gluing still fails because of a second obstruction."
        ),
        "parameter_ranges": {
            "sweep_delta": [0.0, 0.35],
            "sweep_rho": [0.0, 0.95],
            "confirmed_variance_only_neighborhood": report["open_region_grid"],
        },
        "headline_quantities": report,
        "thresholds": thresholds,
        "passes": {
            "correlator_invariance_between_compatible_and_variance_only": bool(
                report["correlator_matrix_match_residual"] <= thresholds["correlator_matrix_match_residual_max"]
            ),
            "compatible_reobservation": bool(
                report["compatible_pair_precision_residual_max"] <= thresholds["compatible_pair_precision_residual_max"]
            ),
            "compatible_linear_completion": bool(
                report["compatible_linear_completion_residual"] <= thresholds["compatible_linear_completion_residual_max"]
            ),
            "variance_only_separation": bool(
                report["variance_only_linear_completion_residual"] >= thresholds["variance_only_linear_completion_residual_min"]
            ),
            "correlator_only_psd_incompatibility": bool(
                report["correlator_only_best_psd_gap"] <= thresholds["correlator_only_best_psd_gap_max"]
            ),
            "open_region_confirmed": bool(report["variance_open_region_confirmed"]),
            "both_obstructions_present": bool(report["both_obstructions_present"]),
        },
        "scripts": {
            "main": "python -m examples.bell_common_gluing.run_main",
            "validate": "python -m examples.bell_common_gluing.validate",
            "audit": "python -m examples.bell_common_gluing.run_all",
        },
    }
    audit["passes"]["all"] = bool(all(audit["passes"].values()))
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
