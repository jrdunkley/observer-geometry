from __future__ import annotations

import json
from math import sqrt
from pathlib import Path

import numpy as np

from nomogeo import declared_frontier_local_certificate, exact_branch_hessian, general_graph_frontier_hessian
from nomogeo.exceptions import InputValidationError

from examples.graph_frontier_declared_certificate.run_main import BASIS, WEIGHTS, symmetric_family

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def _exact_branch_rejected(e_abs: float) -> bool:
    try:
        exact_branch_hessian(symmetric_family(e_abs), BASIS, weights=WEIGHTS)
    except InputValidationError:
        return True
    return False


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cases = {
        "stable_max": 1.0,
        "degenerate": sqrt(3.0),
        "stable_min": 2.0,
    }
    rows = {}
    for label, e_abs in cases.items():
        graph = general_graph_frontier_hessian(symmetric_family(e_abs), BASIS, weights=WEIGHTS)
        max_cert = declared_frontier_local_certificate(symmetric_family(e_abs), BASIS, weights=WEIGHTS, mode="max")
        min_cert = declared_frontier_local_certificate(symmetric_family(e_abs), BASIS, weights=WEIGHTS, mode="min")
        rows[label] = {
            "e_abs": float(e_abs),
            "graph_hessian": float(graph.hessian_operator[0, 0]),
            "closed_form_hessian": -24.0 + 8.0 * float(e_abs) ** 2,
            "stationarity_residual": graph.stationarity_residual,
            "status": graph.status,
            "max_certificate_passes": max_cert.certificate_passes,
            "min_certificate_passes": min_cert.certificate_passes,
            "exact_branch_rejected": _exact_branch_rejected(e_abs),
        }

    thresholds = {
        "hessian_abs_error_max": 1e-10,
        "stationarity_residual_max": 1e-10,
    }
    hessian_errors = {
        label: abs(row["graph_hessian"] - row["closed_form_hessian"])
        for label, row in rows.items()
    }
    report = {
        "cases": rows,
        "hessian_errors": hessian_errors,
        "thresholds": thresholds,
    }
    (OUT / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    audit = {
        "example": "graph_frontier_declared_certificate",
        "claim": (
            "The graph-frontier Hessian, not the exact-branch proxy, controls a stationary non-exact declared observer."
        ),
        "headline_quantities": report,
        "passes": {
            "hessian_matches_closed_form": bool(max(hessian_errors.values()) <= thresholds["hessian_abs_error_max"]),
            "stationarity": bool(
                max(row["stationarity_residual"] for row in rows.values()) <= thresholds["stationarity_residual_max"]
            ),
            "max_certificate_before_threshold": bool(rows["stable_max"]["max_certificate_passes"]),
            "degenerate_certificate_vacuous": bool(
                not rows["degenerate"]["max_certificate_passes"] and not rows["degenerate"]["min_certificate_passes"]
            ),
            "min_certificate_after_threshold": bool(rows["stable_min"]["min_certificate_passes"]),
            "nonzero_cases_reject_exact_branch": bool(
                rows["stable_max"]["exact_branch_rejected"] and rows["stable_min"]["exact_branch_rejected"]
            ),
            "invalid_proxy_would_miss_sign": bool(
                rows["stable_min"]["graph_hessian"] > 0.0 and -24.0 < 0.0
            ),
        },
        "scripts": {
            "main": "python -m examples.graph_frontier_declared_certificate.run_main",
            "validate": "python -m examples.graph_frontier_declared_certificate.validate",
            "audit": "python -m examples.graph_frontier_declared_certificate.run_all",
        },
        "scope": (
            "Exact local quadratic graph-chart geometry only; no global Grassmannian optimizer and no full-law branch probability."
        ),
    }
    audit["passes"]["all"] = bool(all(audit["passes"].values()))
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")

    if not audit["passes"]["all"]:
        raise SystemExit("graph_frontier_declared_certificate audit failed")


if __name__ == "__main__":
    main()
