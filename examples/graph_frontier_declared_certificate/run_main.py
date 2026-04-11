from __future__ import annotations

import csv
import json
from math import sqrt
from pathlib import Path

import numpy as np

from nomogeo import declared_frontier_local_certificate, exact_branch_hessian, general_graph_frontier_hessian
from nomogeo.exceptions import InputValidationError

from examples.common import line_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
BASIS = np.array([[1.0], [0.0]])
WEIGHTS = np.array([0.5, 0.5])


def symmetric_family(e_abs: float) -> list[np.ndarray]:
    e = float(e_abs)
    return [
        np.array([[3.0, e], [e, 1.0]]),
        np.array([[3.0, -e], [-e, 1.0]]),
    ]


def _write_csv(path: Path, rows: list[dict[str, float | int | str | bool]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def exact_branch_applicability(e_abs: float) -> str:
    try:
        exact_branch_hessian(symmetric_family(e_abs), BASIS, weights=WEIGHTS)
    except InputValidationError:
        return "rejected_non_exact_branch"
    return "accepted_exact_branch"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    e_values = np.linspace(0.0, 2.4, 49)
    sweep_rows: list[dict[str, float | int | str | bool]] = []
    hessians: list[float] = []
    invalid_proxy: list[float] = []

    for e_abs in e_values:
        graph = general_graph_frontier_hessian(symmetric_family(float(e_abs)), BASIS, weights=WEIGHTS)
        hessian = float(graph.hessian_operator[0, 0])
        hessians.append(hessian)
        invalid_proxy.append(-24.0)
        sweep_rows.append(
            {
                "e_abs": float(e_abs),
                "graph_hessian": hessian,
                "closed_form_hessian": -24.0 + 8.0 * float(e_abs) ** 2,
                "invalid_exact_branch_proxy": -24.0,
                "stationarity_residual": graph.stationarity_residual,
                "status": graph.status,
                "exact_branch_applicability": exact_branch_applicability(float(e_abs)),
            }
        )

    _write_csv(OUT / "frontier_sweep.csv", sweep_rows)
    line_chart_svg(
        OUT / "hessian_breakpoint.svg",
        e_values,
        [
            ("graph Hessian", np.array(hessians), "#1d3557"),
            ("invalid exact-branch proxy", np.array(invalid_proxy), "#c1121f"),
        ],
        title="Declared Graph-Frontier Hessian Breakpoint",
        x_label="|e|",
        y_label="second variation",
    )

    key_cases = []
    for label, e_abs in [("stable_max", 1.0), ("degenerate", sqrt(3.0)), ("stable_min", 2.0)]:
        graph = general_graph_frontier_hessian(symmetric_family(e_abs), BASIS, weights=WEIGHTS)
        max_cert = declared_frontier_local_certificate(symmetric_family(e_abs), BASIS, weights=WEIGHTS, mode="max")
        min_cert = declared_frontier_local_certificate(symmetric_family(e_abs), BASIS, weights=WEIGHTS, mode="min")
        key_cases.append(
            {
                "label": label,
                "e_abs": float(e_abs),
                "graph_hessian": float(graph.hessian_operator[0, 0]),
                "closed_form_hessian": -24.0 + 8.0 * float(e_abs) ** 2,
                "status": graph.status,
                "stationarity_residual": graph.stationarity_residual,
                "max_certificate_passes": max_cert.certificate_passes,
                "max_certificate_kind": max_cert.certificate_kind,
                "min_certificate_passes": min_cert.certificate_passes,
                "min_certificate_kind": min_cert.certificate_kind,
                "exact_branch_applicability": exact_branch_applicability(e_abs),
            }
        )

    summary = {
        "claim": (
            "A declared observer can be stationary without being an exact branch. "
            "The general graph-frontier Hessian detects the sign change at |e| = sqrt(3), "
            "while the exact-branch proxy is not applicable outside the off-block-zero sector."
        ),
        "threshold_e_abs": sqrt(3.0),
        "key_cases": key_cases,
        "outputs": ["frontier_sweep.csv", "hessian_breakpoint.svg", "summary.json", "audit.json"],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
