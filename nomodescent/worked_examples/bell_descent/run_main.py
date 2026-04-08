from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from nomodescent import (
    AssumptionEntry,
    AssumptionLedger,
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    common_descent_test,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def observers() -> tuple[ObserverSpec, ...]:
    return (
        ObserverSpec(name="00", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="01", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
        ObserverSpec(name="10", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="11", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )


def family_covariances(delta: float, rho: float) -> dict[str, np.ndarray]:
    return {
        "00": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
        "01": np.array([[1.0 + delta, rho * np.sqrt(1.0 + delta)], [rho * np.sqrt(1.0 + delta), 1.0]], dtype=float),
        "10": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
        "11": np.array([[1.0, -rho], [-rho, 1.0]], dtype=float),
    }


def correlator_matrix(family: dict[str, np.ndarray]) -> np.ndarray:
    return np.array(
        [
            [family["00"][0, 1] / np.sqrt(family["00"][0, 0] * family["00"][1, 1]), family["01"][0, 1] / np.sqrt(family["01"][0, 0] * family["01"][1, 1])],
            [family["10"][0, 1] / np.sqrt(family["10"][0, 0] * family["10"][1, 1]), family["11"][0, 1] / np.sqrt(family["11"][0, 0] * family["11"][1, 1])],
        ]
    )


def build_problem(name: str, delta: float, rho: float) -> ProblemSpec:
    family = family_covariances(delta, rho)
    evidence = tuple(
        VisibleEvidenceSpec(name=f"sigma_{key}", observer=key, kind="covariance", matrix=matrix)
        for key, matrix in family.items()
    )
    assumptions = AssumptionLedger(
        entries=(
            AssumptionEntry(label="gaussian", statement="common latent object is Gaussian", exact=True),
            AssumptionEntry(label="finite_dim", statement="latent dimension fixed to four Bell variables", exact=True),
        )
    )
    return ProblemSpec(
        name=name,
        latent_dim=4,
        observers=observers(),
        evidence=evidence,
        assumptions=assumptions,
        goals=(GoalSpec(kind="common_completion"),),
        description="Bell square as a common Gaussian descent problem",
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sample_specs = {
        "compatible": (0.0, 0.6, True),
        "variance_only": (0.22, 0.6, False),
        "correlator_only": (0.0, 0.82, True),
    }

    rows: list[dict[str, object]] = []
    audits: dict[str, object] = {}
    correlators: dict[str, np.ndarray] = {}
    for name, (delta, rho, allow_approx) in sample_specs.items():
        problem = build_problem(name, delta, rho)
        result = common_descent_test(problem, allow_approximate_psd_search=allow_approx, psd_search_grid_size=31)
        family = family_covariances(delta, rho)
        correlators[name] = correlator_matrix(family)
        certificate_kinds = ",".join(cert.kind for cert in result.certificates) if result.certificates else ""
        rows.append(
            {
                "sample": name,
                "delta": delta,
                "rho": rho,
                "classification": result.classification,
                "exact": result.exact,
                "correlator_matches_compatible": bool(np.allclose(correlators[name], correlators["compatible"])) if name != "compatible" else True,
                "linear_residual": result.residuals.get("linear_residual", 0.0),
                "psd_margin_or_search_margin": result.residuals.get("psd_margin", result.residuals.get("best_psd_margin", float("nan"))),
                "certificate_kinds": certificate_kinds,
            }
        )
        audits[name] = {
            "classification": result.classification,
            "exact": result.exact,
            "audit": {
                "theorem_layer": result.audit.theorem_layer,
                "residuals": result.audit.residuals,
                "falsification_route": list(result.audit.falsification_route),
            },
            "certificates": [
                {"kind": cert.kind, "exact": cert.exact, "summary": cert.summary, "details": cert.details}
                for cert in result.certificates
            ],
        }

    _write_csv(OUT / "comparison_table.csv", rows)
    summary = {
        "compatible_vs_variance_only_correlator_residual": float(
            np.linalg.norm(correlators["compatible"] - correlators["variance_only"], ord=np.inf)
        ),
        "sample_classifications": {row["sample"]: row["classification"] for row in rows},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audits, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
