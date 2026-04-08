from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from nomogeo import gaussian_bhattacharyya_distance, observed_covariance, visible_precision
from nomodescent import ObserverSpec, classify_relation, minimal_refinement_search

from worked_examples.common import bar_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def latent_precision() -> np.ndarray:
    return np.array(
        [
            [3.3, 0.9, 0.6, 0.0],
            [0.9, 2.9, 0.4, 0.1],
            [0.6, 0.4, 2.2, 0.0],
            [0.0, 0.1, 0.0, 1.7],
        ],
        dtype=float,
    )


def observers() -> tuple[ObserverSpec, ObserverSpec, list[ObserverSpec]]:
    protocol_a = ObserverSpec(name="protocol_a", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    protocol_b = ObserverSpec(name="protocol_b", matrix=[[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    candidates = [
        ObserverSpec(name="panel_3", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="full_4", matrix=np.eye(4)),
    ]
    return protocol_a, protocol_b, candidates


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    H = latent_precision()
    protocol_a, protocol_b, candidates = observers()

    cov_a = observed_covariance(H, protocol_a.matrix)
    cov_b = observed_covariance(H, protocol_b.matrix)
    phi_a = visible_precision(H, protocol_a.matrix)
    phi_b = visible_precision(H, protocol_b.matrix)
    mismatch = gaussian_bhattacharyya_distance(cov_a, cov_b)

    relation = classify_relation(protocol_a, protocol_b)
    search = minimal_refinement_search(H, [protocol_a, protocol_b], candidates, objective="min_dimension_then_logdet")

    rows = [
        {
            "protocol": "A",
            "visible_logdet": float(np.linalg.slogdet(phi_a)[1]),
            "cov_trace": float(np.trace(cov_a)),
        },
        {
            "protocol": "B",
            "visible_logdet": float(np.linalg.slogdet(phi_b)[1]),
            "cov_trace": float(np.trace(cov_b)),
        },
    ]
    _write_csv(OUT / "protocol_summary.csv", rows)
    _write_csv(
        OUT / "candidate_scores.csv",
        [{"candidate": name, "score": score} for name, score in search.scores.items()],
    )

    bar_chart_svg(
        OUT / "protocol_logdet.svg",
        ["protocol_A", "protocol_B"],
        np.array([np.linalg.slogdet(phi_a)[1], np.linalg.slogdet(phi_b)[1]]),
        title="Observer Mismatch In Visible Precision",
        y_label="log det visible precision",
        color="#2a9d8f",
    )

    summary = {
        "protocol_relation": relation.classification,
        "bhattacharyya_mismatch_score": float(mismatch),
        "refinement_winner": search.winner,
        "refinement_scores": search.scores,
    }
    audit = {
        "theorem_layer": "exact observer relation tests and exact finite-family refinement search",
        "residuals": {
            "protocol_a_through_b_residual": relation.residuals["a_through_b_residual"],
            "protocol_b_through_a_residual": relation.residuals["b_through_a_residual"],
        },
        "falsification_route": list(search.audit.falsification_route),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
