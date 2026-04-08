from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from nomodescent import (
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    classify_qd_relation,
    common_descent_test,
    minimal_refinement_search,
)

from external_comparisons.common import write_json

ROOT = Path(__file__).resolve().parent
SOURCE_ROWS = ROOT.parents[2] / "evidence" / "micro_real_bundles" / "iris_protocol_mismatch" / "raw" / "iris_slice.csv"


def _load_rows() -> np.ndarray:
    rows: list[list[float]] = []
    with SOURCE_ROWS.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append([float(row[key]) for key in ("sepal_length", "sepal_width", "petal_length", "petal_width")])
    return np.asarray(rows, dtype=float)


def _first_canonical_correlation(covariance: np.ndarray) -> float:
    s_xx = covariance[:2, :2]
    s_yy = covariance[2:, 2:]
    s_xy = covariance[:2, 2:]
    whitening_x = _inverse_sqrt_psd(s_xx)
    whitening_y = _inverse_sqrt_psd(s_yy)
    whitened = whitening_x @ s_xy @ whitening_y
    singular_values = np.linalg.svd(whitened, compute_uv=False)
    return float(min(1.0, max(0.0, singular_values[0])))


def _inverse_sqrt_psd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    clipped = np.clip(eigenvalues, 1e-12, None)
    return eigenvectors @ np.diag(1.0 / np.sqrt(clipped)) @ eigenvectors.T


def main() -> None:
    rows = _load_rows()
    covariance = np.cov(rows, rowvar=False, bias=False)
    sepal = ObserverSpec(name="sepal_panel", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    petal = ObserverSpec(name="petal_panel", matrix=[[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
    problem = ProblemSpec(
        name="iris_protocol_mismatch",
        latent_dim=4,
        observers=(sepal, petal),
        evidence=(
            VisibleEvidenceSpec(name="cov_sepal", observer="sepal_panel", kind="covariance", matrix=covariance[:2, :2]),
            VisibleEvidenceSpec(name="cov_petal", observer="petal_panel", kind="covariance", matrix=covariance[2:, 2:]),
        ),
        goals=(GoalSpec(kind="common_completion"),),
    )

    cca_correlation = _first_canonical_correlation(covariance)
    relation = classify_qd_relation(sepal, petal)
    completion = common_descent_test(problem)
    refinement = minimal_refinement_search(
        np.linalg.inv(covariance),
        [sepal, petal],
        [ObserverSpec(name="full_iris", matrix=np.eye(4))],
    )

    summary = {
        "external_method": "canonical correlation analysis",
        "data_source": str(SOURCE_ROWS.relative_to(ROOT.parents[2])),
        "first_canonical_correlation": cca_correlation,
        "qd_relation_classification": relation.classification,
        "qd_completion_classification": completion.classification,
        "qd_refinement_winner": refinement.winner,
        "comparison_relation": "distinguished",
    }
    audit = {
        "relation": {
            "classification": relation.classification,
            "residuals": relation.residuals,
            "theorem_layer": relation.audit.theorem_layer,
        },
        "completion": {
            "classification": completion.classification,
            "residuals": completion.residuals,
            "certificates": [certificate.kind for certificate in completion.certificates],
            "theorem_layer": completion.audit.theorem_layer,
        },
        "refinement": {
            "winner": refinement.winner,
            "scores": refinement.scores,
            "theorem_layer": refinement.audit.theorem_layer,
        },
    }
    comparison = {
        "what_external_method_certifies": "CCA finds a strong cross-panel dependence but does not classify observer relation or decide common completion.",
        "what_qd_adds": "QD certifies that the two panels are exactly non-nested observers and that the current evidence remains exactly underdetermined rather than resolved.",
        "epistemic_boundary": "The covariances are exact sample covariances of the tiny real row slice, so the result is a tiny-bundle comparison rather than a broad empirical claim.",
    }
    write_json(ROOT / "outputs" / "summary.json", summary)
    write_json(ROOT / "outputs" / "audit.json", audit)
    write_json(ROOT / "outputs" / "comparison.json", comparison)


if __name__ == "__main__":
    main()

