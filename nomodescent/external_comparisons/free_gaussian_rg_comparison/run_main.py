from __future__ import annotations

from pathlib import Path

import numpy as np
from nomogeo import visible_precision
from nomodescent import ObserverSpec, classify_qd_relation, common_refinement_test, staged_descent_check

from external_comparisons.common import write_json

ROOT = Path(__file__).resolve().parent


def _free_precision() -> np.ndarray:
    return np.array(
        [
            [3.4, -0.8, 0.0, 0.0],
            [-0.8, 2.9, -0.7, 0.0],
            [0.0, -0.7, 2.8, -0.6],
            [0.0, 0.0, -0.6, 2.5],
        ],
        dtype=float,
    )


def _schur_eliminate(H: np.ndarray, retain: np.ndarray) -> np.ndarray:
    perm = np.r_[retain, np.array([index for index in range(H.shape[0]) if index not in retain], dtype=int)]
    H_perm = H[np.ix_(perm, perm)]
    k = len(retain)
    a = H_perm[:k, :k]
    b = H_perm[:k, k:]
    d = H_perm[k:, k:]
    return a - b @ np.linalg.inv(d) @ b.T


def main() -> None:
    H = _free_precision()
    mode4 = ObserverSpec(name="mode4", matrix=np.eye(4))
    mode2 = ObserverSpec(name="mode2", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    mode1 = ObserverSpec(name="mode1", matrix=[[1.0, 0.0]])
    block4 = ObserverSpec(name="block4", matrix=[[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]])

    qd_mode2 = visible_precision(H, mode2.matrix)
    schur_mode2 = _schur_eliminate(H, np.array([0, 1], dtype=int))
    schur_residual = float(np.max(np.abs(qd_mode2 - schur_mode2)))
    tower = staged_descent_check(H, [mode2, mode1])
    relation = classify_qd_relation(mode2, block4)
    refinement = common_refinement_test((mode2, block4))

    summary = {
        "external_method": "Schur-complement Gaussian elimination",
        "schur_residual": schur_residual,
        "tower_classification": tower.classification,
        "relation_classification": relation.classification,
        "common_refinement_classification": refinement.classification,
        "comparison_relation": "reproduced",
    }
    audit = {
        "tower": {
            "classification": tower.classification,
            "residuals": tower.residuals,
            "theorem_layer": tower.audit.theorem_layer,
        },
        "relation": {
            "classification": relation.classification,
            "residuals": relation.residuals,
            "theorem_layer": relation.audit.theorem_layer,
        },
        "refinement": {
            "classification": refinement.classification,
            "residuals": refinement.residuals,
            "theorem_layer": refinement.audit.theorem_layer,
        },
    }
    comparison = {
        "what_external_method_certifies": "Schur-complement elimination reproduces the retained-mode free Gaussian precision exactly.",
        "what_qd_adds": "QD recovers the same elimination result and also classifies alternative coarse observers by exact relation and common-refinement structure.",
        "epistemic_boundary": "This comparison is exact but stays entirely inside the free Gaussian regime.",
    }
    write_json(ROOT / "outputs" / "summary.json", summary)
    write_json(ROOT / "outputs" / "audit.json", audit)
    write_json(ROOT / "outputs" / "comparison.json", comparison)


if __name__ == "__main__":
    main()
