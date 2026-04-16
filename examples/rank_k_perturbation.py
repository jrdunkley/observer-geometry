from __future__ import annotations

import json

import numpy as np

from nomogeo import rank_k_covariance_perturbation


def main() -> None:
    rng = np.random.default_rng(3201)
    rows = []
    for k in (1, 2, 3):
        n = 7 + k
        m = 3
        raw = rng.normal(size=(n, n))
        sigma0 = raw.T @ raw / n + np.eye(n)
        factor = rng.normal(size=(n, k)) / 4.0
        result = rank_k_covariance_perturbation(sigma0, factor, m)
        rows.append(
            {
                "k": k,
                "visible_dim": m,
                "ambient_dim": n,
                "update_rank": result.update_rank,
                "rank_bound": result.rank_bound,
                "formula_residual": result.formula_residual,
                "singular_values": result.singular_values.tolist(),
            }
        )
    print(json.dumps({"rank_k_covariance_perturbation": rows}, indent=2))


if __name__ == "__main__":
    main()
