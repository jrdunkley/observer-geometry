"""Second-pass probes for the 0.3.2 research plan.

These tests stress the plan itself. They focus on places where a research
target could be over-promoted: near-branch Hessians, variable-rank frontier
selection, and rank-k covariance perturbation generalization.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import weighted_family_frontier_scores


OUT = Path("audit/outputs/0_3_2_second_pass_probe.json")


def _sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def _rank(a: np.ndarray, tol: float = 1.0e-9) -> int:
    return int(np.sum(np.linalg.svd(a, compute_uv=False) > tol))


def _graph_basis(x: np.ndarray) -> np.ndarray:
    graph = np.vstack([np.eye(x.shape[1]), x])
    q, _ = np.linalg.qr(graph)
    return q


def _frontier_score(family: list[np.ndarray], weights: np.ndarray, x: np.ndarray, mu: float) -> float:
    b = _graph_basis(x)
    p = b @ b.T
    visible = 0.0
    leakage = 0.0
    for weight, op in zip(weights, family):
        visible += float(weight * np.linalg.norm(p @ op @ p, ord="fro") ** 2)
        comm = op @ p - p @ op
        leakage += float(weight * 0.5 * np.linalg.norm(comm, ord="fro") ** 2)
    return visible - mu * leakage


def near_branch_stationarity_probe() -> dict:
    """A nearly exact branch is generally not a stationary branch."""
    n, m = 5, 2
    hdim = n - m
    weights = np.array([1.0])
    mu = 0.5
    base = np.diag([1.0, 1.8, 2.2, 2.9, 3.4])
    direction = np.zeros((n, n))
    direction[0, 2] = 1.0
    direction[1, 4] = -0.7
    direction = _sym(direction)
    h = 2.0e-6

    rows = []
    for eps in [0.0, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2]:
        op = base + eps * direction
        gradient = np.zeros((hdim, m))
        for r in range(hdim):
            for c in range(m):
                x = np.zeros((hdim, m))
                x[r, c] = 1.0
                gradient[r, c] = (
                    _frontier_score([op], weights, h * x, mu)
                    - _frontier_score([op], weights, -h * x, mu)
                ) / (2.0 * h)
        off_block = float(np.linalg.norm(op[m:, :m], ord="fro"))
        rows.append(
            {
                "epsilon": eps,
                "off_block_norm": off_block,
                "finite_difference_gradient_norm": float(np.linalg.norm(gradient, ord="fro")),
                "gradient_over_off_block": None if off_block == 0.0 else float(np.linalg.norm(gradient, ord="fro") / off_block),
            }
        )
    return {
        "rows": rows,
        "conclusion": "near-exact branches have first-order stationarity drift; a Hessian-like diagnostic alone is incomplete",
    }


def variable_rank_frontier_probe() -> dict:
    """If full rank is allowed and no dimension cost is charged, full rank wins."""
    rng = np.random.default_rng(33001)
    n = 7
    family = []
    for _ in range(4):
        a = rng.normal(size=(n, n))
        family.append(_sym(a))
    weights = np.array([0.4, 1.1, 0.9, 0.3])
    mu = 0.7
    rows = []
    full = weighted_family_frontier_scores(family, np.eye(n), weights=weights, mu=mu)
    moment_trace = float(np.trace(full.moment_operator))
    for m in range(1, n + 1):
        # Use the top eigenspace of the moment as a strong but declared ladder.
        evals, evecs = np.linalg.eigh(full.moment_operator)
        b = evecs[:, np.argsort(evals)[::-1][:m]]
        result = weighted_family_frontier_scores(family, b, weights=weights, mu=mu)
        rows.append(
            {
                "rank": m,
                "penalized": float(result.penalized_score),
                "visible": float(result.visible_score),
                "leakage": float(result.leakage),
                "captured_curvature": float(result.captured_curvature),
            }
        )

    cost_rows = []
    for dim_cost in [0.0, 0.25 * moment_trace / n, 0.5 * moment_trace / n, moment_trace / n]:
        best = max(rows, key=lambda row: row["penalized"] - dim_cost * row["rank"])
        cost_rows.append(
            {
                "dimension_cost": float(dim_cost),
                "best_rank": int(best["rank"]),
                "best_costed_score": float(best["penalized"] - dim_cost * best["rank"]),
            }
        )
    return {
        "mu": mu,
        "moment_trace_full_score": moment_trace,
        "rank_rows": rows,
        "costed_selection": cost_rows,
        "conclusion": "variable-rank observer selection needs a dimension budget or cost; otherwise full rank is structurally privileged",
    }


def rank_k_covariance_probe() -> dict:
    """Woodbury extension: rank-k covariance update gives rank at most 2k gap update."""
    rng = np.random.default_rng(33002)
    n, m, k = 9, 4, 2
    a = rng.normal(size=(n, n))
    sigma0 = _sym(a.T @ a / n + 1.1 * np.eye(n))
    f = rng.normal(size=(n, k))
    sigma1 = _sym(sigma0 + f @ f.T)
    h0 = np.linalg.inv(sigma0)
    h1 = np.linalg.inv(sigma1)
    phi0 = np.linalg.inv(sigma0[:m, :m])
    phi1 = np.linalg.inv(sigma1[:m, :m])
    direct = _sym((h1[:m, :m] - phi1) - (h0[:m, :m] - phi0))

    u = h0 @ f
    u_v = u[:m, :]
    w_v = phi0 @ f[:m, :]
    beta = np.linalg.inv(np.eye(k) + f.T @ h0 @ f)
    gamma = np.linalg.inv(np.eye(k) + f[:m, :].T @ phi0 @ f[:m, :])
    formula = _sym(w_v @ gamma @ w_v.T - u_v @ beta @ u_v.T)

    hidden_only = np.zeros((n, k))
    hidden_only[m:, :] = rng.normal(size=(n - m, k))
    sigma_h = _sym(sigma0 + hidden_only @ hidden_only.T)
    h_h = np.linalg.inv(sigma_h)
    phi_h = np.linalg.inv(sigma_h[:m, :m])
    hidden_direct = _sym((h_h[:m, :m] - phi_h) - (h0[:m, :m] - phi0))

    block_sigma0 = np.block(
        [
            [np.eye(m), np.zeros((m, n - m))],
            [np.zeros((n - m, m)), 1.3 * np.eye(n - m)],
        ]
    )
    sigma_b = _sym(block_sigma0 + hidden_only @ hidden_only.T)
    block_direct = _sym(
        (np.linalg.inv(sigma_b)[:m, :m] - np.linalg.inv(sigma_b[:m, :m]))
        - (np.linalg.inv(block_sigma0)[:m, :m] - np.linalg.inv(block_sigma0[:m, :m]))
    )

    return {
        "k": k,
        "generic_rank": _rank(direct),
        "generic_rank_bound_2k": 2 * k,
        "formula_residual_max_abs": float(np.max(np.abs(direct - formula))),
        "generic_singular_values": [float(x) for x in np.linalg.svd(direct, compute_uv=False)],
        "hidden_only_correlated_rank": _rank(hidden_direct),
        "hidden_only_block_diagonal_rank": _rank(block_direct),
        "conclusion": "rank-k covariance updates appear to obey the expected rank <= 2k formula, with hidden-only nonzero effects depending on baseline coupling",
    }


def main() -> None:
    results = {
        "near_branch_stationarity": near_branch_stationarity_probe(),
        "variable_rank_frontier": variable_rank_frontier_probe(),
        "rank_k_covariance": rank_k_covariance_probe(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
