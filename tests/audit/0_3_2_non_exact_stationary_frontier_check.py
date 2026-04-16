"""Research checks for non-exact stationary weighted-frontier branches.

This script isolates a failure mode that is not covered by exact-branch
Hessian diagnostics: the weighted-frontier first variation can vanish by
cancellation even when no individual family member preserves the candidate
subspace. It does not modify module source.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


OUT = Path("audit/outputs/0_3_2_non_exact_stationary_frontier_check.json")


def projector_from_graph_scalar(x: float) -> np.ndarray:
    u = np.array([1.0, x], dtype=float)
    u /= np.linalg.norm(u)
    return np.outer(u, u)


def frontier_score_2d(rows: list[dict[str, float]], x: float, mu: float) -> float:
    p = projector_from_graph_scalar(x)
    score = 0.0
    total_weight = 0.0
    for row in rows:
        weight = float(row.get("weight", 1.0))
        a = float(row["a"])
        d = float(row["d"])
        e = float(row["e"])
        op = np.array([[a, e], [e, d]], dtype=float)
        visible = float(np.linalg.norm(p @ op @ p, ord="fro") ** 2)
        leakage = float(0.5 * np.linalg.norm(op @ p - p @ op, ord="fro") ** 2)
        score += weight * (visible - mu * leakage)
        total_weight += weight
    return score / total_weight


def gradient_formula_2d(rows: list[dict[str, float]], mu: float) -> float:
    total = 0.0
    weight_sum = 0.0
    for row in rows:
        weight = float(row.get("weight", 1.0))
        a = float(row["a"])
        d = float(row["d"])
        e = float(row["e"])
        total += weight * 2.0 * e * ((2.0 + mu) * a - mu * d)
        weight_sum += weight
    return total / weight_sum


def hessian_formula_2d(rows: list[dict[str, float]], mu: float) -> float:
    total = 0.0
    weight_sum = 0.0
    for row in rows:
        weight = float(row.get("weight", 1.0))
        a = float(row["a"])
        d = float(row["d"])
        e = float(row["e"])
        s2 = 4.0 * e * e + 2.0 * a * d - 2.0 * a * a
        t2 = d * d - a * a
        total += weight * 2.0 * ((1.0 + mu) * s2 - mu * t2)
        weight_sum += weight
    return total / weight_sum


def exact_branch_proxy_hessian_2d(rows: list[dict[str, float]], mu: float) -> float:
    """Return the second variation obtained by dropping all off-block terms.

    This is not a valid diagnostic when any e != 0. It is included only to
    quantify the mistake an approximate exact-branch extension would make.
    """

    proxy_rows = [{**row, "e": 0.0} for row in rows]
    return hessian_formula_2d(proxy_rows, mu)


def finite_difference_derivatives(rows: list[dict[str, float]], mu: float, h: float = 1.0e-5) -> dict[str, float]:
    f0 = frontier_score_2d(rows, 0.0, mu)
    fp = frontier_score_2d(rows, h, mu)
    fm = frontier_score_2d(rows, -h, mu)
    return {
        "gradient": float((fp - fm) / (2.0 * h)),
        "hessian": float((fp - 2.0 * f0 + fm) / (h * h)),
    }


def grid_verdict(rows: list[dict[str, float]], mu: float, radius: float = 4.0, points: int = 8001) -> dict[str, float]:
    xs = np.linspace(-radius, radius, points)
    vals = np.array([frontier_score_2d(rows, float(x), mu) for x in xs])
    best = int(np.argmax(vals))
    zero_index = int(np.argmin(np.abs(xs)))
    return {
        "radius": float(radius),
        "x_best": float(xs[best]),
        "score_at_zero": float(vals[zero_index]),
        "score_best": float(vals[best]),
        "score_gap_best_minus_zero": float(vals[best] - vals[zero_index]),
    }


def symmetric_pair_sweep() -> list[dict[str, float]]:
    rows = []
    a = 3.0
    d = 1.0
    mu = 0.0
    threshold = math.sqrt(a * (a - d) / 2.0)
    for e in [0.0, 0.5, 1.0, 1.5, 1.70, threshold, 1.75, 2.0, 3.0]:
        family = [{"a": a, "d": d, "e": e}, {"a": a, "d": d, "e": -e}]
        fd = finite_difference_derivatives(family, mu)
        grid = grid_verdict(family, mu)
        hess = hessian_formula_2d(family, mu)
        proxy = exact_branch_proxy_hessian_2d(family, mu)
        rows.append(
            {
                "a": a,
                "d": d,
                "e_abs": float(e),
                "mu": mu,
                "stationarity_gradient_formula": float(gradient_formula_2d(family, mu)),
                "hessian_formula": float(hess),
                "finite_gradient": fd["gradient"],
                "finite_hessian": fd["hessian"],
                "invalid_exact_branch_proxy_hessian": float(proxy),
                "hessian_minus_proxy": float(hess - proxy),
                "sharp_local_sign_threshold_e_abs": float(threshold),
                **grid,
            }
        )
    return rows


def penalty_cancellation_examples() -> list[dict[str, float]]:
    rows = []
    for a, d, e, mu in [(1.0, 3.0, 0.2, 1.0), (1.0, 3.0, 0.02, 1.0), (1.0, 2.0, 0.2, 2.0)]:
        family = [{"a": a, "d": d, "e": e}]
        fd = finite_difference_derivatives(family, mu)
        grid = grid_verdict(family, mu)
        rows.append(
            {
                "a": a,
                "d": d,
                "e": e,
                "mu": mu,
                "stationarity_gradient_formula": float(gradient_formula_2d(family, mu)),
                "hessian_formula": float(hessian_formula_2d(family, mu)),
                "finite_gradient": fd["gradient"],
                "finite_hessian": fd["hessian"],
                "invalid_exact_branch_proxy_hessian": float(exact_branch_proxy_hessian_2d(family, mu)),
                **grid,
            }
        )
    return rows


def random_cancellation_search(seed: int = 39201) -> dict[str, list[dict[str, float]]]:
    rng = np.random.default_rng(seed)
    local_max = []
    local_min = []
    near_degenerate = []
    mu = 0.0
    for _ in range(600):
        a1 = float(rng.uniform(0.5, 5.0))
        d1 = float(rng.uniform(0.1, 4.0))
        a2 = float(rng.uniform(0.5, 5.0))
        d2 = float(rng.uniform(0.1, 4.0))
        e1 = float(rng.choice([-1.0, 1.0]) * rng.uniform(0.05, 2.5))
        factor1 = (2.0 + mu) * a1 - mu * d1
        factor2 = (2.0 + mu) * a2 - mu * d2
        if abs(factor2) < 1.0e-12:
            continue
        e2 = -e1 * factor1 / factor2
        if abs(e2) > 3.5:
            continue
        family = [{"a": a1, "d": d1, "e": e1}, {"a": a2, "d": d2, "e": e2}]
        grad = gradient_formula_2d(family, mu)
        hess = hessian_formula_2d(family, mu)
        proxy = exact_branch_proxy_hessian_2d(family, mu)
        record = {
            "family": family,
            "mu": mu,
            "gradient_formula": float(grad),
            "hessian_formula": float(hess),
            "invalid_exact_branch_proxy_hessian": float(proxy),
            "hessian_minus_proxy": float(hess - proxy),
            **grid_verdict(family, mu, radius=4.0, points=2001),
        }
        if hess < -0.5 and len(local_max) < 4:
            local_max.append(record)
        if hess > 0.5 and len(local_min) < 4:
            local_min.append(record)
        if abs(hess) <= 0.05 and len(near_degenerate) < 4:
            near_degenerate.append(record)
        if len(local_max) >= 4 and len(local_min) >= 4 and len(near_degenerate) >= 4:
            break
    return {
        "stationary_non_exact_local_max_examples": local_max,
        "stationary_non_exact_local_min_examples": local_min,
        "stationary_non_exact_near_degenerate_examples": near_degenerate,
    }


def main() -> None:
    results = {
        "symmetric_pair_sweep": symmetric_pair_sweep(),
        "penalty_cancellation_examples": penalty_cancellation_examples(),
        "random_cancellation_search": random_cancellation_search(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
