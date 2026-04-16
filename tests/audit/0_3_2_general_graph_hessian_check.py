"""Research checks for the general graph-chart frontier Hessian.

The directional second-variation formula was already checked in
`0_3_2_near_branch_variation_check.py`. This script polarizes that formula
into a Hessian matrix, compares it against finite differences, and verifies
that it reduces to the module's exact-branch Hessian in the invariant sector.
It does not modify module source.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import exact_branch_hessian


OUT = Path("audit/outputs/0_3_2_general_graph_hessian_check.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def graph_basis(x: np.ndarray) -> np.ndarray:
    y = np.vstack([np.eye(x.shape[1]), x])
    q, _ = np.linalg.qr(y)
    return q


def frontier_score(family: list[np.ndarray], weights: np.ndarray, x: np.ndarray, mu: float) -> float:
    b = graph_basis(x)
    p = b @ b.T
    visible = 0.0
    leakage = 0.0
    for weight, operator in zip(weights, family):
        visible += float(weight * np.linalg.norm(p @ operator @ p, ord="fro") ** 2)
        commutator = operator @ p - p @ operator
        leakage += float(weight * 0.5 * np.linalg.norm(commutator, ord="fro") ** 2)
    return visible - mu * leakage


def directional_second_formula(family: list[np.ndarray], weights: np.ndarray, m: int, x: np.ndarray, mu: float) -> float:
    value = 0.0
    for weight, operator in zip(weights, family):
        a = operator[:m, :m]
        e = operator[m:, :m]
        d = operator[m:, m:]
        b1 = e.T @ x + x.T @ e
        s2 = (
            float(np.trace(b1 @ b1))
            + 2.0 * float(np.trace(a @ (x.T @ d @ x)))
            - 2.0 * float(np.trace((x.T @ x) @ (a @ a)))
        )
        moment = operator @ operator
        moment_u = moment[:m, :m]
        moment_w = moment[m:, m:]
        t2 = float(np.trace(x.T @ moment_w @ x) - np.trace((x.T @ x) @ moment_u))
        value += weight * ((1.0 + mu) * s2 - mu * t2)
    return 2.0 * value


def tangent_basis(hdim: int, m: int) -> list[np.ndarray]:
    basis = []
    for row in range(hdim):
        for col in range(m):
            x = np.zeros((hdim, m), dtype=float)
            x[row, col] = 1.0
            basis.append(x)
    return basis


def graph_hessian_formula(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float) -> np.ndarray:
    hdim = family[0].shape[0] - m
    basis = tangent_basis(hdim, m)
    dim = len(basis)
    hessian = np.zeros((dim, dim), dtype=float)
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            hessian[i, j] = (
                directional_second_formula(family, weights, m, ei + ej, mu)
                - directional_second_formula(family, weights, m, ei - ej, mu)
            ) / 4.0
    return sym(hessian)


def graph_hessian_finite(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, h: float = 2.0e-5) -> np.ndarray:
    hdim = family[0].shape[0] - m
    basis = tangent_basis(hdim, m)
    dim = len(basis)
    hessian = np.zeros((dim, dim), dtype=float)
    zero = np.zeros((hdim, m), dtype=float)
    f0 = frontier_score(family, weights, zero, mu)
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            if i == j:
                hessian[i, j] = (
                    frontier_score(family, weights, h * ei, mu)
                    - 2.0 * f0
                    + frontier_score(family, weights, -h * ei, mu)
                ) / (h * h)
            else:
                hessian[i, j] = (
                    frontier_score(family, weights, h * (ei + ej), mu)
                    - frontier_score(family, weights, h * (ei - ej), mu)
                    - frontier_score(family, weights, h * (-ei + ej), mu)
                    + frontier_score(family, weights, -h * (ei + ej), mu)
                ) / (4.0 * h * h)
    return sym(hessian)


def exact_branch_reduction(seed: int = 39301) -> dict[str, float | list[float]]:
    rng = np.random.default_rng(seed)
    n = 6
    m = 2
    hdim = n - m
    mu = 0.4
    weights = np.array([0.7, 1.2])
    family = []
    for _ in weights:
        a = sym(rng.normal(size=(m, m)))
        d = sym(rng.normal(size=(hdim, hdim)))
        family.append(np.block([[a, np.zeros((m, hdim))], [np.zeros((hdim, m)), d]]))

    formula = graph_hessian_formula(family, weights, m, mu)
    module = exact_branch_hessian(family, np.eye(n, m), weights=weights, mu=mu).second_variation_operator
    return {
        "max_abs_formula_minus_module": float(np.max(np.abs(formula - module))),
        "formula_eigenvalues": [float(x) for x in np.linalg.eigvalsh(formula)],
        "module_eigenvalues": [float(x) for x in np.linalg.eigvalsh(module)],
    }


def non_exact_finite_difference(seed: int = 39302) -> dict[str, float | list[float]]:
    rng = np.random.default_rng(seed)
    n = 7
    m = 3
    mu = 0.35
    weights = np.array([0.5, 1.1, 0.8])
    family = [sym(rng.normal(size=(n, n))) for _ in weights]
    formula = graph_hessian_formula(family, weights, m, mu)
    finite = graph_hessian_finite(family, weights, m, mu)
    return {
        "max_abs_formula_minus_finite": float(np.max(np.abs(formula - finite))),
        "relative_fro_error": float(np.linalg.norm(formula - finite, ord="fro") / max(1.0, np.linalg.norm(formula, ord="fro"))),
        "formula_eigenvalues": [float(x) for x in np.linalg.eigvalsh(formula)],
    }


def symmetric_pair_signs() -> list[dict[str, float]]:
    rows = []
    m = 1
    mu = 0.0
    weights = np.array([0.5, 0.5])
    for e in [1.0, np.sqrt(3.0), 2.0]:
        family = [
            np.array([[3.0, e], [e, 1.0]], dtype=float),
            np.array([[3.0, -e], [-e, 1.0]], dtype=float),
        ]
        hessian = graph_hessian_formula(family, weights, m, mu)
        proxy = graph_hessian_formula(
            [np.array([[3.0, 0.0], [0.0, 1.0]], dtype=float) for _ in range(2)],
            weights,
            m,
            mu,
        )
        rows.append(
            {
                "e_abs": float(e),
                "formula_hessian": float(hessian[0, 0]),
                "closed_form_hessian": float(-24.0 + 8.0 * e * e),
                "invalid_exact_branch_proxy_hessian": float(proxy[0, 0]),
            }
        )
    return rows


def main() -> None:
    results = {
        "exact_branch_reduction": exact_branch_reduction(),
        "non_exact_finite_difference": non_exact_finite_difference(),
        "symmetric_pair_signs": symmetric_pair_signs(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
