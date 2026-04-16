"""Numerical checks for near-branch first and second variation formulas.

Research-only audit script. It validates the graph-chart formulas recorded in
`audit/0_3_2_theory_research_note.md` and does not modify module source.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


OUT = Path("audit/outputs/0_3_2_near_branch_variation_check.json")


def _sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def _graph_basis(x: np.ndarray) -> np.ndarray:
    y = np.vstack([np.eye(x.shape[1]), x])
    q, _ = np.linalg.qr(y)
    return q


def _score(family: list[np.ndarray], weights: np.ndarray, x: np.ndarray, mu: float) -> float:
    b = _graph_basis(x)
    p = b @ b.T
    visible = 0.0
    leakage = 0.0
    for weight, operator in zip(weights, family):
        visible += float(weight * np.linalg.norm(p @ operator @ p, ord="fro") ** 2)
        commutator = operator @ p - p @ operator
        leakage += float(weight * 0.5 * np.linalg.norm(commutator, ord="fro") ** 2)
    return visible - mu * leakage


def _gradient_formula(family: list[np.ndarray], weights: np.ndarray, visible_dim: int, mu: float) -> np.ndarray:
    m = visible_dim
    hdim = family[0].shape[0] - m
    gradient = np.zeros((hdim, m), dtype=float)
    for weight, operator in zip(weights, family):
        a = operator[:m, :m]
        e = operator[m:, :m]
        d = operator[m:, m:]
        gradient += 2.0 * weight * ((2.0 + mu) * e @ a - mu * d @ e)
    return gradient


def _second_formula(
    family: list[np.ndarray],
    weights: np.ndarray,
    visible_dim: int,
    x: np.ndarray,
    mu: float,
) -> float:
    m = visible_dim
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


def main() -> None:
    rng = np.random.default_rng(35001)
    n, m = 7, 3
    hdim = n - m
    mu = 0.6
    weights = np.array([0.5, 1.3, 0.7])
    family = [_sym(rng.normal(size=(n, n))) for _ in weights]

    gradient = _gradient_formula(family, weights, m, mu)
    grad_h = 1.0e-6
    finite_gradient = np.zeros_like(gradient)
    for row in range(hdim):
        for col in range(m):
            direction = np.zeros((hdim, m), dtype=float)
            direction[row, col] = 1.0
            finite_gradient[row, col] = (
                _score(family, weights, grad_h * direction, mu)
                - _score(family, weights, -grad_h * direction, mu)
            ) / (2.0 * grad_h)

    x = rng.normal(size=(hdim, m))
    x /= np.linalg.norm(x)
    second_formula = _second_formula(family, weights, m, x, mu)
    second_rows = []
    for h in [1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5]:
        finite_second = (
            _score(family, weights, h * x, mu)
            - 2.0 * _score(family, weights, np.zeros_like(x), mu)
            + _score(family, weights, -h * x, mu)
        ) / (h * h)
        second_rows.append(
            {
                "h": h,
                "finite_difference": float(finite_second),
                "formula": float(second_formula),
                "absolute_error": float(abs(finite_second - second_formula)),
                "relative_error": float(abs(finite_second - second_formula) / max(1.0, abs(second_formula))),
            }
        )

    results = {
        "gradient": {
            "formula_norm": float(np.linalg.norm(gradient)),
            "finite_difference_norm": float(np.linalg.norm(finite_gradient)),
            "max_abs_error": float(np.max(np.abs(gradient - finite_gradient))),
            "relative_error": float(np.linalg.norm(gradient - finite_gradient) / max(1.0, np.linalg.norm(gradient))),
        },
        "second_variation": second_rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
