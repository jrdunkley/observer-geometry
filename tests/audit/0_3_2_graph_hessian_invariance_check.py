"""Frame-invariance checks for the general graph-chart frontier Hessian.

Research-only script. It verifies that the graph-Hessian formula is a
geometric declared-observer object, not an artifact of the coordinate split
used in the first derivation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import exact_branch_hessian


OUT = Path("audit/outputs/0_3_2_graph_hessian_invariance_check.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def orthogonal(rng: np.random.Generator, n: int) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(n, n)))
    signs = np.sign(np.diag(r))
    signs[signs == 0.0] = 1.0
    return q @ np.diag(signs)


def blocks(operator: np.ndarray, b: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return sym(b.T @ operator @ b), w.T @ operator @ b, sym(w.T @ operator @ w)


def directional_second_blocks(
    block_family: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    weights: np.ndarray,
    x: np.ndarray,
    mu: float,
) -> float:
    value = 0.0
    for weight, (a, e, d) in zip(weights, block_family):
        b1 = e.T @ x + x.T @ e
        s2 = (
            float(np.trace(b1 @ b1))
            + 2.0 * float(np.trace(a @ (x.T @ d @ x)))
            - 2.0 * float(np.trace((x.T @ x) @ (a @ a)))
        )
        moment_u = a @ a + e.T @ e
        moment_w = d @ d + e @ e.T
        t2 = float(np.trace(x.T @ moment_w @ x) - np.trace((x.T @ x) @ moment_u))
        value += weight * ((1.0 + mu) * s2 - mu * t2)
    return 2.0 * value


def bilinear_second_blocks(
    block_family: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    weights: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    mu: float,
) -> float:
    return 0.25 * (
        directional_second_blocks(block_family, weights, x + y, mu)
        - directional_second_blocks(block_family, weights, x - y, mu)
    )


def gradient_blocks(
    block_family: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    weights: np.ndarray,
    mu: float,
) -> np.ndarray:
    hdim, m = block_family[0][1].shape
    gradient = np.zeros((hdim, m), dtype=float)
    for weight, (a, e, d) in zip(weights, block_family):
        gradient += 2.0 * weight * ((2.0 + mu) * e @ a - mu * d @ e)
    return gradient


def tangent_basis(hdim: int, m: int) -> list[np.ndarray]:
    basis = []
    for row in range(hdim):
        for col in range(m):
            x = np.zeros((hdim, m), dtype=float)
            x[row, col] = 1.0
            basis.append(x)
    return basis


def hessian_matrix_blocks(
    block_family: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    weights: np.ndarray,
    mu: float,
) -> np.ndarray:
    hdim, m = block_family[0][1].shape
    basis = tangent_basis(hdim, m)
    hessian = np.zeros((len(basis), len(basis)), dtype=float)
    for col, x in enumerate(basis):
        for row, y in enumerate(basis):
            hessian[row, col] = bilinear_second_blocks(block_family, weights, y, x, mu)
    return sym(hessian)


def frontier_score_frame(
    family: list[np.ndarray],
    weights: np.ndarray,
    b: np.ndarray,
    w: np.ndarray,
    x: np.ndarray,
    mu: float,
) -> float:
    raw = b + w @ x
    q, _ = np.linalg.qr(raw)
    p = q @ q.T
    score = 0.0
    for weight, operator in zip(weights, family):
        visible = float(np.linalg.norm(p @ operator @ p, ord="fro") ** 2)
        commutator = operator @ p - p @ operator
        leakage = float(0.5 * np.linalg.norm(commutator, ord="fro") ** 2)
        score += weight * (visible - mu * leakage)
    return score


def finite_bilinear(
    family: list[np.ndarray],
    weights: np.ndarray,
    b: np.ndarray,
    w: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    mu: float,
    h: float = 2.0e-5,
) -> float:
    return (
        frontier_score_frame(family, weights, b, w, h * (x + y), mu)
        - frontier_score_frame(family, weights, b, w, h * (x - y), mu)
        - frontier_score_frame(family, weights, b, w, h * (-x + y), mu)
        + frontier_score_frame(family, weights, b, w, -h * (x + y), mu)
    ) / (4.0 * h * h)


def finite_gradient(
    family: list[np.ndarray],
    weights: np.ndarray,
    b: np.ndarray,
    w: np.ndarray,
    hdim: int,
    m: int,
    mu: float,
    h: float = 1.0e-6,
) -> np.ndarray:
    grad = np.zeros((hdim, m), dtype=float)
    zero = np.zeros((hdim, m), dtype=float)
    for row in range(hdim):
        for col in range(m):
            x = np.zeros_like(zero)
            x[row, col] = 1.0
            grad[row, col] = (
                frontier_score_frame(family, weights, b, w, h * x, mu)
                - frontier_score_frame(family, weights, b, w, -h * x, mu)
            ) / (2.0 * h)
    return grad


def random_frame_finite_difference(seed: int = 39401) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = 8
    m = 3
    hdim = n - m
    mu = 0.45
    q = orthogonal(rng, n)
    b = q[:, :m]
    w = q[:, m:]
    weights = np.array([0.4, 1.0, 0.7])
    family = [sym(rng.normal(size=(n, n))) for _ in weights]
    block_family = [blocks(operator, b, w) for operator in family]
    grad_formula = gradient_blocks(block_family, weights, mu)
    grad_finite = finite_gradient(family, weights, b, w, hdim, m, mu)

    max_abs = 0.0
    rels = []
    for _ in range(25):
        x = rng.normal(size=(hdim, m))
        y = rng.normal(size=(hdim, m))
        formula = bilinear_second_blocks(block_family, weights, x, y, mu)
        finite = finite_bilinear(family, weights, b, w, x, y, mu)
        max_abs = max(max_abs, abs(formula - finite))
        rels.append(abs(formula - finite) / max(1.0, abs(formula)))
    return {
        "max_abs_gradient_formula_minus_finite": float(np.max(np.abs(grad_formula - grad_finite))),
        "gradient_relative_fro_error": float(
            np.linalg.norm(grad_formula - grad_finite, ord="fro") / max(1.0, np.linalg.norm(grad_formula, ord="fro"))
        ),
        "max_abs_formula_minus_finite": float(max_abs),
        "max_relative_error": float(max(rels)),
    }


def basis_change_invariance(seed: int = 39402) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = 9
    m = 4
    hdim = n - m
    mu = 0.25
    q = orthogonal(rng, n)
    b = q[:, :m]
    w = q[:, m:]
    ru = orthogonal(rng, m)
    rw = orthogonal(rng, hdim)
    b2 = b @ ru
    w2 = w @ rw
    weights = np.array([0.6, 1.4])
    family = [sym(rng.normal(size=(n, n))) for _ in weights]
    blocks1 = [blocks(operator, b, w) for operator in family]
    blocks2 = [blocks(operator, b2, w2) for operator in family]

    max_abs = 0.0
    for _ in range(60):
        x = rng.normal(size=(hdim, m))
        y = rng.normal(size=(hdim, m))
        # Coordinates transform so that W X B^T = W2 X2 B2^T.
        x2 = rw.T @ x @ ru
        y2 = rw.T @ y @ ru
        h1 = bilinear_second_blocks(blocks1, weights, x, y, mu)
        h2 = bilinear_second_blocks(blocks2, weights, x2, y2, mu)
        max_abs = max(max_abs, abs(h1 - h2))
    return {"max_abs_bilinear_invariance_error": float(max_abs)}


def exact_branch_arbitrary_frame_reduction(seed: int = 39403) -> dict[str, float | list[float]]:
    rng = np.random.default_rng(seed)
    n = 7
    m = 3
    hdim = n - m
    mu = 0.55
    q = orthogonal(rng, n)
    b = q[:, :m]
    w = q[:, m:]
    weights = np.array([0.3, 1.1, 0.9])
    family = []
    block_family = []
    for _ in weights:
        a = sym(rng.normal(size=(m, m)))
        d = sym(rng.normal(size=(hdim, hdim)))
        block = np.block([[a, np.zeros((m, hdim))], [np.zeros((hdim, m)), d]])
        s = np.column_stack([b, w])
        family.append(sym(s @ block @ s.T))
        block_family.append((a, np.zeros((hdim, m)), d))
    formula_in_declared_w = hessian_matrix_blocks(block_family, weights, mu)
    result = exact_branch_hessian(family, b, weights=weights, mu=mu)
    module = result.second_variation_operator
    module_w = result.complement_basis
    block_family_module_w = [blocks(operator, b, module_w) for operator in family]
    formula_in_module_w = hessian_matrix_blocks(block_family_module_w, weights, mu)
    eig_formula_declared_w = np.linalg.eigvalsh(formula_in_declared_w)
    eig_module = np.linalg.eigvalsh(module)
    return {
        "raw_matrix_error_in_different_complement_frames": float(np.max(np.abs(formula_in_declared_w - module))),
        "max_abs_formula_minus_module_in_module_complement_frame": float(np.max(np.abs(formula_in_module_w - module))),
        "max_abs_eigenvalue_error": float(np.max(np.abs(eig_formula_declared_w - eig_module))),
        "formula_eigenvalues": [float(x) for x in eig_formula_declared_w],
        "module_eigenvalues": [float(x) for x in eig_module],
    }


def main() -> None:
    results = {
        "random_frame_finite_difference": random_frame_finite_difference(),
        "basis_change_invariance": basis_change_invariance(),
        "exact_branch_arbitrary_frame_reduction": exact_branch_arbitrary_frame_reduction(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
