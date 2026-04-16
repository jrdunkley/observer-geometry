"""Stress checks for a declared-observer local frontier certificate.

This research-only script combines the general graph gradient/Hessian with the
conservative Hessian-Lipschitz bound from the projector-bounds note. It checks
that the resulting certificate applies beyond exact branches, especially to
stationary non-exact observers.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


OUT = Path("audit/outputs/0_3_2_declared_frontier_certificate_check.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def projector_constants(rho: float) -> tuple[float, float, float]:
    y = math.sqrt(1.0 + rho * rho)
    return (
        2.0 * y + 2.0 * rho * y * y,
        2.0 + 8.0 * rho * y + (8.0 * rho * rho + 2.0) * y * y,
        (24.0 * rho + 48.0 * rho**3) * y * y + 6.0 * (8.0 * rho * rho + 2.0) * y + 12.0 * rho,
    )


def lcert(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, rho: float) -> float:
    n = family[0].shape[0]
    a2 = sum(float(w) * float(np.linalg.norm(a, ord=2) ** 2) for w, a in zip(weights, family))
    g1 = (3.0 * (1.0 + mu) * math.sqrt(m) + mu * math.sqrt(n)) * a2
    g2 = 6.0 * (1.0 + mu) * a2
    g3 = 6.0 * (1.0 + mu) * a2
    c1, c2, c3 = projector_constants(rho)
    return g3 * c1**3 + 3.0 * g2 * c2 * c1 + g1 * c3


def blocks(operator: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return operator[:m, :m], operator[m:, :m], operator[m:, m:]


def gradient(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float) -> np.ndarray:
    hdim = family[0].shape[0] - m
    grad = np.zeros((hdim, m), dtype=float)
    for weight, operator in zip(weights, family):
        u, e, d = blocks(operator, m)
        grad += 2.0 * weight * ((2.0 + mu) * e @ u - mu * d @ e)
    return grad


def directional_second(family: list[np.ndarray], weights: np.ndarray, m: int, x: np.ndarray, mu: float) -> float:
    value = 0.0
    for weight, operator in zip(weights, family):
        u, e, d = blocks(operator, m)
        b = e.T @ x + x.T @ e
        s2 = (
            float(np.trace(b @ b))
            + 2.0 * float(np.trace(u @ (x.T @ d @ x)))
            - 2.0 * float(np.trace((x.T @ x) @ (u @ u)))
        )
        moment = operator @ operator
        t2 = float(np.trace(x.T @ moment[m:, m:] @ x) - np.trace((x.T @ x) @ moment[:m, :m]))
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


def hessian(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float) -> np.ndarray:
    hdim = family[0].shape[0] - m
    basis = tangent_basis(hdim, m)
    out = np.zeros((len(basis), len(basis)), dtype=float)
    for row, x in enumerate(basis):
        for col, y in enumerate(basis):
            out[row, col] = 0.25 * (
                directional_second(family, weights, m, x + y, mu)
                - directional_second(family, weights, m, x - y, mu)
            )
    return sym(out)


def certificate(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, rho: float, mode: str) -> dict[str, float | bool | list[float]]:
    grad = gradient(family, weights, m, mu)
    hess = hessian(family, weights, m, mu)
    eigs = np.linalg.eigvalsh(hess)
    eps = float(np.linalg.norm(grad, ord="fro"))
    if mode == "max":
        lambda_margin = float(max(0.0, -np.max(eigs)))
    elif mode == "min":
        lambda_margin = float(max(0.0, np.min(eigs)))
    else:
        raise ValueError("mode must be 'max' or 'min'")
    lipschitz = lcert(family, weights, m, mu, rho)
    if lambda_margin > 0.0:
        r0 = min(rho, lambda_margin / (2.0 * lipschitz))
        left = 4.0 * eps / lambda_margin
        passes = bool(left < r0)
        displacement = 2.0 * eps / lambda_margin
    else:
        r0 = 0.0
        left = None
        passes = False
        displacement = None
    return {
        "mode": mode,
        "eps": eps,
        "hessian_eigenvalues": [float(x) for x in eigs],
        "lambda_margin": lambda_margin,
        "L_cert": float(lipschitz),
        "r0": float(r0),
        "left_4eps_over_lambda": None if left is None else float(left),
        "displacement_bound_2eps_over_lambda": None if displacement is None else float(displacement),
        "certificate_passes": passes,
    }


def symmetric_non_exact_stationary_cases() -> list[dict[str, object]]:
    rows = []
    weights = np.array([0.5, 0.5])
    m = 1
    mu = 0.0
    rho = 1.0
    for e in [1.0, 1.7, math.sqrt(3.0), 1.75, 2.0]:
        family = [
            np.array([[3.0, e], [e, 1.0]], dtype=float),
            np.array([[3.0, -e], [-e, 1.0]], dtype=float),
        ]
        rows.append(
            {
                "e_abs": e,
                "max_certificate": certificate(family, weights, m, mu, rho, mode="max"),
                "min_certificate": certificate(family, weights, m, mu, rho, mode="min"),
            }
        )
    return rows


def near_stationary_imbalance_cases() -> list[dict[str, object]]:
    rows = []
    weights = np.array([0.5, 0.5])
    m = 1
    mu = 0.0
    rho = 1.0
    for imbalance in [0.0, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3]:
        e_plus = 1.0
        e_minus = -(1.0 - imbalance)
        family = [
            np.array([[3.0, e_plus], [e_plus, 1.0]], dtype=float),
            np.array([[3.0, e_minus], [e_minus, 1.0]], dtype=float),
        ]
        rows.append(
            {
                "imbalance": imbalance,
                "e_plus": e_plus,
                "e_minus": e_minus,
                "max_certificate": certificate(family, weights, m, mu, rho, mode="max"),
            }
        )
    return rows


def random_declared_cases(seed: int = 39501) -> dict[str, list[dict[str, object]]]:
    rng = np.random.default_rng(seed)
    found_max = []
    found_vacuous = []
    m = 2
    n = 5
    weights = np.array([1.0])
    mu = 0.2
    rho = 1.0
    for _ in range(500):
        visible = np.sort(rng.uniform(2.0, 5.0, size=m))[::-1]
        hidden = np.sort(rng.uniform(0.1, 1.5, size=n - m))
        base = np.diag(np.concatenate([visible, hidden]))
        raw = rng.normal(size=(n - m, m))
        raw /= max(1.0, np.linalg.norm(raw, ord="fro"))
        delta = 10.0 ** rng.uniform(-8.0, -3.0)
        op = base.copy()
        op[m:, :m] = delta * raw
        op[:m, m:] = delta * raw.T
        cert = certificate([op], weights, m, mu, rho, mode="max")
        record = {
            "visible": visible.tolist(),
            "hidden": hidden.tolist(),
            "delta": delta,
            "certificate": cert,
        }
        if cert["certificate_passes"] and len(found_max) < 4:
            found_max.append(record)
        if cert["lambda_margin"] > 0 and not cert["certificate_passes"] and len(found_vacuous) < 4:
            found_vacuous.append(record)
        if len(found_max) >= 4 and len(found_vacuous) >= 4:
            break
    return {"nonvacuous_max_examples": found_max, "vacuous_max_examples": found_vacuous}


def main() -> None:
    results = {
        "symmetric_non_exact_stationary_cases": symmetric_non_exact_stationary_cases(),
        "near_stationary_imbalance_cases": near_stationary_imbalance_cases(),
        "random_declared_cases": random_declared_cases(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
