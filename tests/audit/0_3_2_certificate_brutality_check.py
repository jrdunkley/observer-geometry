"""Brutality checks for the 0.3.2 near-branch certificate research.

Research-only script. It stress-checks projector derivative constants and the
practical non-vacuity of the conservative weighted-frontier near-branch
certificate. It does not modify module source.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


OUT = Path("audit/outputs/0_3_2_certificate_brutality_check.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def projector_constants(rho: float) -> tuple[float, float, float]:
    y = math.sqrt(1.0 + rho * rho)
    c1 = 2.0 * y + 2.0 * rho * y * y
    c2 = 2.0 + 8.0 * rho * y + (8.0 * rho * rho + 2.0) * y * y
    c3 = (24.0 * rho + 48.0 * rho**3) * y * y + 6.0 * (8.0 * rho * rho + 2.0) * y + 12.0 * rho
    return c1, c2, c3


def lcert_sharpish(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, rho: float) -> float:
    n = family[0].shape[0]
    a2 = sum(float(w) * float(np.linalg.norm(a, ord=2) ** 2) for w, a in zip(weights, family))
    g1 = (3.0 * (1.0 + mu) * math.sqrt(m) + mu * math.sqrt(n)) * a2
    g2 = 6.0 * (1.0 + mu) * a2
    g3 = 6.0 * (1.0 + mu) * a2
    c1, c2, c3 = projector_constants(rho)
    return g3 * c1**3 + 3.0 * g2 * c2 * c1 + g1 * c3


def lcert_coarse(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float) -> float:
    n = family[0].shape[0]
    a2 = sum(float(w) * float(np.linalg.norm(a, ord=2) ** 2) for w, a in zip(weights, family))
    g1 = (3.0 * (1.0 + mu) * math.sqrt(m) + mu * math.sqrt(n)) * a2
    g2 = 6.0 * (1.0 + mu) * a2
    g3 = 6.0 * (1.0 + mu) * a2
    return g3 * 8.0**3 + 3.0 * g2 * 64.0 * 8.0 + g1 * 512.0


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


def gradient_formula(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float) -> np.ndarray:
    hdim = family[0].shape[0] - m
    gradient = np.zeros((hdim, m), dtype=float)
    for weight, operator in zip(weights, family):
        a = operator[:m, :m]
        e = operator[m:, :m]
        d = operator[m:, m:]
        gradient += 2.0 * weight * ((2.0 + mu) * e @ a - mu * d @ e)
    return gradient


def hessian_matrix_finite(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, h: float = 2.0e-5) -> np.ndarray:
    hdim = family[0].shape[0] - m
    dim = hdim * m
    basis = []
    for row in range(hdim):
        for col in range(m):
            x = np.zeros((hdim, m), dtype=float)
            x[row, col] = 1.0
            basis.append(x)
    hess = np.zeros((dim, dim), dtype=float)
    f0 = frontier_score(family, weights, np.zeros((hdim, m), dtype=float), mu)
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            if i == j:
                hess[i, j] = (
                    frontier_score(family, weights, h * ei, mu)
                    - 2.0 * f0
                    + frontier_score(family, weights, -h * ei, mu)
                ) / (h * h)
            else:
                hess[i, j] = (
                    frontier_score(family, weights, h * (ei + ej), mu)
                    - frontier_score(family, weights, h * (ei - ej), mu)
                    - frontier_score(family, weights, h * (-ei + ej), mu)
                    + frontier_score(family, weights, -h * (ei + ej), mu)
                ) / (4.0 * h * h)
    return sym(hess)


def certificate_summary(family: list[np.ndarray], weights: np.ndarray, m: int, mu: float, rho: float) -> dict:
    grad = gradient_formula(family, weights, m, mu)
    eps = float(np.linalg.norm(grad, ord="fro"))
    hess = hessian_matrix_finite(family, weights, m, mu)
    eig = np.linalg.eigvalsh(hess)
    lambda_margin = float(max(0.0, -np.max(eig)))
    new_l = lcert_sharpish(family, weights, m, mu, rho)
    old_l = lcert_coarse(family, weights, m, mu)
    if lambda_margin > 0.0:
        left = 4.0 * eps / lambda_margin
        new_radius = min(rho, lambda_margin / (2.0 * new_l))
        old_radius = min(rho, lambda_margin / (2.0 * old_l))
        new_pass = bool(left < new_radius)
        old_pass = bool(left < old_radius)
        displacement_bound = 2.0 * eps / lambda_margin
    else:
        left = math.inf
        new_radius = 0.0
        old_radius = 0.0
        new_pass = False
        old_pass = False
        displacement_bound = math.inf
    return {
        "eps": eps,
        "lambda_margin": lambda_margin,
        "hessian_eigenvalues": [float(x) for x in eig],
        "new_L_cert": float(new_l),
        "old_L_cert": float(old_l),
        "L_ratio_old_over_new": float(old_l / new_l) if new_l > 0 else None,
        "left_4eps_over_lambda": float(left),
        "new_certificate_radius": float(new_radius),
        "old_certificate_radius": float(old_radius),
        "new_certificate_passes": new_pass,
        "old_certificate_passes": old_pass,
        "displacement_bound_2eps_over_lambda": float(displacement_bound),
    }


def projector_bound_stress(seed: int = 37002) -> dict:
    rng = np.random.default_rng(seed)
    rows = []
    violations = 0
    for rho in [0.02, 0.05, 0.1, 0.25, 0.5, 1.0]:
        c1, c2, c3 = projector_constants(rho)
        max1 = max2 = max3 = 0.0
        for _ in range(120):
            hdim = int(rng.integers(1, 8))
            m = int(rng.integers(1, 7))
            x = rng.normal(size=(hdim, m))
            norm = np.linalg.norm(x, ord="fro")
            if norm > 0:
                x = x / norm * rho * rng.random()
            z = rng.normal(size=(hdim, m))
            z /= np.linalg.norm(z, ord="fro")

            def p(t: float) -> np.ndarray:
                b = graph_basis(x + t * z)
                return b @ b.T

            h = 2.0e-4
            p0 = p(0.0)
            pp = p(h)
            pm = p(-h)
            p2p = p(2.0 * h)
            p2m = p(-2.0 * h)
            d1 = np.linalg.norm((pp - pm) / (2.0 * h), ord="fro")
            d2 = np.linalg.norm((pp - 2.0 * p0 + pm) / (h * h), ord="fro")
            d3 = np.linalg.norm((p2p - 2.0 * pp + 2.0 * pm - p2m) / (2.0 * h**3), ord="fro")
            max1 = max(max1, float(d1))
            max2 = max(max2, float(d2))
            max3 = max(max3, float(d3))
            if d1 > c1 * 1.02 or d2 > c2 * 1.02 or d3 > c3 * 1.02:
                violations += 1
        rows.append(
            {
                "rho": rho,
                "C1": c1,
                "C2": c2,
                "C3": c3,
                "max_directional_d1": max1,
                "max_directional_d2": max2,
                "max_directional_d3": max3,
            }
        )
    return {"rows": rows, "violations_with_2pct_slack": violations}


def synthetic_certificate_cases() -> dict:
    m = 2
    n = 4
    mu = 0.0
    rho = 1.0
    weights = np.array([1.0])
    base = np.diag([3.0, 2.6, 0.7, 0.5])
    rows = []
    for delta in [0.0, 1.0e-7, 3.0e-7, 1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4]:
        op = base.copy()
        op[2, 0] = delta
        op[0, 2] = delta
        op[3, 1] = -0.8 * delta
        op[1, 3] = -0.8 * delta
        summary = certificate_summary([op], weights, m, mu, rho)
        summary["delta"] = delta
        rows.append(summary)
    return {"rows": rows}


def random_near_branch_search(seed: int = 37003) -> dict:
    rng = np.random.default_rng(seed)
    found_nonvacuous = []
    found_vacuous = []
    m = 2
    n = 5
    mu = 0.2
    weights = np.array([1.0])
    rho = 1.0
    for _ in range(200):
        visible = np.sort(rng.uniform(2.0, 5.0, size=m))[::-1]
        hidden = np.sort(rng.uniform(0.1, 1.5, size=n - m))
        base = np.diag(np.concatenate([visible, hidden]))
        raw = rng.normal(size=(n - m, m))
        raw /= max(1.0, np.linalg.norm(raw, ord="fro"))
        delta = 10.0 ** rng.uniform(-8.0, -3.0)
        op = base.copy()
        op[m:, :m] = delta * raw
        op[:m, m:] = delta * raw.T
        summary = certificate_summary([op], weights, m, mu, rho)
        record = {
            "visible": visible.tolist(),
            "hidden": hidden.tolist(),
            "delta": delta,
            "eps": summary["eps"],
            "lambda_margin": summary["lambda_margin"],
            "new_certificate_passes": summary["new_certificate_passes"],
            "old_certificate_passes": summary["old_certificate_passes"],
            "left_4eps_over_lambda": summary["left_4eps_over_lambda"],
            "new_certificate_radius": summary["new_certificate_radius"],
            "L_ratio_old_over_new": summary["L_ratio_old_over_new"],
        }
        if summary["new_certificate_passes"] and len(found_nonvacuous) < 5:
            found_nonvacuous.append(record)
        if summary["lambda_margin"] > 0 and not summary["new_certificate_passes"] and len(found_vacuous) < 5:
            found_vacuous.append(record)
        if len(found_nonvacuous) >= 5 and len(found_vacuous) >= 5:
            break
    return {"nonvacuous_examples": found_nonvacuous, "vacuous_examples": found_vacuous}


def main() -> None:
    results = {
        "projector_bound_stress": projector_bound_stress(),
        "synthetic_certificate_cases": synthetic_certificate_cases(),
        "random_near_branch_search": random_near_branch_search(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
