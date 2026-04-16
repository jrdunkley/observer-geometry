"""Follow-up numerics for the Claude notes audit.

This script tests claims that are useful as research inspiration but should not
be bulk-imported into the papers:

1. rank-one differential correlations are exactly one-channel in the white
   covariance sector;
2. that one-channel conclusion is not stable under a generic coloured
   covariance background;
3. generic noncommuting perturbation families need not admit exact
   closure-adapted observers of a prescribed intermediate rank;
4. branch or observer selection is robust only with a residual margin;
5. a distribution over local quadratic geometries can carry tail information
   invisible to a single averaged Hessian.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np


OUT = Path("audit/outputs/claude_notes_followup_numerics.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def make_spd(rng: np.random.Generator, n: int, floor: float = 0.7) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a / n + floor * np.eye(n)


def visible_precision_from_cov(sigma: np.ndarray, m: int) -> np.ndarray:
    return np.linalg.inv(sigma[:m, :m])


def hidden_gap_from_precision(h: np.ndarray, sigma: np.ndarray, m: int) -> np.ndarray:
    return sym(h[:m, :m] - visible_precision_from_cov(sigma, m))


def numerical_rank(a: np.ndarray, tol: float = 1.0e-9) -> int:
    return int(np.sum(np.linalg.svd(a, compute_uv=False) > tol))


def rank_one_white_sector(seed: int = 9001) -> dict:
    rng = np.random.default_rng(seed)
    n, m = 9, 4
    eps = 0.37
    f = rng.normal(size=n)
    a = f[:m]
    b = f[m:]
    sigma = np.eye(n) + eps * np.outer(f, f)
    h = np.linalg.inv(sigma)
    gap = hidden_gap_from_precision(h, sigma, m)

    A = float(a @ a)
    B = float(b @ b)
    coeff = eps * eps * B / ((1.0 + eps * A) * (1.0 + eps * (A + B)))
    predicted = coeff * np.outer(a, a)
    fi = float(f @ h @ f)
    fi_pred = float((f @ f) / (1.0 + eps * (f @ f)))
    return {
        "n": n,
        "m": m,
        "eps": eps,
        "gap_formula_max_abs_error": float(np.max(np.abs(gap - predicted))),
        "gap_rank": numerical_rank(gap),
        "gap_eigenvalues": np.linalg.eigvalsh(gap).tolist(),
        "fisher_information": fi,
        "fisher_information_formula": fi_pred,
        "fisher_ceiling": 1.0 / eps,
    }


def rank_one_coloured_background(seed: int = 9002) -> dict:
    rng = np.random.default_rng(seed)
    n, m = 10, 5
    eps = 0.41

    best = None
    for attempt in range(200):
        sigma0 = make_spd(rng, n, floor=0.9)
        f = rng.normal(size=n)
        sigma1 = sigma0 + eps * np.outer(f, f)
        h0 = np.linalg.inv(sigma0)
        h1 = np.linalg.inv(sigma1)
        dg = hidden_gap_from_precision(h1, sigma1, m) - hidden_gap_from_precision(h0, sigma0, m)
        phi0 = visible_precision_from_cov(sigma0, m)
        u = h0 @ f
        u_v = u[:m]
        w_v = phi0 @ f[:m]
        beta = eps / (1.0 + eps * float(f @ h0 @ f))
        gamma = eps / (1.0 + eps * float(f[:m] @ phi0 @ f[:m]))
        predicted = gamma * np.outer(w_v, w_v) - beta * np.outer(u_v, u_v)
        s = np.linalg.svd(dg, compute_uv=False)
        r = int(np.sum(s > 1.0e-8))
        if best is None or r > best["rank"] or (r == best["rank"] and s[1] > best["singular_values"][1]):
            denom = float(np.linalg.norm(u_v) * np.linalg.norm(w_v))
            alignment = abs(float(u_v @ w_v)) / denom if denom > 0 else 0.0
            best = {
                "attempt": attempt,
                "rank": r,
                "singular_values": s.tolist(),
                "fro_norm": float(np.linalg.norm(dg, "fro")),
                "two_rank_formula_max_abs_error": float(np.max(np.abs(dg - predicted))),
                "u_v_w_v_alignment": alignment,
            }
        if r >= 3:
            break

    assert best is not None
    return best


def leakage(Ds: list[np.ndarray], B: np.ndarray) -> float:
    Pi = B @ B.T
    I = np.eye(Pi.shape[0])
    return float(sum(np.linalg.norm((I - Pi) @ D @ B, "fro") ** 2 for D in Ds))


def random_stiefel(rng: np.random.Generator, n: int, m: int) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(n, m)))
    signs = np.sign(np.diag(r))
    signs[signs == 0.0] = 1.0
    return q * signs


def noncommuting_closure_obstruction(seed: int = 9003) -> dict:
    rng = np.random.default_rng(seed)
    n, m = 6, 2
    D1 = np.diag(np.arange(1.0, n + 1.0))
    D2 = sym(rng.normal(size=(n, n)))
    Ds = [D1, D2]

    # Because D1 has simple spectrum, every D1-invariant m-plane is a coordinate
    # span. Exact common invariance therefore requires a zero cross-block of D2
    # for one coordinate subset.
    subset_cross_norms = []
    for subset in itertools.combinations(range(n), m):
        mask = np.zeros(n, dtype=bool)
        mask[list(subset)] = True
        block = D2[np.ix_(~mask, mask)]
        subset_cross_norms.append(float(np.linalg.norm(block, "fro")))
    min_coordinate_cross_norm = min(subset_cross_norms)

    best = math.inf
    for _ in range(30000):
        B = random_stiefel(rng, n, m)
        best = min(best, leakage(Ds, B))

    comm_norm = float(np.linalg.norm(D1 @ D2 - D2 @ D1, "fro"))
    return {
        "n": n,
        "m": m,
        "commutator_norm": comm_norm,
        "min_coordinate_cross_norm_certificate": min_coordinate_cross_norm,
        "random_search_best_leakage": best,
        "exact_common_rank_m_subspace_exists": bool(min_coordinate_cross_norm < 1.0e-10),
        "logic": "D1 has simple spectrum, so exact closure for D1 forces a coordinate eigenspan; all D2 cross-blocks are nonzero.",
    }


def robustness_margin(seed: int = 9004) -> dict:
    rng = np.random.default_rng(seed)
    gap = 0.08
    safe_R = 0.03
    unsafe_R = 0.05
    samples = 50000

    safe_reversals = 0
    unsafe_reversals = 0
    for _ in range(samples):
        e1, e2 = rng.uniform(-safe_R, safe_R, size=2)
        safe_reversals += int((0.0 + e1) > (gap + e2))
        u1, u2 = rng.uniform(-unsafe_R, unsafe_R, size=2)
        unsafe_reversals += int((0.0 + u1) > (gap + u2))

    constructive_reversal = {
        "branch_1_residual": unsafe_R,
        "branch_2_residual": -unsafe_R,
        "reversed_gap": float((gap - unsafe_R) - (0.0 + unsafe_R)),
    }
    return {
        "quadratic_gap": gap,
        "safe_residual_bound": safe_R,
        "safe_condition_gap_gt_2R": bool(gap > 2.0 * safe_R),
        "safe_reversal_count": safe_reversals,
        "unsafe_residual_bound": unsafe_R,
        "unsafe_condition_gap_gt_2R": bool(gap > 2.0 * unsafe_R),
        "unsafe_reversal_count": unsafe_reversals,
        "samples": samples,
        "constructive_reversal": constructive_reversal,
    }


def local_geometry_distribution(seed: int = 9005) -> dict:
    rng = np.random.default_rng(seed)
    n, m = 8, 4
    base = np.eye(n)
    coupling = np.zeros((n, n))
    coupling[:m, m:] = rng.normal(size=(m, n - m))
    coupling[m:, :m] = coupling[:m, m:].T
    coupling = coupling / np.linalg.norm(coupling, 2)

    def clock_for_strength(strength: float) -> float:
        h = base + strength * coupling + 0.8 * strength * strength * np.eye(n)
        sigma = np.linalg.inv(h)
        phi = visible_precision_from_cov(sigma, m)
        T = h[:m, :m]
        return float(np.linalg.slogdet(T)[1] - np.linalg.slogdet(phi)[1])

    strengths_lognormal = rng.lognormal(mean=-2.4, sigma=0.9, size=4000)
    strengths_two_point = np.where(rng.random(4000) < 0.95, 0.04, 0.52)
    # Match means by rescaling the two-point sample.
    strengths_two_point *= float(np.mean(strengths_lognormal) / np.mean(strengths_two_point))

    clocks_ln = np.array([clock_for_strength(x) for x in strengths_lognormal])
    clocks_tp = np.array([clock_for_strength(x) for x in strengths_two_point])
    avg_strength = float(np.mean(strengths_lognormal))
    clock_at_mean = clock_for_strength(avg_strength)

    return {
        "mean_strength_lognormal": avg_strength,
        "mean_strength_two_point": float(np.mean(strengths_two_point)),
        "clock_at_mean_strength": clock_at_mean,
        "mean_clock_lognormal": float(np.mean(clocks_ln)),
        "mean_clock_two_point": float(np.mean(clocks_tp)),
        "p99_clock_lognormal": float(np.quantile(clocks_ln, 0.99)),
        "p99_clock_two_point": float(np.quantile(clocks_tp, 0.99)),
        "clock_kurtosis_lognormal": float(np.mean((clocks_ln - np.mean(clocks_ln)) ** 4) / np.var(clocks_ln) ** 2),
        "clock_kurtosis_two_point": float(np.mean((clocks_tp - np.mean(clocks_tp)) ** 4) / np.var(clocks_tp) ** 2),
    }


def main() -> None:
    results = {
        "rank_one_white_sector": rank_one_white_sector(),
        "rank_one_coloured_background": rank_one_coloured_background(),
        "noncommuting_closure_obstruction": noncommuting_closure_obstruction(),
        "robustness_margin": robustness_margin(),
        "local_geometry_distribution": local_geometry_distribution(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
