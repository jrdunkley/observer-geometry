from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np

from nomogeo import canonical_lift, dv_bridge, hidden_projector, local_visible_calculus, visible_precision


OUT = Path(__file__).resolve().parent / "closure_adapted_outputs"
OUT.mkdir(exist_ok=True)

RNG = np.random.default_rng(20260409)


def sym(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def make_spd(dim: int, cond_target: float = 10.0) -> np.ndarray:
    Q, _ = np.linalg.qr(RNG.standard_normal((dim, dim)))
    eigs = np.geomspace(1.0, cond_target, dim)
    return sym(Q @ np.diag(eigs) @ Q.T)


def random_stiefel(n: int, m: int) -> np.ndarray:
    Q, _ = np.linalg.qr(RNG.standard_normal((n, m)))
    return Q[:, :m]


def sqrt_spd(H: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eigh(sym(H))
    return sym(V @ np.diag(np.sqrt(w)) @ V.T)


def inv_sqrt_spd(H: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eigh(sym(H))
    return sym(V @ np.diag(1.0 / np.sqrt(w)) @ V.T)


def observer_from_B(H: np.ndarray, B: np.ndarray) -> np.ndarray:
    return B.T @ sqrt_spd(H)


def whitened_delta(H: np.ndarray, Delta: np.ndarray) -> np.ndarray:
    return sym(inv_sqrt_spd(H) @ Delta @ inv_sqrt_spd(H))


def closure_scores(H: np.ndarray, family: list[np.ndarray], B: np.ndarray) -> dict[str, float]:
    s_total = 0.0
    h_total = 0.0
    projector = np.eye(B.shape[0]) - B @ B.T
    for Delta in family:
        W = whitened_delta(H, Delta)
        V = B.T @ W @ B
        Q = B.T @ W @ projector @ W @ B
        s_total += float(np.linalg.norm(V, ord="fro") ** 2)
        h_total += 2.0 * float(np.trace(Q))
    eta = h_total / (s_total + h_total) if (s_total + h_total) > 0.0 else 0.0
    return {"S": s_total, "H": h_total, "eta": eta}


def theorem_normal_form_audit(num_trials: int = 50) -> dict[str, float]:
    phi_errors = []
    lift_errors = []
    proj_errors = []
    v_errors = []
    q_errors = []
    for _ in range(num_trials):
        n = int(RNG.integers(4, 10))
        m = int(RNG.integers(2, n))
        H = make_spd(n, cond_target=float(10.0 ** RNG.uniform(0.0, 4.0)))
        B = random_stiefel(n, m)
        C = observer_from_B(H, B)
        Delta = sym(RNG.standard_normal((n, n)))
        W = whitened_delta(H, Delta)

        phi = visible_precision(H, C)
        lift = canonical_lift(H, C)
        projector = hidden_projector(H, C)
        local = local_visible_calculus(H, C, Delta)

        phi_errors.append(float(np.linalg.norm(phi - np.eye(m), ord="fro")))
        lift_errors.append(float(np.linalg.norm(lift - inv_sqrt_spd(H) @ B, ord="fro")))
        proj_expected = np.eye(n) - inv_sqrt_spd(H) @ B @ B.T @ sqrt_spd(H)
        proj_errors.append(float(np.linalg.norm(projector - proj_expected, ord="fro")))
        v_errors.append(float(np.linalg.norm(local.V - B.T @ W @ B, ord="fro")))
        q_expected = B.T @ W @ (np.eye(n) - B @ B.T) @ W @ B
        q_errors.append(float(np.linalg.norm(local.Q - q_expected, ord="fro")))

    return {
        "num_trials": num_trials,
        "max_phi_error": float(max(phi_errors)),
        "max_lift_error": float(max(lift_errors)),
        "max_projector_error": float(max(proj_errors)),
        "max_V_error": float(max(v_errors)),
        "max_Q_error": float(max(q_errors)),
    }


def closure_criterion_audit() -> dict[str, float]:
    n, m = 7, 3
    H = make_spd(n, cond_target=30.0)
    U, _ = np.linalg.qr(RNG.standard_normal((n, n)))
    eigs = np.array([4.0, 2.0, 1.2, -0.5, -1.0, -2.5, 0.7], dtype=float)
    W = U @ np.diag(eigs) @ U.T
    B_adapted = U[:, :m]
    C_adapted = observer_from_B(H, B_adapted)
    Delta = sqrt_spd(H) @ W @ sqrt_spd(H)
    local_adapted = local_visible_calculus(H, C_adapted, Delta)

    B_random = random_stiefel(n, m)
    C_random = observer_from_B(H, B_random)
    local_random = local_visible_calculus(H, C_random, Delta)

    offplane = (np.eye(n) - B_adapted @ B_adapted.T) @ W @ B_adapted
    return {
        "adapted_Q_norm": float(np.linalg.norm(local_adapted.Q, ord="fro")),
        "adapted_offplane_norm": float(np.linalg.norm(offplane, ord="fro")),
        "random_Q_norm": float(np.linalg.norm(local_random.Q, ord="fro")),
        "random_eta": float(closure_scores(H, [Delta], B_random)["eta"]),
    }


def commuting_family_demo() -> dict[str, object]:
    n, m = 8, 3
    H = make_spd(n, cond_target=25.0)
    U, _ = np.linalg.qr(RNG.standard_normal((n, n)))

    lambdas = np.array(
        [
            [3.5, 1.0, 0.0, -0.5, 0.0, 2.0, -1.0, 0.0],
            [0.0, 2.8, 0.2, -1.0, 0.0, 1.2, 0.0, -0.3],
            [1.5, 0.0, 2.4, 0.0, -1.2, 0.0, 0.5, -0.8],
        ],
        dtype=float,
    )
    family = []
    for row in lambdas:
        W = U @ np.diag(row) @ U.T
        family.append(sym(sqrt_spd(H) @ W @ sqrt_spd(H)))

    mu = np.sum(lambdas**2, axis=0)
    top_idx = list(np.argsort(mu)[::-1][:m])
    B_star = U[:, top_idx]
    scores_star = closure_scores(H, family, B_star)

    exhaustive = []
    for subset in itertools.combinations(range(n), m):
        B = U[:, list(subset)]
        scores = closure_scores(H, family, B)
        exhaustive.append({"subset": list(subset), **scores})
    best_exact = min((row for row in exhaustive if row["H"] < 1e-10), key=lambda row: -row["S"])

    random_rows = []
    for _ in range(300):
        B = random_stiefel(n, m)
        random_rows.append(closure_scores(H, family, B))

    return {
        "top_indices_by_mu": top_idx,
        "mu_sorted": [float(x) for x in np.sort(mu)[::-1]],
        "adapted_scores": {k: float(v) for k, v in scores_star.items()},
        "best_exact_subset": best_exact,
        "random_eta_median": float(np.median([row["eta"] for row in random_rows])),
        "random_eta_min": float(np.min([row["eta"] for row in random_rows])),
        "random_S_median": float(np.median([row["S"] for row in random_rows])),
    }


def bridge_demo() -> dict[str, object]:
    n, m = 6, 2
    H0 = make_spd(n, cond_target=20.0)
    J = RNG.standard_normal((n, n))
    J = J - J.T
    bridge = dv_bridge(H0, J)
    W = whitened_delta(H0, bridge.delta_dv)
    evals, evecs = np.linalg.eigh(W)
    order = np.argsort(np.abs(evals))[::-1]
    B_adapted = evecs[:, order[:m]]
    C_adapted = observer_from_B(H0, B_adapted)
    local_adapted = local_visible_calculus(H0, C_adapted, bridge.delta_dv)

    random_metrics = []
    eps_demo = 2e-1
    for _ in range(200):
        B = random_stiefel(n, m)
        C = observer_from_B(H0, B)
        local = local_visible_calculus(H0, C, bridge.delta_dv)
        score = closure_scores(H0, [bridge.delta_dv], B)
        phi_eps = visible_precision(H0 + (eps_demo * eps_demo) * bridge.delta_dv, C)
        residual = phi_eps - visible_precision(H0, C) - (eps_demo * eps_demo) * local.V
        random_metrics.append(
            {
                "Q_norm": float(np.linalg.norm(local.Q, ord="fro")),
                "V_norm": float(np.linalg.norm(local.V, ord="fro")),
                "eta": score["eta"],
                "S": score["S"],
                "quartic_residual_scaled": float(np.linalg.norm(residual, ord="fro") / (eps_demo**4)),
            }
        )

    eps_grid = np.array([2e-1, 1e-1, 5e-2, 2e-2, 1e-2, 5e-3], dtype=float)
    residuals = []
    for eps in eps_grid:
        phi_eps = visible_precision(H0 + (eps * eps) * bridge.delta_dv, C_adapted)
        residual = phi_eps - np.eye(m) - (eps * eps) * local_adapted.V
        residuals.append(float(np.linalg.norm(residual, ord="fro")))

    slope, _ = np.polyfit(np.log(eps_grid), np.log(residuals), 1)

    return {
        "adapted_Q_norm": float(np.linalg.norm(local_adapted.Q, ord="fro")),
        "adapted_V_norm": float(np.linalg.norm(local_adapted.V, ord="fro")),
        "adapted_eta": float(closure_scores(H0, [bridge.delta_dv], B_adapted)["eta"]),
        "adapted_residual_scaled_at_eps_0p2": float(residuals[0] / (eps_grid[0] ** 4)),
        "random_eta_median": float(np.median([row["eta"] for row in random_metrics])),
        "random_eta_min": float(np.min([row["eta"] for row in random_metrics])),
        "random_Q_norm_median": float(np.median([row["Q_norm"] for row in random_metrics])),
        "random_V_norm_max": float(np.max([row["V_norm"] for row in random_metrics])),
        "random_quartic_residual_scaled_median_at_eps_0p2": float(
            np.median([row["quartic_residual_scaled"] for row in random_metrics])
        ),
        "eps_grid": [float(x) for x in eps_grid],
        "residuals_after_eps2V": residuals,
        "residual_loglog_slope": float(slope),
        "top_abs_eigenvalues": [float(abs(evals[i])) for i in order[:m]],
    }


def clean(obj: object) -> object:
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean(v) for v in obj]
    return obj


def main() -> None:
    results = {
        "paper_summary": {
            "main_implication": "observer synthesis reduces to subspace selection after H-whitening",
            "module_gap": "repo has the exact kernel but no closure-adapted observer construction layer",
        },
        "normal_form_audit": theorem_normal_form_audit(),
        "closure_criterion_audit": closure_criterion_audit(),
        "commuting_family_demo": commuting_family_demo(),
        "bridge_demo": bridge_demo(),
    }
    out_path = OUT / "closure_adapted_audit_results.json"
    cleaned = clean(results)
    out_path.write_text(json.dumps(cleaned, indent=2), encoding="utf-8")
    print(json.dumps(cleaned, indent=2))
    print(f"\nSaved results to {out_path}")


if __name__ == "__main__":
    main()
