from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import dv_bridge, local_visible_calculus, visible_precision

from papers.closure_adapted_audit import (
    RNG,
    clean,
    closure_scores,
    make_spd,
    observer_from_B,
    random_stiefel,
    sqrt_spd,
    sym,
    whitened_delta,
)


OUT = Path(__file__).resolve().parent / "closure_adapted_outputs"
OUT.mkdir(exist_ok=True)


def sample_random_scores(H: np.ndarray, family: list[np.ndarray], m: int, draws: int = 500) -> dict[str, float]:
    rows = []
    n = H.shape[0]
    for _ in range(draws):
        B = random_stiefel(n, m)
        rows.append(closure_scores(H, family, B))
    return {
        "eta_median": float(np.median([row["eta"] for row in rows])),
        "eta_min": float(np.min([row["eta"] for row in rows])),
        "S_median": float(np.median([row["S"] for row in rows])),
        "S_max": float(np.max([row["S"] for row in rows])),
    }


def clinical_panel_demo() -> dict[str, object]:
    n, m = 12, 3
    H = make_spd(n, cond_target=40.0)
    U, _ = np.linalg.qr(RNG.standard_normal((n, n)))

    # Three disease signatures and one treatment response, all living on the same
    # latent pathway basis. This is the commuting sweet spot of the paper.
    lambdas = np.array(
        [
            [3.2, 1.9, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.7, 0.0, 1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 2.1, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [1.9, 1.1, 2.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    family = [sym(sqrt_spd(H) @ U @ np.diag(row) @ U.T @ sqrt_spd(H)) for row in lambdas]
    mu = np.sum(lambdas**2, axis=0)
    top_idx = list(np.argsort(mu)[::-1][:m])
    B_star = U[:, top_idx]
    scores_star = closure_scores(H, family, B_star)
    random_summary = sample_random_scores(H, family, m, draws=800)

    local_checks = []
    C_star = observer_from_B(H, B_star)
    for Delta in family:
        local = local_visible_calculus(H, C_star, Delta)
        local_checks.append(
            {
                "Q_norm": float(np.linalg.norm(local.Q, ord="fro")),
                "V_norm": float(np.linalg.norm(local.V, ord="fro")),
            }
        )

    return {
        "domain": "biomarker panel design",
        "ambient_dim": n,
        "panel_dim": m,
        "task_count": len(family),
        "top_indices": top_idx,
        "mu_sorted": [float(x) for x in np.sort(mu)[::-1][:6]],
        "adapted_scores": {k: float(v) for k, v in scores_star.items()},
        "adapted_local_checks": local_checks,
        "random_summary": random_summary,
        "claim": (
            "A three-channel panel can exactly close a four-task clinical family living in a 12D latent space. "
            "Random three-channel observers leak heavily."
        ),
    }


def bridge_monitor_demo() -> dict[str, object]:
    n, m = 8, 2
    H0 = make_spd(n, cond_target=25.0)
    J = RNG.standard_normal((n, n))
    J = J - J.T
    bridge = dv_bridge(H0, J)
    W = whitened_delta(H0, bridge.delta_dv)
    evals, evecs = np.linalg.eigh(W)
    order = np.argsort(np.abs(evals))[::-1]
    B_star = evecs[:, order[:m]]
    C_star = observer_from_B(H0, B_star)
    local_star = local_visible_calculus(H0, C_star, bridge.delta_dv)
    score_star = closure_scores(H0, [bridge.delta_dv], B_star)
    random_summary = sample_random_scores(H0, [bridge.delta_dv], m, draws=800)

    eps = 0.2
    phi_star_eps = visible_precision(H0 + (eps * eps) * bridge.delta_dv, C_star)
    quartic_residual_star = float(
        np.linalg.norm(phi_star_eps - np.eye(m) - (eps * eps) * local_star.V, ord="fro") / (eps**4)
    )

    return {
        "domain": "nonequilibrium monitor design",
        "ambient_dim": n,
        "observer_dim": m,
        "top_abs_eigenvalues": [float(abs(evals[i])) for i in order[:m]],
        "adapted_scores": {k: float(v) for k, v in score_star.items()},
        "adapted_Q_norm": float(np.linalg.norm(local_star.Q, ord="fro")),
        "adapted_V_norm": float(np.linalg.norm(local_star.V, ord="fro")),
        "adapted_quartic_residual_scaled": quartic_residual_star,
        "random_summary": random_summary,
        "claim": (
            "A two-channel bridge-adapted observer exposes the visible DV curvature while driving quartic hidden birth "
            "essentially to zero."
        ),
    }


def approximate_noncommuting_demo() -> dict[str, object]:
    n, m = 10, 3
    H = make_spd(n, cond_target=30.0)
    U, _ = np.linalg.qr(RNG.standard_normal((n, n)))

    # Two tasks share most of their spectral structure but not all of it.
    W1 = U @ np.diag([3.5, 2.2, 1.8, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) @ U.T
    mix = np.eye(n)
    theta = 0.35
    mix[:4, :4] = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0, 0.0],
            [np.sin(theta), np.cos(theta), 0.0, 0.0],
            [0.0, 0.0, np.cos(theta), -np.sin(theta)],
            [0.0, 0.0, np.sin(theta), np.cos(theta)],
        ]
    )
    W2 = mix @ U @ np.diag([3.2, 2.0, 1.6, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) @ U.T @ mix.T
    family = [sym(sqrt_spd(H) @ W1 @ sqrt_spd(H)), sym(sqrt_spd(H) @ W2 @ sqrt_spd(H))]

    rows = []
    for _ in range(4000):
        B = random_stiefel(n, m)
        rows.append({"B": B, **closure_scores(H, family, B)})
    rows.sort(key=lambda row: (row["eta"], -row["S"]))
    best = rows[0]
    median_eta = float(np.median([row["eta"] for row in rows]))
    median_S = float(np.median([row["S"] for row in rows]))

    return {
        "domain": "approximate task-matched sensing",
        "ambient_dim": n,
        "observer_dim": m,
        "best_eta_found": float(best["eta"]),
        "best_S_found": float(best["S"]),
        "random_eta_median": median_eta,
        "random_S_median": median_S,
        "eta_improvement_factor": float(median_eta / max(best["eta"], 1e-15)),
        "S_improvement_factor": float(best["S"] / max(median_S, 1e-15)),
        "claim": (
            "Even when exact closure is impossible, the leakage objective still finds observers with dramatically "
            "lower hidden birth than generic choices."
        ),
    }


def main() -> None:
    results = {
        "capability_shift": {
            "before": "audit a chosen observer",
            "after": "synthesise an observer from a task family when the paper's conditions apply",
        },
        "use_cases": {
            "clinical_panel_demo": clinical_panel_demo(),
            "bridge_monitor_demo": bridge_monitor_demo(),
            "approximate_noncommuting_demo": approximate_noncommuting_demo(),
        },
    }
    out_path = OUT / "closure_adapted_use_cases.json"
    cleaned = clean(results)
    out_path.write_text(json.dumps(cleaned, indent=2), encoding="utf-8")
    print(json.dumps(cleaned, indent=2))
    print(f"\nSaved results to {out_path}")


if __name__ == "__main__":
    main()
