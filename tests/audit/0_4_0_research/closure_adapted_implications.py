from __future__ import annotations

import itertools
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


def family_from_common_basis(H: np.ndarray, U: np.ndarray, lambdas: np.ndarray) -> list[np.ndarray]:
    H_half = sqrt_spd(H)
    family = []
    for row in lambdas:
        family.append(sym(H_half @ U @ np.diag(row) @ U.T @ H_half))
    return family


def sampled_random_summary(H: np.ndarray, family: list[np.ndarray], m: int, draws: int) -> dict[str, float]:
    rows = [closure_scores(H, family, random_stiefel(H.shape[0], m)) for _ in range(draws)]
    return {
        "eta_median": float(np.median([row["eta"] for row in rows])),
        "eta_min": float(np.min([row["eta"] for row in rows])),
        "S_median": float(np.median([row["S"] for row in rows])),
        "S_max": float(np.max([row["S"] for row in rows])),
    }


def commuting_solver_stress(trials: int = 40, random_draws: int = 400) -> dict[str, object]:
    exact_matches = 0
    top_mu_beats_random_max = 0
    adapted_eta = []
    random_eta_medians = []
    signal_gain_vs_random_median = []
    signal_gain_vs_random_max = []
    representative = None

    for trial_idx in range(trials):
        n = int(RNG.integers(7, 10))
        m = int(RNG.integers(2, 4))
        task_count = int(RNG.integers(3, 6))
        H = make_spd(n, cond_target=float(10.0 ** RNG.uniform(0.4, 1.7)))
        U, _ = np.linalg.qr(RNG.standard_normal((n, n)))
        lambdas = RNG.normal(scale=1.5, size=(task_count, n))
        # Force a real low-rank shared task family so exact closure exists and is nontrivial.
        lambdas[:, m + 2 :] = 0.0
        family = family_from_common_basis(H, U, lambdas)

        mu = np.sum(lambdas**2, axis=0)
        top_idx = tuple(int(i) for i in np.argsort(mu)[::-1][:m])
        top_set = tuple(sorted(top_idx))
        B_star = U[:, list(top_idx)]
        scores_star = closure_scores(H, family, B_star)

        exhaustive = []
        for subset in itertools.combinations(range(n), m):
            B = U[:, list(subset)]
            scores = closure_scores(H, family, B)
            exhaustive.append({"subset": tuple(int(i) for i in subset), **scores})
        best_exact = min((row for row in exhaustive if abs(row["H"]) < 1e-10), key=lambda row: -row["S"])
        if tuple(sorted(best_exact["subset"])) == top_set:
            exact_matches += 1

        random_summary = sampled_random_summary(H, family, m, random_draws)
        if scores_star["S"] > random_summary["S_max"]:
            top_mu_beats_random_max += 1
        adapted_eta.append(float(abs(scores_star["eta"])))
        random_eta_medians.append(float(random_summary["eta_median"]))
        signal_gain_vs_random_median.append(float(scores_star["S"] / max(random_summary["S_median"], 1e-15)))
        signal_gain_vs_random_max.append(float(scores_star["S"] / max(random_summary["S_max"], 1e-15)))

        if representative is None or signal_gain_vs_random_median[-1] > representative["gain_vs_random_median"]:
            representative = {
                "n": n,
                "m": m,
                "task_count": task_count,
                "top_idx": list(top_idx),
                "mu_top": [float(x) for x in np.sort(mu)[::-1][: min(6, n)]],
                "adapted_scores": {k: float(v) for k, v in scores_star.items()},
                "best_exact_subset": list(best_exact["subset"]),
                "random_summary": random_summary,
                "gain_vs_random_median": signal_gain_vs_random_median[-1],
            }

    return {
        "trials": trials,
        "exact_match_rate": float(exact_matches / trials),
        "adapted_eta_max": float(max(adapted_eta)),
        "random_eta_median_over_trials": float(np.median(random_eta_medians)),
        "signal_gain_vs_random_median_median": float(np.median(signal_gain_vs_random_median)),
        "signal_gain_vs_random_max_median": float(np.median(signal_gain_vs_random_max)),
        "beats_sampled_random_max_rate": float(top_mu_beats_random_max / trials),
        "representative_case": representative,
    }


def bridge_stress(trials: int = 50, random_draws: int = 400) -> dict[str, object]:
    adapted_q = []
    adapted_eta = []
    adapted_scaled_residual = []
    random_eta_medians = []
    random_s_medians = []
    gain_vs_random_median = []
    gain_vs_random_max = []
    adapted_beats_random_max = 0
    representative = None

    for _ in range(trials):
        n = int(RNG.integers(6, 10))
        m = 2
        H0 = make_spd(n, cond_target=float(10.0 ** RNG.uniform(0.3, 1.5)))
        J = RNG.standard_normal((n, n))
        J = J - J.T
        bridge = dv_bridge(H0, J)
        W = whitened_delta(H0, bridge.delta_dv)
        evals, evecs = np.linalg.eigh(W)
        order = np.argsort(np.abs(evals))[::-1]
        B_star = evecs[:, order[:m]]
        C_star = observer_from_B(H0, B_star)
        local_star = local_visible_calculus(H0, C_star, bridge.delta_dv)
        scores_star = closure_scores(H0, [bridge.delta_dv], B_star)

        eps = 0.2
        phi_eps = visible_precision(H0 + (eps * eps) * bridge.delta_dv, C_star)
        residual = phi_eps - np.eye(m) - (eps * eps) * local_star.V
        scaled_residual = float(np.linalg.norm(residual, ord="fro") / (eps**4))

        random_summary = sampled_random_summary(H0, [bridge.delta_dv], m, random_draws)
        if scores_star["S"] > random_summary["S_max"]:
            adapted_beats_random_max += 1

        adapted_q.append(float(np.linalg.norm(local_star.Q, ord="fro")))
        adapted_eta.append(float(abs(scores_star["eta"])))
        adapted_scaled_residual.append(scaled_residual)
        random_eta_medians.append(float(random_summary["eta_median"]))
        random_s_medians.append(float(random_summary["S_median"]))
        gain_vs_random_median.append(float(scores_star["S"] / max(random_summary["S_median"], 1e-15)))
        gain_vs_random_max.append(float(scores_star["S"] / max(random_summary["S_max"], 1e-15)))

        if representative is None or gain_vs_random_median[-1] > representative["gain_vs_random_median"]:
            representative = {
                "ambient_dim": n,
                "top_abs_eigenvalues": [float(abs(evals[i])) for i in order[:m]],
                "adapted_scores": {k: float(v) for k, v in scores_star.items()},
                "adapted_Q_norm": adapted_q[-1],
                "adapted_scaled_residual": scaled_residual,
                "random_summary": random_summary,
                "gain_vs_random_median": gain_vs_random_median[-1],
            }

    return {
        "trials": trials,
        "adapted_Q_max": float(max(adapted_q)),
        "adapted_eta_max": float(max(adapted_eta)),
        "adapted_scaled_residual_median": float(np.median(adapted_scaled_residual)),
        "random_eta_median_over_trials": float(np.median(random_eta_medians)),
        "signal_gain_vs_random_median_median": float(np.median(gain_vs_random_median)),
        "signal_gain_vs_random_max_median": float(np.median(gain_vs_random_max)),
        "beats_sampled_random_max_rate": float(adapted_beats_random_max / trials),
        "representative_case": representative,
    }


def approximate_compatibility_sweep(
    thetas: tuple[float, ...] = (0.0, 0.08, 0.16, 0.28, 0.42, 0.65),
    trials_per_theta: int = 12,
    search_draws: int = 2500,
) -> dict[str, object]:
    rows = []
    for theta in thetas:
        designed_etas = []
        designed_signals = []
        median_etas = []
        commutator_norms = []
        for _ in range(trials_per_theta):
            n, m = 10, 3
            H = make_spd(n, cond_target=float(10.0 ** RNG.uniform(0.4, 1.5)))
            U, _ = np.linalg.qr(RNG.standard_normal((n, n)))

            W1 = U @ np.diag([3.4, 2.4, 1.9, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) @ U.T
            rot = np.eye(n)
            rot[:4, :4] = np.array(
                [
                    [np.cos(theta), -np.sin(theta), 0.0, 0.0],
                    [np.sin(theta), np.cos(theta), 0.0, 0.0],
                    [0.0, 0.0, np.cos(theta), -np.sin(theta)],
                    [0.0, 0.0, np.sin(theta), np.cos(theta)],
                ]
            )
            W2 = rot @ U @ np.diag([3.1, 2.1, 1.8, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) @ U.T @ rot.T
            family = [sym(sqrt_spd(H) @ W1 @ sqrt_spd(H)), sym(sqrt_spd(H) @ W2 @ sqrt_spd(H))]

            aggregate = sym(W1 @ W1 + W2 @ W2)
            evals, evecs = np.linalg.eigh(aggregate)
            order = np.argsort(evals)[::-1]
            B_design = evecs[:, order[:m]]
            designed = closure_scores(H, family, B_design)

            scores = [closure_scores(H, family, random_stiefel(n, m)) for _ in range(search_draws)]
            designed_etas.append(float(designed["eta"]))
            designed_signals.append(float(designed["S"]))
            median_etas.append(float(np.median([row["eta"] for row in scores])))
            commutator_norms.append(float(np.linalg.norm(W1 @ W2 - W2 @ W1, ord="fro")))

        rows.append(
            {
                "theta": float(theta),
                "commutator_norm_median": float(np.median(commutator_norms)),
                "designed_eta_median": float(np.median(designed_etas)),
                "designed_eta_max": float(np.max(designed_etas)),
                "designed_signal_median": float(np.median(designed_signals)),
                "random_eta_median": float(np.median(median_etas)),
            }
        )

    return {
        "ambient_dim": 10,
        "observer_dim": 3,
        "task_count": 2,
        "rows": rows,
        "interpretation": (
            "For the aggregate-mode design heuristic, eta is exactly zero in the commuting limit and rises as the commutator grows. "
            "This makes eta a practical compatibility diagnostic even when exact closure is unavailable."
        ),
    }


def main() -> None:
    results = {
        "commuting_solver_stress": commuting_solver_stress(),
        "bridge_stress": bridge_stress(),
        "approximate_compatibility_sweep": approximate_compatibility_sweep(),
    }
    out_path = OUT / "closure_adapted_implications.json"
    cleaned = clean(results)
    out_path.write_text(json.dumps(cleaned, indent=2), encoding="utf-8")
    print(json.dumps(cleaned, indent=2))
    print(f"\nSaved results to {out_path}")


if __name__ == "__main__":
    main()
