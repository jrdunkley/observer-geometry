"""Research checks for previously shelved 0.3.2 feature candidates.

This is audit-only. It tests whether several older withheld ideas now have
safe semantics under the layered doctrine:

- rank-k covariance/Fisher perturbations;
- declared-ladder dimension-cost comparisons;
- guarded fibre-dominance diagnostics;
- affine-hidden stage sign and branch reversal;
- ensemble finite-candidate residual margins.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from nomogeo import residual_margin_ordering, staged_affine_hidden_elimination, variable_precision_affine_hidden_reduction


OUT = Path("audit/outputs/0_3_2_shelved_feature_recovery_check.json")


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def rank_k_covariance_update(seed: int = 39601) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    rows = []
    for n, m, k in [(8, 3, 1), (9, 4, 2), (10, 4, 3)]:
        raw = rng.normal(size=(n, n))
        sigma0 = raw.T @ raw / n + 0.9 * np.eye(n)
        f = rng.normal(size=(n, k))
        sigma1 = sigma0 + f @ f.T
        h0 = np.linalg.inv(sigma0)
        h1 = np.linalg.inv(sigma1)
        phi0 = np.linalg.inv(sigma0[:m, :m])
        phi1 = np.linalg.inv(sigma1[:m, :m])
        direct = sym((h1[:m, :m] - phi1) - (h0[:m, :m] - phi0))

        fv = f[:m, :]
        u = h0[:m, :] @ f
        w = phi0 @ fv
        beta = np.linalg.inv(np.eye(k) + f.T @ h0 @ f)
        gamma = np.linalg.inv(np.eye(k) + fv.T @ phi0 @ fv)
        formula = sym(w @ gamma @ w.T - u @ beta @ u.T)
        singular = np.linalg.svd(formula, compute_uv=False)
        tol = 1.0e-9 * max(1.0, float(np.max(singular)))
        rows.append(
            {
                "n": n,
                "visible_dim": m,
                "k": k,
                "formula_residual": float(np.max(np.abs(direct - formula))),
                "rank_formula": int(np.sum(singular > tol)),
                "rank_bound_2k": 2 * k,
                "singular_values": [float(x) for x in singular],
            }
        )
    return {"rows": rows}


def declared_ladder_dimension_cost() -> dict[str, object]:
    candidates = [
        {"name": "small_clean", "score": 4.0, "dim": 1},
        {"name": "medium", "score": 5.8, "dim": 2},
        {"name": "large", "score": 6.4, "dim": 4},
    ]
    breakpoints: list[float] = []
    for i, left in enumerate(candidates):
        for right in candidates[i + 1 :]:
            denom = left["dim"] - right["dim"]
            if denom != 0:
                c = (left["score"] - right["score"]) / denom
                if c >= 0.0:
                    breakpoints.append(float(c))
    grid = sorted({0.0, *breakpoints, 10.0})
    probes = sorted({0.0, 10.0, *[max(0.0, c - 1.0e-6) for c in breakpoints], *[c + 1.0e-6 for c in breakpoints]})
    rows = []
    for cost in probes:
        adjusted = [(item["score"] - cost * item["dim"], item["name"]) for item in candidates]
        best_score = max(score for score, _name in adjusted)
        winners = sorted(name for score, name in adjusted if abs(score - best_score) <= 1.0e-10)
        rows.append({"cost": float(cost), "winners": winners, "adjusted_scores": {name: float(score) for score, name in adjusted}})
    return {"candidates": candidates, "breakpoints": grid, "probe_rows": rows}


def centered_norm(values: np.ndarray) -> float:
    centered = values - float(np.mean(values))
    return float(np.linalg.norm(centered))


def guarded_fibre_dominance() -> dict[str, object]:
    v = np.linspace(-1.0, 1.0, 101)
    floor = 1.0e-6
    rows = []
    for base_action in ["flat", "quartic"]:
      for coupling_amp in [0.0, 1.0e-4, 1.0e-2, 1.0]:
        action = np.zeros_like(v) if base_action == "flat" else 0.08 * v**4
        coupling = coupling_amp * (0.2 + v**2)
        precision = np.exp(0.7 * v)
        variational = action - 0.5 * coupling**2 / precision
        fibre = 0.5 * np.log(precision)
        denom = centered_norm(variational)
        numer = centered_norm(fibre)
        ratio = None if denom < floor else numer / denom
        rows.append(
            {
                "base_action": base_action,
                "coupling_amp": coupling_amp,
                "variational_centered_norm": denom,
                "fibre_centered_norm": numer,
                "ratio": None if ratio is None else float(ratio),
                "ratio_defined": ratio is not None,
            }
        )
    return {"denominator_floor": floor, "rows": rows}


def affine_hidden_branch_and_stage_signs() -> dict[str, object]:
    # Branch flip: same variational action, different precision volume.
    branch = variable_precision_affine_hidden_reduction(
        np.array([0.0, 0.0]),
        np.zeros((2, 1)),
        np.array([[[0.1]], [[10.0]]]),
    )
    # Stage signs: same D, coupling controls whether log-det or saddle term wins.
    d = np.diag([4.0, 2.0])
    rows = []
    for j0 in [0.0, 1.0, 3.0]:
        stage = staged_affine_hidden_elimination(0.0, np.array([j0, 0.0]), d, [0])
        rows.append(
            {
                "j0": j0,
                "action_shift": float(stage.action_shift),
                "shift_sign": int(np.sign(stage.action_shift)),
                "logdet_half": float(0.5 * np.log(4.0)),
                "saddle_half": float(0.5 * j0 * j0 / 4.0),
            }
        )
    return {
        "branch_visible_action": [float(x) for x in np.ravel(branch.visible_action)],
        "branch_variational_action": [float(x) for x in np.ravel(branch.variational_action)],
        "branch_fibre_volume": [float(x) for x in np.ravel(branch.fibre_volume)],
        "stage_rows": rows,
    }


def ensemble_residual_margin(seed: int = 39602) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    scenarios = []
    for label, n, centers in [
        ("plausible_but_not_certified", 80, {"A": 0.72, "B": 0.66, "C": 0.60}),
        ("large_gap_certified", 800, {"A": 0.78, "B": 0.50, "C": 0.42}),
    ]:
        samples = {name: np.clip(rng.normal(center, 0.04, size=n), 0.0, 1.0) for name, center in centers.items()}
        delta = 0.05
        m = len(samples)
        hoeffding = math.sqrt(math.log(2.0 * m / delta) / (2.0 * n))
        means = {name: float(np.mean(values)) for name, values in samples.items()}
        winner = max(means, key=means.get)
        comparisons = []
        robust_all = True
        for name in samples:
            if name == winner:
                continue
            gap = means[winner] - means[name]
            cert = residual_margin_ordering(gap, hoeffding)
            robust_all = robust_all and bool(cert.robust)
            comparisons.append(
                {
                    "winner": winner,
                    "other": name,
                    "empirical_gap": float(gap),
                    "per_candidate_bound": float(hoeffding),
                    "required_gap": float(cert.required_gap),
                    "robust": bool(cert.robust),
                }
            )
        scenarios.append(
            {
                "label": label,
                "sample_count": n,
                "means": means,
                "winner": winner,
                "hoeffding_bound": float(hoeffding),
                "robust_all": robust_all,
                "comparisons": comparisons,
            }
        )
    return {"scenarios": scenarios}


def main() -> None:
    results = {
        "rank_k_covariance_update": rank_k_covariance_update(),
        "declared_ladder_dimension_cost": declared_ladder_dimension_cost(),
        "guarded_fibre_dominance": guarded_fibre_dominance(),
        "affine_hidden_branch_and_stage_signs": affine_hidden_branch_and_stage_signs(),
        "ensemble_residual_margin": ensemble_residual_margin(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
