"""Quarantined 0.3.2 research probes for external Claude notes.

These probes do not modify the module. They test claims that are tempting to
promote from the notes but need exact scope discipline before becoming 0.3.2
planning items.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import (
    exact_branch_hessian,
    gaussian_data_processing_contraction,
    rank_one_covariance_perturbation,
    variable_precision_affine_hidden_reduction,
    weighted_family_frontier_scores,
)
from nomogeo.exceptions import InputValidationError


OUT = Path("audit/outputs/0_3_2_research_probe.json")


def _sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def _rank(a: np.ndarray, tol: float = 1.0e-9) -> int:
    return int(np.sum(np.linalg.svd(a, compute_uv=False) > tol))


def _make_correlated_covariance(m: int, h: int, rho: float = 0.18) -> np.ndarray:
    rng = np.random.default_rng(32001)
    cross = rng.normal(size=(m, h))
    cross /= np.linalg.norm(cross, ord=2)
    sigma = np.block(
        [
            [np.eye(m), rho * cross],
            [rho * cross.T, 1.4 * np.eye(h)],
        ]
    )
    return _sym(sigma)


def rank_one_boundary_probe() -> dict:
    """Separate the actual rank criterion from a hidden-only folk theorem."""
    rng = np.random.default_rng(32002)
    n, m, eps = 8, 3, 0.37
    hidden_dim = n - m

    f_hidden = np.zeros(n)
    f_hidden[m:] = rng.normal(size=hidden_dim)

    block_cov = np.block(
        [
            [np.diag([1.1, 1.3, 1.7]), np.zeros((m, hidden_dim))],
            [np.zeros((hidden_dim, m)), np.diag(np.linspace(0.9, 1.6, hidden_dim))],
        ]
    )
    hidden_block = rank_one_covariance_perturbation(block_cov, f_hidden, m, eps)

    correlated_cov = _make_correlated_covariance(m, hidden_dim)
    hidden_correlated = rank_one_covariance_perturbation(correlated_cov, f_hidden, m, eps)

    f_mixed = rng.normal(size=n)
    # Block-diagonal covariance makes u_V and w_V collinear whenever only one
    # visible coordinate is active, even though the signal has visible mass.
    f_mixed[:m] = 0.0
    f_mixed[1] = 0.8
    mixed_aligned = rank_one_covariance_perturbation(block_cov, f_mixed, m, eps)

    random_cov = None
    random_two = None
    for _ in range(500):
        a = rng.normal(size=(n, n))
        sigma = _sym(a.T @ a / n + 0.8 * np.eye(n))
        f = rng.normal(size=n)
        res = rank_one_covariance_perturbation(sigma, f, m, eps)
        if res.update_rank == 2:
            random_cov = sigma
            random_two = res
            break
    if random_cov is None or random_two is None:
        raise RuntimeError("failed to find a rank-two coloured example")

    def pack(res) -> dict:
        return {
            "rank": int(res.update_rank),
            "one_channel": bool(res.one_channel),
            "formula_residual": float(res.formula_residual),
            "alignment": float(res.direction_alignment),
            "signal_visible_norm": float(np.linalg.norm(res.covariance_perturbed[:m, :m] - res.covariance_base[:m, :m])),
            "singular_values": [float(x) for x in res.singular_values],
        }

    return {
        "hidden_only_block_diagonal": pack(hidden_block),
        "hidden_only_correlated_background": pack(hidden_correlated),
        "visible_signal_block_aligned": pack(mixed_aligned),
        "generic_coloured_background": pack(random_two),
        "conclusion": (
            "one-channel is controlled by collinearity or vanishing of u_V=(H0 f)_V "
            "and w_V=Phi0 f_V, not simply by zero visible projection of f"
        ),
    }


def near_branch_failure_probe() -> dict:
    """Show why a non-throwing near-branch diagnostic would be useful."""
    B = np.eye(4, 2)
    base = np.diag([1.0, 1.6, 2.4, 3.2])
    family = []
    for eps in [0.0, 1.0e-10, 1.0e-7, 1.0e-5, 1.0e-3]:
        op = base.copy()
        op[0, 2] = eps
        op[2, 0] = eps
        op[1, 3] = -0.5 * eps
        op[3, 1] = -0.5 * eps
        off_block = float(np.linalg.norm(op[2:, :2], ord="fro"))
        try:
            result = exact_branch_hessian([op], B, mu=0.4)
            status = result.status
            returned_off = float(result.off_block_norm)
        except InputValidationError as exc:
            status = f"rejected: {exc}"
            returned_off = None
        family.append({"epsilon": eps, "manual_off_block_norm": off_block, "status": status, "returned_off_block_norm": returned_off})
    return {
        "grid": family,
        "conclusion": "current exact branch Hessian correctly rejects non-exact branches, but discards a useful continuous near-branch scale",
    }


def fibre_dominance_normalisation_probe() -> dict:
    """Demonstrate that a fibre/variational ratio needs a declared norm."""
    x = np.linspace(-1.0, 1.0, 51)
    D = np.zeros((x.size, 2, 2))
    D[:, 0, 0] = 1.0 + 0.6 * x**2
    D[:, 1, 1] = 1.2 + 0.2 * np.cos(2.0 * x)
    A = np.zeros_like(x)

    ratios = []
    for amp in [0.0, 1.0e-4, 1.0e-2, 1.0]:
        J = np.zeros((x.size, 2))
        J[:, 0] = amp * x
        J[:, 1] = 0.3 * amp * x**2
        result = variable_precision_affine_hidden_reduction(A, J, D)
        var_centred = result.variational_action - float(np.mean(result.variational_action))
        fib_centred = result.fibre_volume - float(np.mean(result.fibre_volume))
        denom = float(np.linalg.norm(var_centred))
        numer = float(np.linalg.norm(fib_centred))
        ratios.append(
            {
                "coupling_amplitude": amp,
                "centred_l2_variational_norm": denom,
                "centred_l2_fibre_norm": numer,
                "ratio": None if denom == 0.0 else numer / denom,
            }
        )
    return {
        "ratios": ratios,
        "conclusion": "fibre dominance is meaningful only after choosing a centred variation norm and guarding small denominators",
    }


def divergence_contraction_probe() -> dict:
    """Compare Gaussian KL and Hellinger contraction ratios over perturbation size."""
    rng = np.random.default_rng(32003)
    n, m, k = 7, 5, 2
    C1 = np.eye(m, n)
    D = rng.normal(size=(k, m))
    q, _ = np.linalg.qr(D.T)
    D = q[:, :k].T
    a = rng.normal(size=(n, n))
    cov_base = _sym(a.T @ a / n + np.eye(n))
    direction = _sym(rng.normal(size=(n, n)))
    rows = []
    for scale in [1.0e-3, 1.0e-2, 1.0e-1, 0.5, 1.0]:
        cov_pert = _sym(cov_base + scale * direction)
        eig = np.linalg.eigvalsh(cov_pert)
        if float(np.min(eig)) <= 0.05:
            cov_pert += (0.06 - float(np.min(eig))) * np.eye(n)
        H1 = np.linalg.inv(cov_base)
        H2 = np.linalg.inv(cov_pert)
        res = gaussian_data_processing_contraction(H1, H2, C1, D)
        kl_ratio = res.forward_kl_coarse / res.forward_kl_fine if res.forward_kl_fine > 0.0 else None
        hell_ratio = res.hellinger_sq_coarse / res.hellinger_sq_fine if res.hellinger_sq_fine > 0.0 else None
        rows.append(
            {
                "scale": scale,
                "forward_kl_fine": float(res.forward_kl_fine),
                "forward_kl_coarse": float(res.forward_kl_coarse),
                "kl_contraction_ratio": None if kl_ratio is None else float(kl_ratio),
                "hellinger_contraction_ratio": None if hell_ratio is None else float(hell_ratio),
                "ratio_gap": None if kl_ratio is None or hell_ratio is None else float(abs(kl_ratio - hell_ratio)),
            }
        )
    return {
        "rows": rows,
        "conclusion": "KL and Hellinger contraction ratios agree only as a small-perturbation heuristic, not as an invariant fragility score",
    }


def declared_ladder_probe() -> dict:
    """Show frontier selection is selection over a declared observer ladder."""
    rng = np.random.default_rng(32004)
    n = 8
    family = []
    # A rank-two signal concentrated in coordinates 5 and 6.
    for strength in [0.8, 1.1, 1.5]:
        v = np.zeros(n)
        v[5] = strength
        v[6] = 0.7 * strength
        noise = _sym(0.04 * rng.normal(size=(n, n)))
        family.append(_sym(np.outer(v, v) + noise))

    rows_coordinate = []
    rows_rotated = []
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    for m in range(1, n):
        Bc = np.eye(n, m)
        Br = q[:, :m]
        rc = weighted_family_frontier_scores(family, Bc, mu=0.25)
        rr = weighted_family_frontier_scores(family, Br, mu=0.25)
        rows_coordinate.append({"m": m, "penalized": float(rc.penalized_score), "leakage": float(rc.leakage), "visible": float(rc.visible_score)})
        rows_rotated.append({"m": m, "penalized": float(rr.penalized_score), "leakage": float(rr.leakage), "visible": float(rr.visible_score)})
    coord_best = max(rows_coordinate, key=lambda row: row["penalized"])
    rot_best = max(rows_rotated, key=lambda row: row["penalized"])
    return {
        "coordinate_ladder": rows_coordinate,
        "rotated_ladder": rows_rotated,
        "coordinate_best": coord_best,
        "rotated_best": rot_best,
        "conclusion": "frontier scores support declared-ladder comparison; global observer selection needs an optimizer and a search domain",
    }


def main() -> None:
    results = {
        "rank_one_boundary": rank_one_boundary_probe(),
        "near_branch_failure": near_branch_failure_probe(),
        "fibre_dominance_normalisation": fibre_dominance_normalisation_probe(),
        "divergence_contraction": divergence_contraction_probe(),
        "declared_ladder_frontier": declared_ladder_probe(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
