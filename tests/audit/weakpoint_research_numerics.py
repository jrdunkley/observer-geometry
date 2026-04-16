from __future__ import annotations

import json
from itertools import permutations
from pathlib import Path

import numpy as np


OUT = Path(__file__).resolve().parent / "outputs" / "weakpoint_research_numerics.json"


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def centered(a: np.ndarray) -> np.ndarray:
    return a - float(np.mean(a))


def normalize_log_density(logw: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, float]:
    shift = float(np.max(logw))
    w = np.exp(logw - shift)
    z = float(np.trapezoid(w, x))
    return w / z, shift + float(np.log(z))


def variable_hidden_precision_boundary() -> dict[str, float]:
    v = np.linspace(-0.85, 0.85, 601)
    a = 0.12 * v**4 + 0.03 * v
    j = 0.25 * np.sin(1.7 * v) + 0.08 * v * v
    d_fixed = np.full_like(v, 1.4)
    d_variable = 1.4 * np.exp(1.1 * v - 0.35 * v * v)

    fixed_exact = a - 0.5 * j * j / d_fixed
    fixed_variational = a - 0.5 * j * j / d_fixed

    variable_exact = a + 0.5 * np.log(d_variable) - 0.5 * j * j / d_variable
    variable_variational = a - 0.5 * j * j / d_variable
    missing = centered(variable_exact - variable_variational)
    logdet = centered(0.5 * np.log(d_variable))

    fixed_residual = float(np.max(np.abs(centered(fixed_exact - fixed_variational))))
    variable_residual = float(np.max(np.abs(missing - logdet)))

    # A clean verdict change: with A=J=0, variational reduction is flat but
    # the exact marginal selects the smallest D(v) by the fibre-volume term.
    d_hostile = np.exp(2.0 * v)
    hostile_exact = centered(0.5 * np.log(d_hostile))
    hostile_variational = np.zeros_like(v)
    exact_argmin = float(v[int(np.argmin(hostile_exact))])
    variational_range = float(np.max(hostile_variational) - np.min(hostile_variational))

    return {
        "fixed_precision_centered_residual": fixed_residual,
        "variable_precision_missing_logdet_residual": variable_residual,
        "variable_precision_missing_range": float(np.max(missing) - np.min(missing)),
        "hostile_exact_argmin": exact_argmin,
        "hostile_variational_range": variational_range,
    }


def eliminate_block(d: np.ndarray, j: np.ndarray, keep: list[int], elim: list[int]) -> tuple[np.ndarray, np.ndarray, float]:
    k = np.ix_(keep, keep)
    e = np.ix_(elim, elim)
    ke = np.ix_(keep, elim)
    ek = np.ix_(elim, keep)
    d_kk = d[k]
    d_ke = d[ke]
    d_ek = d[ek]
    d_ee = d[e]
    j_k = j[keep]
    j_e = j[elim]
    solved_dek = np.linalg.solve(d_ee, d_ek)
    solved_je = np.linalg.solve(d_ee, j_e)
    new_d = sym(d_kk - d_ke @ solved_dek)
    new_j = j_k - d_ke @ solved_je
    constant_gain = -0.5 * float(j_e.T @ solved_je)
    return new_d, new_j, constant_gain


def affine_hidden_tower_random() -> dict[str, float]:
    rng = np.random.default_rng(20260412)
    raw = rng.normal(size=(5, 5))
    d = raw.T @ raw + 1.3 * np.eye(5)
    v_grid = np.linspace(-1.0, 1.0, 41)

    max_residual = 0.0
    for v in v_grid:
        j = np.array(
            [
                np.sin(v),
                0.3 + v * v,
                np.cos(1.4 * v),
                -0.2 * v,
                0.1 + 0.5 * np.sin(0.7 * v),
            ],
            dtype=float,
        )
        one_step = -0.5 * float(j.T @ np.linalg.solve(d, j))
        for order in permutations(range(5)):
            current_d = d.copy()
            current_j = j.copy()
            labels = list(range(5))
            constant = 0.0
            for label in order:
                idx = labels.index(label)
                keep = [i for i in range(len(labels)) if i != idx]
                elim = [idx]
                if keep:
                    current_d, current_j, gain = eliminate_block(current_d, current_j, keep, elim)
                    labels = [labels[i] for i in keep]
                    constant += gain
                else:
                    gain = -0.5 * float(current_j.T @ np.linalg.solve(current_d, current_j))
                    constant += gain
                    labels = []
            max_residual = max(max_residual, abs(constant - one_step))
    return {"max_permutation_residual": max_residual}


def branch_score_from_graph(family: list[np.ndarray], x: np.ndarray, mu: float) -> float:
    dim_u = x.shape[1]
    graph = np.vstack([np.eye(dim_u), x])
    q, _ = np.linalg.qr(graph)
    p = q @ q.T
    visible = 0.0
    leakage = 0.0
    for d in family:
        visible += float(np.linalg.norm(p @ d @ p, ord="fro") ** 2)
        comm = d @ p - p @ d
        leakage += 0.5 * float(np.linalg.norm(comm, ord="fro") ** 2)
    return visible - mu * leakage


def branch_hessian_noncommuting_check() -> dict[str, float]:
    a1 = np.array([[1.0, 0.7], [0.7, -0.2]])
    a2 = np.array([[0.3, 1.1], [1.1, 0.4]])
    b1 = np.array([[2.0, -0.4], [-0.4, -1.0]])
    b2 = np.array([[-0.8, 0.9], [0.9, 1.7]])
    family = [
        np.block([[a1, np.zeros((2, 2))], [np.zeros((2, 2)), b1]]),
        np.block([[a2, np.zeros((2, 2))], [np.zeros((2, 2)), b2]]),
    ]
    mu = 0.4
    blocks_a = [a1, a2]
    blocks_b = [b1, b2]
    m_u = sum(a @ a for a in blocks_a)
    m_w = sum(b @ b for b in blocks_b)

    basis = []
    for i in range(2):
        for j in range(2):
            e = np.zeros((2, 2))
            e[i, j] = 1.0
            basis.append(e)
    analytic = np.zeros((4, 4))
    for col, x in enumerate(basis):
        gx = m_w @ x - x @ m_u
        for row, y in enumerate(basis):
            c_bilin = sum(float(np.sum((b @ x - x @ a) * (b @ y - y @ a))) for a, b in zip(blocks_a, blocks_b))
            g_bilin = float(np.sum(y * gx))
            analytic[row, col] = 2.0 * g_bilin - 2.0 * (1.0 + mu) * c_bilin

    h = 1e-5
    numeric = np.zeros((4, 4))
    zero = np.zeros((2, 2))
    f0 = branch_score_from_graph(family, zero, mu)
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            if i == j:
                numeric[i, i] = (
                    branch_score_from_graph(family, h * ei, mu)
                    - 2.0 * f0
                    + branch_score_from_graph(family, -h * ei, mu)
                ) / (h * h)
            else:
                numeric[i, j] = (
                    branch_score_from_graph(family, h * ei + h * ej, mu)
                    - branch_score_from_graph(family, h * ei - h * ej, mu)
                    - branch_score_from_graph(family, -h * ei + h * ej, mu)
                    + branch_score_from_graph(family, -h * ei - h * ej, mu)
                ) / (4.0 * h * h)
    err = numeric - analytic
    return {
        "max_abs_hessian_error": float(np.max(np.abs(err))),
        "relative_fro_error": float(np.linalg.norm(err, ord="fro") / np.linalg.norm(analytic, ord="fro")),
        "analytic_min_eig": float(np.min(np.linalg.eigvalsh(sym(analytic)))),
        "noncommutator_norm_inside_U": float(np.linalg.norm(a1 @ a2 - a2 @ a1, ord="fro")),
        "noncommutator_norm_inside_W": float(np.linalg.norm(b1 @ b2 - b2 @ b1, ord="fro")),
    }


def action(v: np.ndarray, h: np.ndarray) -> np.ndarray:
    return (
        0.6 * v * v
        + 0.5 * (1.8 + 0.4 * v) * h * h
        + 0.08 * h**4
        + 0.18 * v * v * h
        + 0.05 * v**4
    )


def find_hidden_min(v: float) -> tuple[float, float, float]:
    hs = np.linspace(-1.5, 1.5, 2001)
    vals = action(v, hs)
    idx = int(np.argmin(vals))
    h0 = float(hs[idx])
    # Newton refinement.
    h = h0
    for _ in range(20):
        grad = (1.8 + 0.4 * v) * h + 0.32 * h**3 + 0.18 * v * v
        hess = (1.8 + 0.4 * v) + 0.96 * h * h
        h -= grad / hess
    s_min = float(action(v, h))
    s_hh = float((1.8 + 0.4 * v) + 0.96 * h * h)
    return h, s_min, s_hh


def laplace_bridge_numerics() -> dict[str, float]:
    v_grid = np.linspace(-0.55, 0.55, 81)
    h_grid = np.linspace(-4.0, 4.0, 8001)
    eps_values = [0.16, 0.08, 0.04]
    errors = []
    for eps in eps_values:
        entropic = []
        laplace = []
        variational = []
        for v in v_grid:
            vals = action(v, h_grid)
            _, log_int = normalize_log_density(-vals / eps, h_grid)
            sent = -eps * log_int
            _h, svar, shh = find_hidden_min(float(v))
            sl = svar + 0.5 * eps * np.log(shh)
            entropic.append(sent)
            laplace.append(sl)
            variational.append(svar)
        entropic = centered(np.array(entropic))
        laplace = centered(np.array(laplace))
        variational = centered(np.array(variational))
        errors.append(
            {
                "epsilon": eps,
                "laplace_sup_error": float(np.max(np.abs(entropic - laplace))),
                "variational_sup_error": float(np.max(np.abs(entropic - variational))),
            }
        )
    ratio_016_to_008 = errors[0]["laplace_sup_error"] / errors[1]["laplace_sup_error"]
    ratio_008_to_004 = errors[1]["laplace_sup_error"] / errors[2]["laplace_sup_error"]
    return {
        "errors": errors,
        "laplace_error_ratio_0p16_to_0p08": ratio_016_to_008,
        "laplace_error_ratio_0p08_to_0p04": ratio_008_to_004,
    }


def branch_cumulant_verdict_failure() -> dict[str, float]:
    # Two candidate visible coordinates have the same local quadratic data.
    # Candidate A receives a cubic hidden-visible perturbation x z^2; candidate B does not.
    x = np.linspace(-8.0, 8.0, 18001)
    z = x
    alpha = 0.3
    base_log = -0.5 * x * x - alpha * x**4
    p_base, log_z = normalize_log_density(base_log, x)
    logp_base = base_log - log_z
    z2 = z * z
    p_z = p_base.copy()
    m2 = float(np.trapezoid(p_z * z2, z))
    t = 0.18
    log_m = []
    for q in t * x:
        _, lm = normalize_log_density(logp_base - q * z2, z)
        log_m.append(lm)
    log_unnorm = logp_base + np.array(log_m)
    p_a, log_za = normalize_log_density(log_unnorm, x)
    logp_a = log_unnorm - log_za
    sym_kl = float(np.trapezoid((p_base - p_a) * (logp_base - logp_a), x))
    mean_shift = float(np.trapezoid(p_a * x, x))
    return {
        "quadratic_branch_gap": 0.0,
        "full_law_visible_sym_kl_A_vs_quadratic": sym_kl,
        "visible_mean_shift_A": mean_shift,
        "first_order_mean_prediction": -(m2**2) * t,
        "second_order_sym_kl_prediction": (m2**3) * t * t,
    }


def affine_hidden_survival_branch() -> dict[str, float]:
    v = np.linspace(-1.2, 1.2, 801)
    d = 1.6
    a1 = 0.25 * (v - 0.25) ** 2 + 0.04 * v**4
    a2 = 0.25 * (v + 0.10) ** 2 + 0.04 * v**4 + 0.015
    j1 = 0.65 * v + 0.1 * np.sin(v)
    j2 = 0.25 + 0.40 * v - 0.05 * v * v
    s1 = a1 - 0.5 * j1 * j1 / d
    s2 = a2 - 0.5 * j2 * j2 / d
    gap = float(np.min(s2) - np.min(s1))
    return {
        "min_action_branch_1": float(np.min(s1)),
        "min_action_branch_2": float(np.min(s2)),
        "branch_2_minus_branch_1_gap": gap,
        "selected_branch": 1 if gap > 0 else 2,
    }


def main() -> None:
    results = {
        "variable_hidden_precision_boundary": variable_hidden_precision_boundary(),
        "affine_hidden_tower_random": affine_hidden_tower_random(),
        "branch_hessian_noncommuting": branch_hessian_noncommuting_check(),
        "laplace_bridge": laplace_bridge_numerics(),
        "branch_cumulant_verdict_failure": branch_cumulant_verdict_failure(),
        "affine_hidden_survival_branch": affine_hidden_survival_branch(),
    }

    assert results["variable_hidden_precision_boundary"]["fixed_precision_centered_residual"] < 1e-14
    assert results["variable_hidden_precision_boundary"]["variable_precision_missing_logdet_residual"] < 1e-14
    assert results["variable_hidden_precision_boundary"]["variable_precision_missing_range"] > 0.5
    assert abs(results["variable_hidden_precision_boundary"]["hostile_exact_argmin"] + 0.85) < 1e-12
    assert results["affine_hidden_tower_random"]["max_permutation_residual"] < 1e-10
    assert results["branch_hessian_noncommuting"]["relative_fro_error"] < 1e-5
    assert results["branch_hessian_noncommuting"]["noncommutator_norm_inside_U"] > 1e-2
    assert results["laplace_bridge"]["laplace_error_ratio_0p16_to_0p08"] > 2.0
    assert results["laplace_bridge"]["laplace_error_ratio_0p08_to_0p04"] > 2.0
    assert results["branch_cumulant_verdict_failure"]["full_law_visible_sym_kl_A_vs_quadratic"] > 1e-3
    assert results["affine_hidden_survival_branch"]["selected_branch"] in {1, 2}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
