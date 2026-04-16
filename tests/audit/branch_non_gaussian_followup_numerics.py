from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "branch_non_gaussian_followup_numerics.json"


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def projector_from_graph(x: float) -> np.ndarray:
    v = np.array([1.0, x], dtype=float)
    v = v / np.linalg.norm(v)
    return np.outer(v, v)


def rank_one_scores(alpha: np.ndarray, beta: np.ndarray, x: float) -> tuple[float, float]:
    p = projector_from_graph(x)
    visible = 0.0
    leakage = 0.0
    for a, b in zip(alpha, beta):
        d = np.diag([float(a), float(b)])
        visible += float(np.linalg.norm(p @ d @ p, ord="fro") ** 2)
        comm = d @ p - p @ d
        leakage += 0.5 * float(np.linalg.norm(comm, ord="fro") ** 2)
    return visible, leakage


def rank_one_formula(alpha: np.ndarray, beta: np.ndarray, mu: float, x: float) -> float:
    a = float(alpha @ alpha)
    b = float(beta @ beta)
    c = float(alpha @ beta)
    d = float((beta - alpha) @ (beta - alpha))
    y = x * x
    return (a + (2.0 * c - mu * d) * y + b * y * y) / ((1.0 + y) ** 2)


def rank_one_hessian_formula(alpha: np.ndarray, beta: np.ndarray, mu: float) -> float:
    a = float(alpha @ alpha)
    b = float(beta @ beta)
    d = float((beta - alpha) @ (beta - alpha))
    return 2.0 * (b - a - (1.0 + mu) * d)


def finite_second_derivative(f, h: float = 1e-5) -> float:
    return float((f(h) - 2.0 * f(0.0) + f(-h)) / (h * h))


def branch_operator_gap(t: float, mu: float, branch: str) -> np.ndarray:
    d1 = np.diag([t, -t, 3.0, -3.0])
    d2 = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 2.0, 0.0],
        ],
        dtype=float,
    )
    if branch == "U":
        u = [0, 1]
        w = [2, 3]
    elif branch == "V":
        u = [2, 3]
        w = [0, 1]
    else:
        raise ValueError(branch)
    family = [d1, d2]
    a_blocks = [d[np.ix_(u, u)] for d in family]
    b_blocks = [d[np.ix_(w, w)] for d in family]
    m_u = sum(a @ a for a in a_blocks)
    m_w = sum(b @ b for b in b_blocks)

    basis = []
    for i in range(2):
        for j in range(2):
            e = np.zeros((2, 2), dtype=float)
            e[i, j] = 1.0
            basis.append(e)

    mat = np.zeros((4, 4), dtype=float)
    for col, x in enumerate(basis):
        gx = m_w @ x - x @ m_u
        for row, y in enumerate(basis):
            c_bilin = sum(float(np.sum((b @ x - x @ a) * (b @ y - y @ a))) for a, b in zip(a_blocks, b_blocks))
            g_bilin = float(np.sum(y * gx))
            mat[row, col] = (1.0 + mu) * c_bilin - g_bilin
    return sym(mat)


def branch_scores_4x4(t: float) -> tuple[float, float]:
    s_u = 2.0 * t * t + 2.0
    s_v = 26.0
    return s_u, s_v


def normalized_density_on_grid(logw: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, float]:
    shift = float(np.max(logw))
    w = np.exp(logw - shift)
    z_scaled = float(np.trapezoid(w, x))
    p = w / z_scaled
    log_z = shift + np.log(z_scaled)
    return p, log_z


def quartic_density_stats(h: float, alpha: float, x: np.ndarray) -> dict[str, float | np.ndarray]:
    logw = -0.5 * h * x * x - alpha * x**4
    p, log_z = normalized_density_on_grid(logw, x)
    second = float(np.trapezoid(p * x * x, x))
    fourth = float(np.trapezoid(p * x**4, x))
    return {"p": p, "log_z": log_z, "second": second, "var_x2": fourth - second * second}


def same_shape_quartic_checks() -> dict[str, float]:
    x = np.linspace(-8.0, 8.0, 24001)
    alpha = 0.17
    h0 = 1.1
    h1 = 1.38
    s0 = quartic_density_stats(h0, alpha, x)
    s1 = quartic_density_stats(h1, alpha, x)
    p0 = s0["p"]
    p1 = s1["p"]
    logp0 = np.log(p0)
    logp1 = np.log(p1)
    d_sym_numeric = float(np.trapezoid((p0 - p1) * (logp0 - logp1), x))
    d_sym_formula = 0.5 * (h1 - h0) * (float(s0["second"]) - float(s1["second"]))

    eps = 0.025
    se = quartic_density_stats(h0 + eps, alpha, x)
    pe = se["p"]
    d_sym_eps = float(np.trapezoid((p0 - pe) * (logp0 - np.log(pe)), x))
    d_sym_quad = 0.25 * eps * eps * float(s0["var_x2"])
    return {
        "sym_kl_numeric": d_sym_numeric,
        "sym_kl_formula": d_sym_formula,
        "sym_kl_abs_error": abs(d_sym_numeric - d_sym_formula),
        "small_eps_numeric": d_sym_eps,
        "small_eps_quadratic_prediction": d_sym_quad,
        "small_eps_relative_error": abs(d_sym_eps - d_sym_quad) / d_sym_eps,
    }


def affine_hidden_check() -> dict[str, float]:
    v_grid = np.linspace(-2.5, 2.5, 301)
    h_grid = np.linspace(-8.0, 8.0, 12001)
    d = 1.7

    def a_of_v(v: np.ndarray) -> np.ndarray:
        return 0.2 * v**4 + 0.35 * v**2 + 0.07 * v

    def j_of_v(v: np.ndarray) -> np.ndarray:
        return 0.6 * np.sin(v) + 0.15 * v * v - 0.05 * v

    numeric_neg = []
    formula_neg = []
    for v in v_grid:
        log_integrand = -float(a_of_v(np.array([v]))[0]) - 0.5 * d * h_grid * h_grid - float(j_of_v(np.array([v]))[0]) * h_grid
        _, log_int = normalized_density_on_grid(log_integrand, h_grid)
        numeric_neg.append(-log_int)
        formula_neg.append(float(a_of_v(np.array([v]))[0]) - 0.5 * float(j_of_v(np.array([v]))[0]) ** 2 / d)
    numeric_neg = np.array(numeric_neg)
    formula_neg = np.array(formula_neg)
    residual = numeric_neg - formula_neg
    residual = residual - float(np.mean(residual))
    return {
        "max_centered_visible_action_error": float(np.max(np.abs(residual))),
        "rms_centered_visible_action_error": float(np.sqrt(np.mean(residual**2))),
    }


def variational_quartic_check() -> dict[str, float]:
    phi = 1.3
    d = 2.4
    t = -0.9
    u = 0.7
    v = np.linspace(-0.08, 0.08, 161)
    h_star = -0.5 * t * v * v / d
    s_var = 0.5 * phi * v * v + 0.5 * d * h_star * h_star + 0.5 * t * v * v * h_star + (u / 24.0) * v**4
    y = s_var - 0.5 * phi * v * v
    coeff = float(np.linalg.lstsq((v**4)[:, None], y, rcond=None)[0][0])
    expected = u / 24.0 - (t * t) / (8.0 * d)
    return {
        "fitted_v4_coefficient": coeff,
        "expected_v4_coefficient": expected,
        "abs_error": abs(coeff - expected),
    }


def cubic_fibre_cumulant_check() -> dict[str, float]:
    x = np.linspace(-8.0, 8.0, 22001)
    z = x
    alpha = 0.3
    base = quartic_density_stats(1.0, alpha, x)
    p_x = base["p"]
    p_z = base["p"]
    logp_x = -0.5 * x * x - alpha * x**4 - float(base["log_z"])
    logp_z = -0.5 * z * z - alpha * z**4 - float(base["log_z"])
    m2 = float(base["second"])
    coeff = m2**3
    t = 0.16

    # M(q)=E exp(-q Z^2), evaluated for q=t*x.
    z2 = z * z
    log_m_q = []
    for q in t * x:
        _, log_m = normalized_density_on_grid(logp_z - q * z2, z)
        log_m_q.append(log_m)
    log_m_q = np.array(log_m_q)
    log_unnorm = logp_x + log_m_q
    p_t, log_z_vis = normalized_density_on_grid(log_unnorm, x)
    logp_t = log_unnorm - log_z_vis
    d_sym = float(np.trapezoid((p_x - p_t) * (logp_x - logp_t), x))
    mean_t = float(np.trapezoid(p_t * x, x))

    pred_sym = coeff * t * t
    pred_mean = -(m2**2) * t
    return {
        "m2": m2,
        "sym_kl_numeric": d_sym,
        "sym_kl_cumulant_prediction": pred_sym,
        "sym_kl_relative_error": abs(d_sym - pred_sym) / d_sym,
        "visible_mean_numeric": mean_t,
        "visible_mean_first_order_prediction": pred_mean,
        "visible_mean_abs_error": abs(mean_t - pred_mean),
    }


def hidden_fibre_bifurcation_check() -> dict[str, float]:
    alpha = 0.4
    t = 0.7
    x_c = -1.0 / (2.0 * t)

    def correction(x: float) -> float:
        q = 1.0 + 2.0 * t * x
        if q >= 0.0:
            return 0.0
        return -(q * q) / (16.0 * alpha)

    h = 1e-4
    left_second = (correction(x_c - h) - 2.0 * correction(x_c - 2.0 * h) + correction(x_c - 3.0 * h)) / (h * h)
    right_second = (correction(x_c + 3.0 * h) - 2.0 * correction(x_c + 2.0 * h) + correction(x_c + h)) / (h * h)
    expected_jump = -(t * t) / (2.0 * alpha)
    return {
        "x_c": x_c,
        "left_second_derivative": float(left_second),
        "right_second_derivative": float(right_second),
        "expected_left_second_derivative": expected_jump,
        "second_derivative_jump": float(left_second - right_second),
    }


def leakage_identity_check(rng: np.random.Generator) -> dict[str, float]:
    n = 5
    m = 2
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    p = q[:, :m] @ q[:, :m].T
    family = []
    for _ in range(4):
        a = rng.normal(size=(n, n))
        family.append(sym(a))
    s = 0.0
    l = 0.0
    total = 0.0
    det_hidden = 0.0
    for d in family:
        s += float(np.linalg.norm(p @ d @ p, ord="fro") ** 2)
        comm = d @ p - p @ d
        l += 0.5 * float(np.linalg.norm(comm, ord="fro") ** 2)
        total += float(np.trace(p @ d @ d))
        det_hidden += float(np.linalg.norm(comm, ord="fro") ** 2)
    return {
        "S": s,
        "L_trace_leakage": l,
        "captured_trace": total,
        "captured_identity_abs_error": abs((s + l) - total),
        "determinant_hidden_minus_2L_abs_error": abs(det_hidden - 2.0 * l),
    }


def main() -> None:
    rng = np.random.default_rng(20260411)
    alpha = np.array([0.3, -1.4, 0.8])
    beta = np.array([1.2, 0.2, -0.4])
    mu = 0.65
    xs = np.linspace(-2.0, 2.0, 51)
    phase_errors = []
    for x in xs:
        visible, leakage = rank_one_scores(alpha, beta, float(x))
        phase_errors.append(abs((visible - mu * leakage) - rank_one_formula(alpha, beta, mu, float(x))))
    hess_numeric = finite_second_derivative(lambda z: rank_one_formula(alpha, beta, mu, z))
    hess_formula = rank_one_hessian_formula(alpha, beta, mu)

    t_cross = float(np.sqrt(12.0))
    branch_gap = {}
    for branch in ("U", "V"):
        branch_gap[branch] = {}
        for mu_value in (0.0, 1.0, 5.0, 20.0):
            eigs = np.linalg.eigvalsh(branch_operator_gap(t_cross, mu_value, branch))
            branch_gap[branch][str(mu_value)] = {
                "min_selector_gap": float(np.min(eigs)),
                "eigenvalues": [float(x) for x in eigs],
            }
    s_u_left, s_v_left = branch_scores_4x4(t_cross - 0.15)
    s_u_cross, s_v_cross = branch_scores_4x4(t_cross)
    s_u_right, s_v_right = branch_scores_4x4(t_cross + 0.15)

    results = {
        "leakage_identity": leakage_identity_check(rng),
        "rank_one_phase_portrait": {
            "max_formula_error": float(np.max(phase_errors)),
            "hessian_numeric": hess_numeric,
            "hessian_formula": hess_formula,
            "hessian_abs_error": abs(hess_numeric - hess_formula),
        },
        "four_by_four_branch_crossing": {
            "t_cross": t_cross,
            "scores_left": {"S_U": s_u_left, "S_V": s_v_left},
            "scores_cross": {"S_U": s_u_cross, "S_V": s_v_cross},
            "scores_right": {"S_U": s_u_right, "S_V": s_v_right},
            "selector_gap_eigenvalues": branch_gap,
        },
        "same_shape_quartic_family": same_shape_quartic_checks(),
        "affine_hidden_exact_sector": affine_hidden_check(),
        "variational_quartic_jet": variational_quartic_check(),
        "cubic_fibre_cumulant_pathology": cubic_fibre_cumulant_check(),
        "hidden_fibre_bifurcation": hidden_fibre_bifurcation_check(),
    }

    assert results["leakage_identity"]["captured_identity_abs_error"] < 1e-10
    assert results["leakage_identity"]["determinant_hidden_minus_2L_abs_error"] < 1e-10
    assert results["rank_one_phase_portrait"]["max_formula_error"] < 1e-12
    assert results["rank_one_phase_portrait"]["hessian_abs_error"] < 1e-5
    assert abs(s_u_cross - s_v_cross) < 1e-12
    assert s_u_left < s_v_left and s_u_right > s_v_right
    assert results["same_shape_quartic_family"]["sym_kl_abs_error"] < 1e-10
    assert results["same_shape_quartic_family"]["small_eps_relative_error"] < 0.02
    assert results["affine_hidden_exact_sector"]["max_centered_visible_action_error"] < 1e-9
    assert results["variational_quartic_jet"]["abs_error"] < 1e-12
    assert results["cubic_fibre_cumulant_pathology"]["sym_kl_relative_error"] < 0.20
    assert results["hidden_fibre_bifurcation"]["right_second_derivative"] == 0.0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
