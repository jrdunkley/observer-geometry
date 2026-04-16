"""Planning stress tests for 0.3.1 theory candidates.

This script deliberately does not touch the module implementation.  It probes
the theorem claims that are candidates for later implementation authorization:

- intrinsic vs ceiling-mediated ensemble separation;
- variable-precision affine-hidden elimination and tower composition;
- weighted-family closure/frontier identities and exact-branch Hessian signs;
- typed residual margin and branch-drift bounds;
- support-event charge scaling;
- positive packet lift and rank-one obstruction;
- empirical micro-real bundles as finite weighted families.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit" / "outputs" / "0_3_1_planning_stress_tests.json"


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def herm(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.conj().T)


def make_spd(rng: np.random.Generator, n: int, floor: float = 0.5) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return sym(a.T @ a / n + floor * np.eye(n))


def phi_c(h: np.ndarray, c: np.ndarray) -> np.ndarray:
    return sym(np.linalg.inv(c @ np.linalg.solve(h, c.T)))


def graph_projector(x: np.ndarray) -> np.ndarray:
    graph = np.vstack([np.eye(x.shape[1]), x])
    q, _ = np.linalg.qr(graph)
    return q @ q.T


def weighted_score(family: list[np.ndarray], weights: np.ndarray, x: np.ndarray, mu: float) -> float:
    p = graph_projector(x)
    visible = 0.0
    leakage = 0.0
    for w, a in zip(weights, family):
        visible += float(w * np.linalg.norm(p @ a @ p, ord="fro") ** 2)
        comm = a @ p - p @ a
        leakage += float(w * 0.5 * np.linalg.norm(comm, ord="fro") ** 2)
    return visible - mu * leakage


def normalize_grid_log_weight(logw: np.ndarray, grid: np.ndarray) -> tuple[np.ndarray, float]:
    shift = float(np.max(logw))
    w = np.exp(logw - shift)
    z = float(np.trapezoid(w, grid))
    return w / z, shift + float(np.log(z))


def finite_difference_second(f: Callable[[float], float], x: float, h: float = 1e-5) -> float:
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)


@dataclass(frozen=True)
class StressResult:
    passed: bool
    metrics: dict[str, object]


def no_canonical_ceiling_stress() -> StressResult:
    h = np.array([[2.0, 1.0], [1.0, 2.0]])
    c = np.array([[1.0, 0.0]])
    shears = np.linspace(-2.0, 2.0, 21)
    phis = []
    blocks = []
    for b in shears:
        s = np.array([[1.0, 0.0], [b, 1.0]])
        hb = s.T @ h @ s
        phis.append(float(phi_c(hb, c)[0, 0]))
        blocks.append(float(hb[0, 0]))
    metrics = {
        "phi_range": float(max(phis) - min(phis)),
        "coordinate_ceiling_block_range": float(max(blocks) - min(blocks)),
        "min_coordinate_block": float(min(blocks)),
        "max_coordinate_block": float(max(blocks)),
    }
    return StressResult(
        passed=metrics["phi_range"] < 1e-12 and metrics["coordinate_ceiling_block_range"] > 1.0,
        metrics=metrics,
    )


def intrinsic_covariance_stress(seed: int = 3101) -> StressResult:
    rng = np.random.default_rng(seed)
    max_latent_error = 0.0
    max_visible_error = 0.0
    for _ in range(100):
        n = 6
        m = 3
        h = make_spd(rng, n, floor=1.0)
        c = rng.normal(size=(m, n))
        while np.linalg.matrix_rank(c) < m:
            c = rng.normal(size=(m, n))
        s = rng.normal(size=(n, n))
        while abs(np.linalg.det(s)) < 0.1:
            s = rng.normal(size=(n, n))
        u = rng.normal(size=(m, m))
        while abs(np.linalg.det(u)) < 0.1:
            u = rng.normal(size=(m, m))
        phi0 = phi_c(h, c)
        max_latent_error = max(max_latent_error, float(np.linalg.norm(phi_c(s.T @ h @ s, c @ s) - phi0, ord="fro")))
        moved = phi_c(h, u @ c)
        expected = np.linalg.inv(u).T @ phi0 @ np.linalg.inv(u)
        max_visible_error = max(max_visible_error, float(np.linalg.norm(moved - expected, ord="fro")))
    metrics = {
        "trials": 100,
        "max_latent_basis_error": max_latent_error,
        "max_visible_basis_error": max_visible_error,
    }
    return StressResult(passed=max_latent_error < 1e-10 and max_visible_error < 1e-10, metrics=metrics)


def variable_precision_elimination_stress() -> StressResult:
    v_grid = np.linspace(-0.75, 0.75, 301)
    h_grid = np.linspace(-10.0, 10.0, 20001)

    def a(v: np.ndarray | float) -> np.ndarray | float:
        return 0.08 * np.asarray(v) ** 4 + 0.03 * np.asarray(v)

    def j(v: np.ndarray | float) -> np.ndarray | float:
        return 0.2 + 0.45 * np.sin(1.3 * np.asarray(v)) + 0.1 * np.asarray(v) ** 2

    def d(v: np.ndarray | float) -> np.ndarray | float:
        return 1.4 * np.exp(0.8 * np.asarray(v) - 0.2 * np.asarray(v) ** 2)

    exact_formula = a(v_grid) + 0.5 * np.log(d(v_grid)) - 0.5 * j(v_grid) ** 2 / d(v_grid)
    numeric = []
    for v in v_grid:
        action = a(v) + 0.5 * d(v) * h_grid**2 + j(v) * h_grid
        _, log_int = normalize_grid_log_weight(-action, h_grid)
        numeric.append(-log_int)
    numeric = np.asarray(numeric)
    centered_residual = (numeric - np.mean(numeric)) - (exact_formula - np.mean(exact_formula))

    # Hostile branch: variational action is exactly flat but full law selects
    # the smaller hidden fibre precision through the log-det term.
    d_left = 0.1
    d_right = 10.0
    variational_gap = 0.0
    exact_gap_left_minus_right = 0.5 * np.log(d_left / d_right)

    # Scalar derivative finite difference for the log-det and quadratic terms.
    x0 = 0.17
    d1 = float(d(x0) * (0.8 - 0.4 * x0))
    d2 = float(d(x0) * ((0.8 - 0.4 * x0) ** 2 - 0.4))
    j0 = float(j(x0))
    j1 = float(0.45 * 1.3 * np.cos(1.3 * x0) + 0.2 * x0)
    j2 = float(-0.45 * 1.3 * 1.3 * np.sin(1.3 * x0) + 0.2)
    a2 = float(0.96 * x0 * x0)
    k0 = 1.0 / float(d(x0))
    second_logdet_formula = 0.5 * (k0 * d2 - k0 * d1 * k0 * d1)
    second_var_formula = -j2 * k0 * j0 - j1 * k0 * j1 + 2.0 * j1 * k0 * d1 * k0 * j0 + 0.5 * j0 * k0 * d2 * k0 * j0 - j0 * k0 * d1 * k0 * d1 * k0 * j0
    second_formula = a2 + second_logdet_formula + second_var_formula
    second_numeric = finite_difference_second(lambda x: float(a(x) + 0.5 * np.log(d(x)) - 0.5 * j(x) ** 2 / d(x)), x0)

    metrics = {
        "centered_numeric_integral_residual": float(np.max(np.abs(centered_residual))),
        "hostile_variational_gap": variational_gap,
        "hostile_exact_left_minus_right_gap": float(exact_gap_left_minus_right),
        "second_derivative_formula": float(second_formula),
        "second_derivative_numeric": float(second_numeric),
        "second_derivative_abs_error": float(abs(second_formula - second_numeric)),
    }
    return StressResult(
        passed=metrics["centered_numeric_integral_residual"] < 1e-9
        and exact_gap_left_minus_right < -1.0
        and metrics["second_derivative_abs_error"] < 1e-5,
        metrics=metrics,
    )


def eliminate_block(d: np.ndarray, j: np.ndarray, keep: list[int], elim: list[int]) -> tuple[np.ndarray, np.ndarray, float]:
    d_kk = d[np.ix_(keep, keep)]
    d_ke = d[np.ix_(keep, elim)]
    d_ek = d[np.ix_(elim, keep)]
    d_ee = d[np.ix_(elim, elim)]
    j_k = j[keep]
    j_e = j[elim]
    solved_dek = np.linalg.solve(d_ee, d_ek)
    solved_je = np.linalg.solve(d_ee, j_e)
    new_d = sym(d_kk - d_ke @ solved_dek)
    new_j = j_k - d_ke @ solved_je
    gain = 0.5 * float(np.linalg.slogdet(d_ee)[1]) - 0.5 * float(j_e.T @ solved_je)
    return new_d, new_j, gain


def variable_precision_tower_stress(seed: int = 3102) -> StressResult:
    rng = np.random.default_rng(seed)
    max_residual = 0.0
    for _ in range(100):
        n = 6
        d = make_spd(rng, n, floor=1.2)
        j = rng.normal(size=n)
        one_step = 0.5 * float(np.linalg.slogdet(d)[1]) - 0.5 * float(j.T @ np.linalg.solve(d, j))
        # A hostile staged order, not just contiguous blocks.
        order = [4, 1, 5, 0, 3, 2]
        labels = list(range(n))
        current_d = d.copy()
        current_j = j.copy()
        staged = 0.0
        for label in order:
            idx = labels.index(label)
            keep = [i for i in range(len(labels)) if i != idx]
            if keep:
                current_d, current_j, gain = eliminate_block(current_d, current_j, keep, [idx])
                labels = [labels[i] for i in keep]
                staged += gain
            else:
                staged += 0.5 * float(np.linalg.slogdet(current_d)[1]) - 0.5 * float(current_j.T @ np.linalg.solve(current_d, current_j))
        max_residual = max(max_residual, abs(staged - one_step))
    metrics = {"trials": 100, "max_staged_vs_one_step_residual": max_residual}
    return StressResult(passed=max_residual < 1e-10, metrics=metrics)


def weighted_family_energy_stress(seed: int = 3103) -> StressResult:
    rng = np.random.default_rng(seed)
    max_split_error = 0.0
    max_stationarity_error = 0.0
    for _ in range(100):
        n = 7
        m = 3
        k = 5
        family = [sym(rng.normal(size=(n, n))) for _ in range(k)]
        weights = rng.uniform(0.1, 2.0, size=k)
        q, _ = np.linalg.qr(rng.normal(size=(n, m)))
        p = q @ q.T
        m_nu = sum(w * a @ a for w, a in zip(weights, family))
        lhs = float(np.trace(p @ m_nu))
        s_val = sum(float(w * np.linalg.norm(p @ a @ p, ord="fro") ** 2) for w, a in zip(weights, family))
        l_val = sum(float(w * 0.5 * np.linalg.norm(a @ p - p @ a, ord="fro") ** 2) for w, a in zip(weights, family))
        max_split_error = max(max_split_error, abs(lhs - (s_val + l_val)))

        # If all family members are block diagonal in P, stationarity must be exact.
        blocks_u = [sym(rng.normal(size=(m, m))) for _ in range(k)]
        blocks_w = [sym(rng.normal(size=(n - m, n - m))) for _ in range(k)]
        exact_family = [np.block([[au, np.zeros((m, n - m))], [np.zeros((n - m, m)), aw]]) for au, aw in zip(blocks_u, blocks_w)]
        p0 = np.diag([1.0] * m + [0.0] * (n - m))
        m_exact = sum(w * a @ a for w, a in zip(weights, exact_family))
        bracket_source = m_exact - sum(weights[i] * (exact_family[i] @ (exact_family[i] @ p0 - p0 @ exact_family[i]) - (exact_family[i] @ p0 - p0 @ exact_family[i]) @ exact_family[i]) for i in range(k))
        max_stationarity_error = max(max_stationarity_error, float(np.linalg.norm(p0 @ bracket_source - bracket_source @ p0, ord="fro")))
    metrics = {
        "trials": 100,
        "max_energy_split_error": max_split_error,
        "max_exact_branch_stationarity_error": max_stationarity_error,
    }
    return StressResult(passed=max_split_error < 1e-10 and max_stationarity_error < 1e-10, metrics=metrics)


def branch_hessian_status_stress() -> StressResult:
    # U is one-dimensional, W is two-dimensional. The two W directions are set
    # to opposite Hessian statuses, creating a clean saddle witness.
    family = [np.diag([1.0, 0.0, 3.0])]
    weights = np.array([1.0])
    mu = 0.0
    a = 1.0
    bs = np.array([0.0, 3.0])
    h_contract_diag = (1.0 + mu) * (bs - a) ** 2 - (bs**2 - a**2)
    analytic_second_diag = -2.0 * h_contract_diag
    h = 2e-5
    numeric_second_diag = []
    for row in range(2):
        x = np.zeros((2, 1))
        x[row, 0] = 1.0
        f0 = weighted_score(family, weights, np.zeros((2, 1)), mu)
        numeric = (weighted_score(family, weights, h * x, mu) - 2.0 * f0 + weighted_score(family, weights, -h * x, mu)) / (h * h)
        numeric_second_diag.append(float(numeric))
    numeric_second_diag = np.asarray(numeric_second_diag)
    metrics = {
        "contract_hessian_diag": h_contract_diag.tolist(),
        "analytic_second_variation_diag": analytic_second_diag.tolist(),
        "numeric_second_variation_diag": numeric_second_diag.tolist(),
        "max_second_variation_error": float(np.max(np.abs(numeric_second_diag - analytic_second_diag))),
        "status": "mixed" if np.min(h_contract_diag) < 0.0 < np.max(h_contract_diag) else "not_mixed",
    }
    return StressResult(
        passed=metrics["status"] == "mixed" and metrics["max_second_variation_error"] < 1e-5,
        metrics=metrics,
    )


def residual_certification_stress() -> StressResult:
    quadratic_gap = 0.08
    safe_r = 0.039
    threshold_r = 0.04
    unsafe_r = 0.041
    safe_constructive_gap = quadratic_gap - 2.0 * safe_r
    threshold_constructive_gap = quadratic_gap - 2.0 * threshold_r
    unsafe_constructive_gap = quadratic_gap - 2.0 * unsafe_r

    # Branch drift: J_model(x)=0.5 lambda x^2, residual r(x)=eta x.
    lam = 3.0
    eta = 0.07
    true_shift = eta / lam
    theorem_bound = 2.0 * eta / lam
    metrics = {
        "safe_worst_case_gap": safe_constructive_gap,
        "threshold_worst_case_gap": threshold_constructive_gap,
        "unsafe_worst_case_gap": unsafe_constructive_gap,
        "drift_true_shift": true_shift,
        "drift_theorem_bound": theorem_bound,
        "drift_bound_ratio": true_shift / theorem_bound,
    }
    return StressResult(
        passed=safe_constructive_gap > 0.0
        and abs(threshold_constructive_gap) < 1e-14
        and unsafe_constructive_gap < 0.0
        and true_shift <= theorem_bound,
        metrics=metrics,
    )


def support_event_charge_stress() -> StressResult:
    m_sup = 3
    rank = 2
    s = np.geomspace(1e-5, 1e-1, 80)
    tau = -0.5 * m_sup * rank * np.log(s) + 0.17 * s
    slope, intercept = np.polyfit(-np.log(s), tau, 1)
    charge = 0.5 * m_sup * rank
    metrics = {
        "m_sup": m_sup,
        "rank": rank,
        "fitted_charge": float(slope),
        "theory_charge": float(charge),
        "charge_abs_error": float(abs(slope - charge)),
        "intercept": float(intercept),
    }
    return StressResult(passed=metrics["charge_abs_error"] < 1e-3, metrics=metrics)


def packet_lift_stress(seed: int = 3104) -> StressResult:
    rng = np.random.default_rng(seed)
    max_formula_error = 0.0
    max_floor_support_error = 0
    monotone_violations = 0
    for _ in range(100):
        m = 4
        hidden = 5
        r = 3
        eps = 0.2
        v = rng.normal(size=(m, r)) + 1j * rng.normal(size=(m, r))
        w = rng.normal(size=(hidden, r)) + 1j * rng.normal(size=(hidden, r))
        psi = np.vstack([v, w])
        h = eps * np.eye(m + hidden, dtype=complex) + psi @ psi.conj().T
        a = h[:m, :m]
        b = h[:m, m:]
        d = h[m:, m:]
        schur = herm(a - b @ np.linalg.solve(d, b.conj().T))
        formula = herm(eps * np.eye(m, dtype=complex) + v @ np.linalg.solve(np.eye(r) + (w.conj().T @ w) / eps, v.conj().T))
        max_formula_error = max(max_formula_error, float(np.linalg.norm(schur - formula, ord="fro")))
        correction = formula - eps * np.eye(m, dtype=complex)
        support_rank = int(np.sum(np.linalg.svd(correction, compute_uv=False) > 1e-9))
        max_floor_support_error = max(max_floor_support_error, abs(support_rank - np.linalg.matrix_rank(v, tol=1e-9)))

        low_w = 0.2 * w
        high_w = 2.0 * w
        low = herm(v @ np.linalg.solve(np.eye(r) + (low_w.conj().T @ low_w) / eps, v.conj().T))
        high = herm(v @ np.linalg.solve(np.eye(r) + (high_w.conj().T @ high_w) / eps, v.conj().T))
        if np.min(np.linalg.eigvalsh(low - high)) < -1e-9:
            monotone_violations += 1

    v1 = rng.normal(size=(m, 1)) + 1j * rng.normal(size=(m, 1))
    w1 = rng.normal(size=(hidden, 1)) + 1j * rng.normal(size=(hidden, 1))
    rank_one_correction = v1 @ np.linalg.solve(np.eye(1) + (w1.conj().T @ w1) / 0.2, v1.conj().T)
    rank_one_rank = int(np.sum(np.linalg.svd(rank_one_correction, compute_uv=False) > 1e-9))
    metrics = {
        "trials": 100,
        "max_formula_error": max_formula_error,
        "max_floor_support_rank_error": int(max_floor_support_error),
        "monotone_screening_violations": monotone_violations,
        "rank_one_correction_rank": rank_one_rank,
    }
    return StressResult(
        passed=max_formula_error < 1e-10
        and max_floor_support_error == 0
        and monotone_violations == 0
        and rank_one_rank <= 1,
        metrics=metrics,
    )


def read_numeric_csv(path: Path, skip_first: bool = True) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            labels.append(row[0])
            rows.append([float(x) for x in row[1 if skip_first else 0 :]])
    return labels, np.asarray(rows, dtype=float)


def empirical_covariance_reference(data: np.ndarray, ridge: float = 0.1) -> np.ndarray:
    centered = data - np.mean(data, axis=0, keepdims=True)
    cov = np.cov(centered, rowvar=False, bias=False)
    return sym(cov + ridge * np.eye(cov.shape[0]))


def empirical_weighted_family_stress() -> StressResult:
    bundle_paths = {
        "iris": ROOT / "evidence" / "micro_real_bundles" / "iris_protocol_mismatch" / "raw" / "iris_slice.csv",
        "leaderboard": ROOT / "evidence" / "micro_real_bundles" / "leaderboard_benchmark_slice" / "raw" / "score_slice.csv",
    }
    metrics: dict[str, object] = {}
    passed = True
    for name, path in bundle_paths.items():
        labels, data = read_numeric_csv(path)
        h = empirical_covariance_reference(data)
        centered = data - np.mean(data, axis=0, keepdims=True)
        family = [np.outer(row, row) - h for row in centered]
        weights = np.full(len(family), 1.0 / len(family))
        m = min(2, h.shape[0] - 1)
        eigvals, eigvecs = np.linalg.eigh(sum(w * a @ a for w, a in zip(weights, family)))
        q = eigvecs[:, np.argsort(eigvals)[-m:]]
        p = q @ q.T
        m_nu = sum(w * a @ a for w, a in zip(weights, family))
        lhs = float(np.trace(p @ m_nu))
        s_val = sum(float(w * np.linalg.norm(p @ a @ p, ord="fro") ** 2) for w, a in zip(weights, family))
        l_val = sum(float(w * 0.5 * np.linalg.norm(a @ p - p @ a, ord="fro") ** 2) for w, a in zip(weights, family))
        err = abs(lhs - (s_val + l_val))
        metrics[name] = {
            "dimension": int(h.shape[0]),
            "samples": int(len(family)),
            "top_M_energy_split_error": err,
            "top_M_visibility": s_val,
            "top_M_leakage": l_val,
        }
        passed = passed and err < 1e-9

    bell_path = ROOT / "evidence" / "micro_real_bundles" / "bell_counts_bundle" / "raw" / "counts_by_context.csv"
    rows: list[list[float]] = []
    with bell_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            counts = np.array([float(row["n_pp"]), float(row["n_pm"]), float(row["n_mp"]), float(row["n_mm"])])
            rows.append((counts / np.sum(counts)).tolist())
    data = np.asarray(rows, dtype=float)
    h = empirical_covariance_reference(data, ridge=0.05)
    centered = data - np.mean(data, axis=0, keepdims=True)
    family = [np.outer(row, row) - h for row in centered]
    weights = np.full(len(family), 1.0 / len(family))
    p = np.diag([1.0, 1.0, 0.0, 0.0])
    m_nu = sum(w * a @ a for w, a in zip(weights, family))
    lhs = float(np.trace(p @ m_nu))
    s_val = sum(float(w * np.linalg.norm(p @ a @ p, ord="fro") ** 2) for w, a in zip(weights, family))
    l_val = sum(float(w * 0.5 * np.linalg.norm(a @ p - p @ a, ord="fro") ** 2) for w, a in zip(weights, family))
    err = abs(lhs - (s_val + l_val))
    metrics["bell"] = {
        "dimension": int(h.shape[0]),
        "contexts": int(len(family)),
        "fixed_context_projector_energy_split_error": err,
        "fixed_context_visibility": s_val,
        "fixed_context_leakage": l_val,
    }
    passed = passed and err < 1e-9
    return StressResult(passed=passed, metrics=metrics)


def main() -> None:
    tests = {
        "no_canonical_ceiling": no_canonical_ceiling_stress(),
        "intrinsic_covariance": intrinsic_covariance_stress(),
        "variable_precision_elimination": variable_precision_elimination_stress(),
        "variable_precision_tower": variable_precision_tower_stress(),
        "weighted_family_energy": weighted_family_energy_stress(),
        "branch_hessian_status": branch_hessian_status_stress(),
        "residual_certification": residual_certification_stress(),
        "support_event_charge": support_event_charge_stress(),
        "packet_lift": packet_lift_stress(),
        "empirical_weighted_family": empirical_weighted_family_stress(),
    }
    result = {name: asdict(value) for name, value in tests.items()}
    result["overall_passed"] = all(value.passed for value in tests.values())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["overall_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
