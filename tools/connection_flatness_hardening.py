from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as la

from nomogeo import canonical_lift, local_visible_calculus, visible_precision

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tools" / "outputs"


@dataclass(frozen=True)
class ConnectionFlatnessHardeningSummary:
    seed: int
    synthetic_trials: int
    max_factorisation_residual: float
    max_phi_recovery_residual: float
    max_metric_split_residual: float
    max_q_from_dk_residual: float
    max_transition_residual: float
    max_hidden_gauge_residual: float
    max_visible_gauge_residual: float
    max_current_variation: float
    max_q_from_current_residual: float
    max_geodesic_deviation_t3_scaled: float
    max_horizontal_projection_residual: float
    max_fixed_k_free_sector_q: float
    min_fixed_k_nonhorizontal_cross_norm: float
    conditioning_acceptance: dict[str, int]
    conditioning_max_metric_residual: dict[str, float]
    empirical_transition_residuals: dict[str, float]
    empirical_current_variation: dict[str, float]
    empirical_q_from_current_residual: dict[str, float]


def _random_spd(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a + n * np.eye(n, dtype=float)


def _random_surjective(rng: np.random.Generator, m: int, n: int) -> np.ndarray:
    while True:
        c = rng.normal(size=(m, n))
        if np.linalg.matrix_rank(c) == m:
            return c


def _random_symmetric(rng: np.random.Generator, n: int, scale: float = 1.0) -> np.ndarray:
    a = rng.normal(size=(n, n))
    s = 0.5 * (a + a.T)
    return scale * s / max(np.linalg.norm(s, ord=2), 1e-15)


def _symmetrize(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def _sqrt_spd(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(la.sqrtm(matrix).real, dtype=float)


def _fr_geodesic(h0: np.ndarray, a0: np.ndarray, t: float) -> np.ndarray:
    h_half = _sqrt_spd(h0)
    return np.asarray(h_half @ la.expm(t * a0) @ h_half, dtype=float)


def _fixed_adapted_basis(c: np.ndarray) -> np.ndarray:
    right_inverse = c.T @ np.linalg.inv(c @ c.T)
    null = la.null_space(c, rcond=1e-12)
    return np.concatenate([right_inverse, null], axis=1)


def _coords_from_h(c: np.ndarray, u: np.ndarray, h: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    h_tilde = u.T @ h @ u
    m = c.shape[0]
    h_vv = h_tilde[:m, :m]
    h_vh = h_tilde[:m, m:]
    h_hh = h_tilde[m:, m:]
    phi = h_vv - h_vh @ np.linalg.solve(h_hh, h_vh.T)
    k = h_vh @ np.linalg.inv(h_hh)
    return phi, h_hh, k, h_tilde


def _h_from_coords(u: np.ndarray, phi: np.ndarray, r: np.ndarray, k: np.ndarray) -> np.ndarray:
    m = phi.shape[0]
    h = r.shape[0]
    t_k = np.block([[np.eye(m), k], [np.zeros((h, m)), np.eye(h)]])
    diag = np.block([[phi, np.zeros((m, h))], [np.zeros((h, m)), r]])
    h_tilde = t_k @ diag @ t_k.T
    u_inv = np.linalg.inv(u)
    return u_inv.T @ h_tilde @ u_inv


def _delta_from_coord_variation(
    u: np.ndarray,
    phi: np.ndarray,
    r: np.ndarray,
    k: np.ndarray,
    dphi: np.ndarray,
    dr: np.ndarray,
    dk: np.ndarray,
) -> np.ndarray:
    m = phi.shape[0]
    h = r.shape[0]
    t_k = np.block([[np.eye(m), k], [np.zeros((h, m)), np.eye(h)]])
    dt_k = np.block([[np.zeros((m, m)), dk], [np.zeros((h, m)), np.zeros((h, h))]])
    diag = np.block([[phi, np.zeros((m, h))], [np.zeros((h, m)), r]])
    ddiag = np.block([[dphi, np.zeros((m, h))], [np.zeros((h, m)), dr]])
    dh_tilde = dt_k @ diag @ t_k.T + t_k @ ddiag @ t_k.T + t_k @ diag @ dt_k.T
    u_inv = np.linalg.inv(u)
    return u_inv.T @ dh_tilde @ u_inv


def _fisher_metric(h: np.ndarray, d1: np.ndarray, d2: np.ndarray) -> float:
    h_inv = np.linalg.inv(h)
    return float(0.5 * np.trace(h_inv @ d1 @ h_inv @ d2))


def _read_numeric_csv(path: Path, skip_first: bool = True) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        _header = next(reader)
        for row in reader:
            labels.append(row[0])
            rows.append([float(x) for x in row[1 if skip_first else 0 :]])
    return labels, np.asarray(rows, dtype=float)


def _empirical_reference(data: np.ndarray, ridge: float = 1e-1) -> np.ndarray:
    centered = data - np.mean(data, axis=0, keepdims=True)
    cov = np.cov(centered, rowvar=False, bias=False)
    cov = 0.5 * (cov + cov.T)
    return cov + ridge * np.eye(cov.shape[0], dtype=float)


def _bell_probabilities() -> np.ndarray:
    path = ROOT / "evidence" / "micro_real_bundles" / "bell_counts_bundle" / "raw" / "counts_by_context.csv"
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append([float(row["n_pp"]), float(row["n_pm"]), float(row["n_mp"]), float(row["n_mm"])])
    counts = np.asarray(rows, dtype=float)
    return counts / np.sum(counts, axis=1, keepdims=True)


def _transition_residual(h: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> float:
    u1 = _fixed_adapted_basis(c1)
    u2 = _fixed_adapted_basis(c2)
    phi1, r1, k1, _ = _coords_from_h(c1, u1, h)
    phi2, r2, k2, _ = _coords_from_h(c2, u2, h)

    m = c1.shape[0]
    transform = np.linalg.inv(u1) @ u2
    a, b = transform[:m, :m], transform[:m, m:]
    c, d = transform[m:, :m], transform[m:, m:]
    c_hat = c + k1.T @ a
    d_hat = d + k1.T @ b

    a2 = a.T @ phi1 @ a + c_hat.T @ r1 @ c_hat
    b2 = a.T @ phi1 @ b + c_hat.T @ r1 @ d_hat
    d2 = b.T @ phi1 @ b + d_hat.T @ r1 @ d_hat
    phi2_formula = a2 - b2 @ np.linalg.inv(d2) @ b2.T
    r2_formula = d2
    k2_formula = b2 @ np.linalg.inv(d2)
    return max(
        float(np.linalg.norm(phi2 - phi2_formula, ord="fro")),
        float(np.linalg.norm(r2 - r2_formula, ord="fro")),
        float(np.linalg.norm(k2 - k2_formula, ord="fro")),
    )


def _current_diagnostics(h0: np.ndarray, c: np.ndarray, a0: np.ndarray) -> tuple[float, float]:
    u = _fixed_adapted_basis(c)
    ts = np.linspace(-0.25, 0.25, 11)
    currents: list[np.ndarray] = []
    max_q_residual = 0.0
    for t in ts:
        step = 1e-4
        h_minus = _fr_geodesic(h0, a0, t - step)
        h_center = _fr_geodesic(h0, a0, t)
        h_plus = _fr_geodesic(h0, a0, t + step)
        phi, r, k, _ = _coords_from_h(c, u, h_center)
        _phi_m, _r_m, k_minus, _ = _coords_from_h(c, u, h_minus)
        _phi_p, _r_p, k_plus, _ = _coords_from_h(c, u, h_plus)
        k_prime = (k_plus - k_minus) / (2.0 * step)
        j = np.linalg.inv(phi) @ k_prime @ r
        currents.append(j)
        local = local_visible_calculus(h_center, c, _symmetrize((h_plus - h_minus) / (2.0 * step)))
        q_from_j = phi @ j @ np.linalg.inv(r) @ j.T @ phi
        max_q_residual = max(max_q_residual, float(np.linalg.norm(local.Q - q_from_j, ord="fro")))

    anchor = currents[len(currents) // 2]
    max_j_var = max(float(np.linalg.norm(current - anchor, ord="fro")) for current in currents)
    return max_j_var, max_q_residual


def run_connection_flatness_hardening(seed: int = 20260409, synthetic_trials: int = 40) -> ConnectionFlatnessHardeningSummary:
    rng = np.random.default_rng(seed)
    max_factorisation = 0.0
    max_phi = 0.0
    max_metric = 0.0
    max_q = 0.0
    max_transition = 0.0
    max_hidden_gauge = 0.0
    max_visible_gauge = 0.0
    max_current_variation = 0.0
    max_q_from_current = 0.0
    max_geodesic_deviation = 0.0
    max_horizontal_projection = 0.0
    max_fixed_k_q = 0.0
    min_fixed_k_cross = float("inf")

    for idx in range(synthetic_trials):
        n = 6 if idx % 2 == 0 else 5
        m = 1 if idx % 4 == 0 else (n - 1 if idx % 4 == 1 else 2)
        h_dim = n - m
        h = _random_spd(rng, n)
        c = _random_surjective(rng, m, n)
        u = _fixed_adapted_basis(c)
        phi, r, k, _ = _coords_from_h(c, u, h)

        h_reconstructed = _h_from_coords(u, phi, r, k)
        max_factorisation = max(max_factorisation, float(np.linalg.norm(h - h_reconstructed, ord="fro")))
        max_phi = max(max_phi, float(np.linalg.norm(visible_precision(h, c) - phi, ord="fro")))

        dphi = _random_symmetric(rng, m)
        dr = _random_symmetric(rng, h_dim) if h_dim else np.zeros((0, 0), dtype=float)
        dk = rng.normal(size=(m, h_dim)) if h_dim else np.zeros((m, 0), dtype=float)
        delta = _delta_from_coord_variation(u, phi, r, k, dphi, dr, dk)
        local = local_visible_calculus(h, c, delta)
        q_expected = dk @ r @ dk.T if h_dim else np.zeros((m, m), dtype=float)
        max_q = max(max_q, float(np.linalg.norm(local.Q - q_expected, ord="fro")))
        metric_exact = _fisher_metric(h, delta, delta)
        metric_coords = (
            0.5 * np.trace(np.linalg.solve(phi, dphi) @ np.linalg.solve(phi, dphi))
            + 0.5 * np.trace(np.linalg.solve(r, dr) @ np.linalg.solve(r, dr))
            + np.trace(np.linalg.solve(phi, dk) @ r @ dk.T)
        )
        max_metric = max(max_metric, abs(metric_exact - float(metric_coords)))

        c2 = _random_surjective(rng, m, n)
        max_transition = max(max_transition, _transition_residual(h, c, c2))

        if h_dim:
            s = rng.normal(size=(h_dim, m))
            t = rng.normal(size=(h_dim, h_dim))
            while abs(np.linalg.det(t)) < 0.2:
                t = rng.normal(size=(h_dim, h_dim))
            hidden_gauge = np.block([[np.eye(m), np.zeros((m, h_dim))], [s, t]])
            u2 = u @ hidden_gauge
            phi2, r2, k2, _ = _coords_from_h(c, u2, h)
            max_hidden_gauge = max(
                max_hidden_gauge,
                float(np.linalg.norm(phi2 - phi, ord="fro")),
                float(np.linalg.norm(r2 - (t.T @ r @ t), ord="fro")),
                float(np.linalg.norm(k2 - ((k + s.T) @ np.linalg.inv(t.T)), ord="fro")),
            )

            k_fixed = rng.normal(size=(m, h_dim))
            phi0 = _random_spd(rng, m)
            r0 = _random_spd(rng, h_dim)
            a_phi = _random_symmetric(rng, m, scale=0.12)
            a_r = _random_symmetric(rng, h_dim, scale=0.08)
            phi_half = _sqrt_spd(phi0)
            r_half = _sqrt_spd(r0)
            cross_norms: list[float] = []
            for t_val in (-0.2, -0.1, 0.0, 0.1, 0.2):
                phi_t = phi_half @ la.expm(t_val * a_phi) @ phi_half
                r_t = r_half @ la.expm(t_val * a_r) @ r_half
                h_t = _h_from_coords(u, np.asarray(phi_t, dtype=float), np.asarray(r_t, dtype=float), k_fixed)
                step = 1e-4
                h_minus = _h_from_coords(
                    u,
                    np.asarray(phi_half @ la.expm((t_val - step) * a_phi) @ phi_half, dtype=float),
                    np.asarray(r_half @ la.expm((t_val - step) * a_r) @ r_half, dtype=float),
                    k_fixed,
                )
                h_plus = _h_from_coords(
                    u,
                    np.asarray(phi_half @ la.expm((t_val + step) * a_phi) @ phi_half, dtype=float),
                    np.asarray(r_half @ la.expm((t_val + step) * a_r) @ r_half, dtype=float),
                    k_fixed,
                )
                local_t = local_visible_calculus(h_t, c, _symmetrize((h_plus - h_minus) / (2.0 * step)))
                max_fixed_k_q = max(max_fixed_k_q, float(np.linalg.norm(local_t.Q, ord="fro")))
                h_tilde_prime = u.T @ ((h_plus - h_minus) / (2.0 * step)) @ u
                cross_norms.append(float(np.linalg.norm(h_tilde_prime[:m, m:], ord="fro")))
                phi_ref = np.asarray(phi_t, dtype=float)
                max_horizontal_projection = max(max_horizontal_projection, float(np.linalg.norm(visible_precision(h_t, c) - phi_ref, ord="fro")))
            min_fixed_k_cross = min(min_fixed_k_cross, max(cross_norms))

        g = rng.normal(size=(m, m))
        while abs(np.linalg.det(g)) < 0.2:
            g = rng.normal(size=(m, m))
        compatible = np.block([[np.linalg.inv(g), np.zeros((m, h_dim))], [np.zeros((h_dim, m)), np.eye(h_dim)]])
        u_g = u @ compatible
        phi_g, r_g, k_g, _ = _coords_from_h(g @ c, u_g, h)
        max_visible_gauge = max(
            max_visible_gauge,
            float(np.linalg.norm(phi_g - (np.linalg.inv(g).T @ phi @ np.linalg.inv(g)), ord="fro")),
            float(np.linalg.norm(r_g - r, ord="fro")),
            float(np.linalg.norm(k_g - (np.linalg.inv(g).T @ k), ord="fro")),
        )

        a0 = _random_symmetric(rng, n, scale=0.12)
        max_j_var, max_q_residual = _current_diagnostics(h, c, a0)
        max_current_variation = max(max_current_variation, max_j_var)
        max_q_from_current = max(max_q_from_current, max_q_residual)

        step = 0.02
        delta0 = _sqrt_spd(h) @ a0 @ _sqrt_spd(h)
        local0 = local_visible_calculus(h, c, delta0)
        phi0 = visible_precision(h, c)
        phi_half0 = _sqrt_spd(phi0)
        a_bar = np.linalg.inv(phi_half0) @ local0.V @ np.linalg.inv(phi_half0)
        phi_ref = np.asarray(phi_half0 @ la.expm(step * a_bar) @ phi_half0, dtype=float)
        phi_t = visible_precision(_fr_geodesic(h, a0, step), c)
        deviation = phi_t - phi_ref + 0.5 * (step**2) * local0.Q
        max_geodesic_deviation = max(max_geodesic_deviation, float(np.linalg.norm(deviation, ord="fro") / (step**3)))

    conditioning_acceptance: dict[str, int] = {}
    conditioning_max_metric_residual: dict[str, float] = {}
    for exponent in (4, 6, 8):
        accepted = 0
        worst_metric = 0.0
        for _ in range(12):
            n = 6
            m = 2
            spectrum = np.geomspace(1.0, 10.0**exponent, n)
            q, _ = np.linalg.qr(rng.normal(size=(n, n)))
            h = q @ np.diag(spectrum) @ q.T
            c = _random_surjective(rng, m, n)
            try:
                u = _fixed_adapted_basis(c)
                phi, r, k, _ = _coords_from_h(c, u, h)
                dphi = _random_symmetric(rng, m)
                dr = _random_symmetric(rng, n - m)
                dk = rng.normal(size=(m, n - m))
                delta = _delta_from_coord_variation(u, phi, r, k, dphi, dr, dk)
                metric_exact = _fisher_metric(h, delta, delta)
                metric_coords = (
                    0.5 * np.trace(np.linalg.solve(phi, dphi) @ np.linalg.solve(phi, dphi))
                    + 0.5 * np.trace(np.linalg.solve(r, dr) @ np.linalg.solve(r, dr))
                    + np.trace(np.linalg.solve(phi, dk) @ r @ dk.T)
                )
                worst_metric = max(worst_metric, abs(metric_exact - float(metric_coords)))
                accepted += 1
            except Exception:
                pass
        conditioning_acceptance[f"cond_1e{exponent}"] = accepted
        conditioning_max_metric_residual[f"cond_1e{exponent}"] = worst_metric

    empirical_transition_residuals: dict[str, float] = {}
    empirical_current_variation: dict[str, float] = {}
    empirical_q_from_current_residual: dict[str, float] = {}

    _iris_labels, iris_data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "iris_protocol_mismatch" / "raw" / "iris_slice.csv")
    iris_h = _empirical_reference(iris_data)
    empirical_transition_residuals["iris"] = _transition_residual(
        iris_h,
        np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )
    iris_a0 = _random_symmetric(rng, 4, scale=0.1)
    iris_current = _current_diagnostics(iris_h, np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]]), iris_a0)
    empirical_current_variation["iris"] = iris_current[0]
    empirical_q_from_current_residual["iris"] = iris_current[1]

    _leader_labels, leader_data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "leaderboard_benchmark_slice" / "raw" / "score_slice.csv")
    leader_h = _empirical_reference(leader_data)
    empirical_transition_residuals["leaderboard"] = _transition_residual(
        leader_h,
        np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
        np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )
    leader_current = _current_diagnostics(leader_h, np.array([[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 1.0, 0.0]]), _random_symmetric(rng, 4, scale=0.1))
    empirical_current_variation["leaderboard"] = leader_current[0]
    empirical_q_from_current_residual["leaderboard"] = leader_current[1]

    bell_probs = _bell_probabilities()
    bell_h = _empirical_reference(bell_probs, ridge=5e-2)
    empirical_transition_residuals["bell"] = _transition_residual(
        bell_h,
        np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )
    bell_current = _current_diagnostics(bell_h, np.array([[1.0, -1.0, 0.0, 0.0], [0.0, 0.0, 1.0, -1.0]]), _random_symmetric(rng, 4, scale=0.06))
    empirical_current_variation["bell"] = bell_current[0]
    empirical_q_from_current_residual["bell"] = bell_current[1]

    return ConnectionFlatnessHardeningSummary(
        seed=seed,
        synthetic_trials=synthetic_trials,
        max_factorisation_residual=max_factorisation,
        max_phi_recovery_residual=max_phi,
        max_metric_split_residual=max_metric,
        max_q_from_dk_residual=max_q,
        max_transition_residual=max_transition,
        max_hidden_gauge_residual=max_hidden_gauge,
        max_visible_gauge_residual=max_visible_gauge,
        max_current_variation=max_current_variation,
        max_q_from_current_residual=max_q_from_current,
        max_geodesic_deviation_t3_scaled=max_geodesic_deviation,
        max_horizontal_projection_residual=max_horizontal_projection,
        max_fixed_k_free_sector_q=max_fixed_k_q,
        min_fixed_k_nonhorizontal_cross_norm=min_fixed_k_cross,
        conditioning_acceptance=conditioning_acceptance,
        conditioning_max_metric_residual=conditioning_max_metric_residual,
        empirical_transition_residuals=empirical_transition_residuals,
        empirical_current_variation=empirical_current_variation,
        empirical_q_from_current_residual=empirical_q_from_current_residual,
    )


def write_connection_flatness_hardening_outputs(summary: ConnectionFlatnessHardeningSummary) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "connection_flatness_hardening_summary.json"
    path.write_text(json.dumps(asdict(summary), indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    summary = run_connection_flatness_hardening()
    write_connection_flatness_hardening_outputs(summary)
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
