from __future__ import annotations

import json
from dataclasses import asdict, dataclass

import numpy as np
import scipy.linalg as la

from nomogeo import (
    clock,
    dv_bridge,
    gaussian_data_processing_contraction,
    hidden_contraction,
    hidden_load,
    load_from_hidden_contraction,
    local_visible_calculus,
    transport_hidden_load,
    visible_geometry,
    visible_from_hidden_load,
    visible_precision,
)
from nomogeo.exceptions import SupportError


@dataclass(frozen=True)
class SweepSummary:
    seed: int
    dims: list[int]
    draws_per_dim: int
    max_lift_projector_residual: float
    max_h_projector_symmetry_residual: float
    max_energy_split_residual: float
    max_tower_law_residual: float
    max_curvature_split_residual: float
    max_local_t3_normalised_residual: float
    max_hidden_rank_mismatch: int
    max_hidden_clock_residual: float
    max_hidden_loewner_violation: float
    max_inverse_lambda_roundtrip_rel_error: float
    max_inverse_x_roundtrip_rel_error: float
    max_inverse_order_violation: float
    max_transport_two_step_residual: float
    max_downstairs_associativity_residual: float
    max_long_chain_transport_residual: float
    max_long_chain_clock_residual: float
    max_local_hidden_birth_residual: float
    max_dv_visible_residual: float
    max_dv_hidden_normalized_residual: float
    max_divergence_contraction_violation: float


def _random_spd(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a + n * np.eye(n, dtype=float)


def _random_surjective(rng: np.random.Generator, m: int, n: int) -> np.ndarray:
    while True:
        c = rng.normal(size=(m, n))
        if np.linalg.matrix_rank(c) == m:
            return c


def _random_psd(rng: np.random.Generator, n: int, rank: int | None = None, scale: float = 1.0) -> np.ndarray:
    target_rank = n if rank is None else rank
    factor = rng.normal(size=(n, target_rank))
    return scale * (factor @ factor.T)


def _random_skew(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a - a.T


def _orthogonal(rng: np.random.Generator, n: int) -> np.ndarray:
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return q


def _null_space(matrix: np.ndarray) -> np.ndarray:
    return la.null_space(matrix, rcond=1e-12)


def _fro(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, ord="fro"))


def _kappa(H: np.ndarray, C: np.ndarray) -> float:
    phi = visible_precision(H, C)
    sign, value = np.linalg.slogdet(phi)
    if sign <= 0:
        raise RuntimeError("non-positive visible precision determinant in sweep")
    return float(-value)


def _psd_violation(matrix: np.ndarray) -> float:
    if matrix.size == 0:
        return 0.0
    return float(max(0.0, -np.min(np.linalg.eigvalsh(matrix))))


def run_validation_sweep(seed: int = 20260408, dims: list[int] | None = None, draws_per_dim: int = 4) -> SweepSummary:
    rng = np.random.default_rng(seed)
    dimensions = dims if dims is not None else list(range(2, 9))

    metrics: dict[str, float] = {
        "max_lift_projector_residual": 0.0,
        "max_h_projector_symmetry_residual": 0.0,
        "max_energy_split_residual": 0.0,
        "max_tower_law_residual": 0.0,
        "max_curvature_split_residual": 0.0,
        "max_local_t3_normalised_residual": 0.0,
        "max_hidden_rank_mismatch": 0.0,
        "max_hidden_clock_residual": 0.0,
        "max_hidden_loewner_violation": 0.0,
        "max_inverse_lambda_roundtrip_rel_error": 0.0,
        "max_inverse_x_roundtrip_rel_error": 0.0,
        "max_inverse_order_violation": 0.0,
        "max_transport_two_step_residual": 0.0,
        "max_downstairs_associativity_residual": 0.0,
        "max_long_chain_transport_residual": 0.0,
        "max_long_chain_clock_residual": 0.0,
        "max_local_hidden_birth_residual": 0.0,
        "max_dv_visible_residual": 0.0,
        "max_dv_hidden_normalized_residual": 0.0,
        "max_divergence_contraction_violation": 0.0,
    }

    for n in dimensions:
        for _ in range(draws_per_dim):
            m = int(rng.integers(1, n + 1))
            H = _random_spd(rng, n)
            C = _random_surjective(rng, m, n)
            geometry = visible_geometry(H, C)
            phi = geometry.phi
            L = geometry.lift
            P = geometry.projector
            N = _null_space(C)

            projector_residual = max(
                _fro(C @ L - np.eye(m)),
                _fro(C @ P),
                _fro(P @ P - P),
                _fro(P @ L),
                _fro(L @ C + P - np.eye(n)),
            )
            metrics["max_lift_projector_residual"] = max(metrics["max_lift_projector_residual"], projector_residual)
            metrics["max_h_projector_symmetry_residual"] = max(
                metrics["max_h_projector_symmetry_residual"],
                _fro(P.T @ H - H @ P),
            )

            y = rng.normal(size=m)
            z = N @ rng.normal(size=N.shape[1]) if N.shape[1] else np.zeros(n, dtype=float)
            x = L @ y + z
            energy_residual = abs(float(x.T @ H @ x - (y.T @ phi @ y + z.T @ H @ z)))
            if N.shape[1]:
                energy_residual = max(energy_residual, _fro(L.T @ H @ N))
            metrics["max_energy_split_residual"] = max(metrics["max_energy_split_residual"], energy_residual)

            stage_dims = sorted({n, max(1, n - 1), max(1, n - 2), 1}, reverse=True)
            maps = [_random_surjective(rng, stage_dims[i + 1], stage_dims[i]) for i in range(len(stage_dims) - 1)]
            composed = maps[-1]
            for map_ in reversed(maps[:-1]):
                composed = composed @ map_
            staged = H
            for map_ in maps:
                staged = visible_precision(staged, map_)
            metrics["max_tower_law_residual"] = max(
                metrics["max_tower_law_residual"],
                _fro(visible_precision(H, composed) - staged),
            )

            raw_delta = rng.normal(size=(n, n))
            Delta = 0.01 * (raw_delta + raw_delta.T) / np.linalg.norm(raw_delta + raw_delta.T, ord=2)
            local = local_visible_calculus(H, C, Delta)
            step = 1e-5
            curvature_fd = (_kappa(H + step * Delta, C) - 2.0 * _kappa(H, C) + _kappa(H - step * Delta, C)) / (step * step)
            metrics["max_curvature_split_residual"] = max(
                metrics["max_curvature_split_residual"],
                abs(float(curvature_fd - local.det_split)),
            )

            ts = np.array([1e-1, 5e-2, 2e-2, 1e-2, 5e-3], dtype=float)
            residuals = []
            for t in ts:
                phi_t = visible_precision(H + t * Delta, C)
                residual = phi_t - local.phi - t * local.V + (t * t) * local.Q
                residuals.append(_fro(residual))
            for t, residual in zip(ts, residuals, strict=True):
                metrics["max_local_t3_normalised_residual"] = max(
                    metrics["max_local_t3_normalised_residual"],
                    float(residual / (t**3)),
                )

            support_rank = int(rng.integers(0, n + 1))
            basis = _orthogonal(rng, n)[:, :support_rank]
            values = np.linspace(1.0, 2.0, support_rank)
            T = basis @ np.diag(values) @ basis.T if support_rank else np.zeros((n, n), dtype=float)
            lambda_reduced = _random_psd(rng, support_rank) if support_rank else np.zeros((0, 0), dtype=float)
            X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
            load = hidden_load(T, X)
            gap_rank = int(np.linalg.matrix_rank(T - X, tol=max(load.metadata.rank_tol, 1e-12)))
            metrics["max_hidden_rank_mismatch"] = max(metrics["max_hidden_rank_mismatch"], abs(load.rank - gap_rank))
            lambda_roundtrip = hidden_load(T, X).reduced_lambda
            x_roundtrip = visible_from_hidden_load(T, lambda_roundtrip, lambda_representation="reduced")
            lambda_denom = max(1.0, float(np.linalg.norm(lambda_reduced, ord="fro")))
            x_denom = max(1.0, float(np.linalg.norm(X, ord="fro")))
            metrics["max_inverse_lambda_roundtrip_rel_error"] = max(
                metrics["max_inverse_lambda_roundtrip_rel_error"],
                float(np.linalg.norm(lambda_roundtrip - lambda_reduced, ord="fro") / lambda_denom),
            )
            metrics["max_inverse_x_roundtrip_rel_error"] = max(
                metrics["max_inverse_x_roundtrip_rel_error"],
                float(np.linalg.norm(x_roundtrip - X, ord="fro") / x_denom),
            )
            if support_rank:
                t_s = load.support_basis.T @ T @ load.support_basis
                x_s = load.support_basis.T @ X @ load.support_basis
                clock_residual = abs(float(np.linalg.slogdet(t_s)[1] - np.linalg.slogdet(x_s)[1] - clock(load.reduced_lambda)))
                metrics["max_hidden_clock_residual"] = max(metrics["max_hidden_clock_residual"], clock_residual)

                extra = _random_psd(rng, support_rank)
                x_small = visible_from_hidden_load(T, load.reduced_lambda, lambda_representation="reduced")
                x_big = visible_from_hidden_load(T, load.reduced_lambda + extra, lambda_representation="reduced")
                lambda_small = hidden_load(T, x_small).reduced_lambda
                lambda_big = hidden_load(T, x_big).reduced_lambda
                metrics["max_hidden_loewner_violation"] = max(
                    metrics["max_hidden_loewner_violation"],
                    _psd_violation(x_small - x_big),
                    _psd_violation(lambda_big - lambda_small),
                )
                metrics["max_inverse_order_violation"] = max(
                    metrics["max_inverse_order_violation"],
                    _psd_violation(x_small - x_big),
                )

                two_a = _random_psd(rng, support_rank, scale=0.05)
                two_b = _random_psd(rng, support_rank, scale=0.05)
                transport = transport_hidden_load(two_a, two_b)
                factor_transport = load_from_hidden_contraction(hidden_contraction(two_b) @ hidden_contraction(two_a))
                metrics["max_transport_two_step_residual"] = max(
                    metrics["max_transport_two_step_residual"],
                    _fro(transport - factor_transport),
                )

                assoc_a = _random_psd(rng, support_rank, scale=0.05)
                assoc_b = _random_psd(rng, support_rank, scale=0.05)
                assoc_c = _random_psd(rng, support_rank, scale=0.05)
                left_factor = hidden_contraction(assoc_c) @ (hidden_contraction(assoc_b) @ hidden_contraction(assoc_a))
                right_factor = (hidden_contraction(assoc_c) @ hidden_contraction(assoc_b)) @ hidden_contraction(assoc_a)
                metrics["max_downstairs_associativity_residual"] = max(
                    metrics["max_downstairs_associativity_residual"],
                    _fro(left_factor - right_factor),
                    _fro(load_from_hidden_contraction(left_factor) - load_from_hidden_contraction(right_factor)),
                )

                chain = [_random_psd(rng, support_rank, scale=0.03) for _ in range(12)]
                factor_total = np.eye(support_rank, dtype=float)
                for load_i in chain:
                    factor_total = hidden_contraction(load_i) @ factor_total
                total_from_chain = load_from_hidden_contraction(factor_total)
                pi_total = factor_total.T @ factor_total
                total_from_pi = hidden_load(np.eye(support_rank), pi_total, support_mode="ambient").reduced_lambda
                metrics["max_long_chain_transport_residual"] = max(
                    metrics["max_long_chain_transport_residual"],
                    _fro(total_from_chain - total_from_pi),
                )
                metrics["max_long_chain_clock_residual"] = max(
                    metrics["max_long_chain_clock_residual"],
                    abs(float(clock(total_from_chain) - sum(clock(load_i) for load_i in chain))),
                )

            V = np.diag(np.linspace(1.1, 1.7, max(1, min(3, n))))
            Q = _random_psd(rng, V.shape[0], scale=0.2)
            evals, evecs = np.linalg.eigh(V)
            inv_sqrt_v = (evecs * (1.0 / np.sqrt(evals))) @ evecs.T
            A_birth = inv_sqrt_v @ Q @ inv_sqrt_v
            birth_residuals = []
            for t in (1e-1, 5e-2, 2e-2, 1e-2, 5e-3):
                X_t = t * V - (t * t) * Q
                T_t = t * V
                lambda_t = hidden_load(T_t, X_t, support_mode="ambient").reduced_lambda
                birth_residuals.append(_fro(lambda_t / t - A_birth))
            metrics["max_local_hidden_birth_residual"] = max(
                metrics["max_local_hidden_birth_residual"],
                min(birth_residuals),
            )

            H0 = _random_spd(rng, n)
            J1 = _random_skew(rng, n)
            C_bridge = _random_surjective(rng, max(1, min(3, n)), n)
            delta2 = dv_bridge(H0, J1).delta_dv
            local_bridge = local_visible_calculus(H0, C_bridge, delta2)
            eps = 2e-2
            bridge_eps = dv_bridge(H0, eps * J1)
            X_eps = visible_precision(bridge_eps.h_dv, C_bridge) - visible_precision(H0, C_bridge)
            metrics["max_dv_visible_residual"] = max(
                metrics["max_dv_visible_residual"],
                _fro(X_eps - (eps * eps) * local_bridge.V + (eps**4) * local_bridge.Q),
            )
            if local_bridge.active_support.shape[1]:
                v_eigs, v_vecs = np.linalg.eigh(local_bridge.V)
                keep = v_eigs > 1e-10
                v_basis = v_vecs[:, keep]
                v_s = v_basis.T @ local_bridge.V @ v_basis
                q_s = v_basis.T @ local_bridge.Q @ v_basis
                evals, evecs = np.linalg.eigh(v_s)
                inv_sqrt = (evecs * (1.0 / np.sqrt(evals))) @ evecs.T
                a_dv = inv_sqrt @ q_s @ inv_sqrt
                dv_residuals = []
                for eps_try in (2e-2, 1e-2, 5e-3, 2e-3):
                    bridge_try = dv_bridge(H0, eps_try * J1)
                    x_try = visible_precision(bridge_try.h_dv, C_bridge) - visible_precision(H0, C_bridge)
                    try:
                        lambda_eps = hidden_load((eps_try * eps_try) * local_bridge.V, x_try, support_mode="auto").reduced_lambda
                        if lambda_eps.shape == a_dv.shape:
                            dv_residuals.append(_fro(lambda_eps / (eps_try * eps_try) - a_dv))
                    except SupportError:
                        continue
                if not dv_residuals:
                    raise RuntimeError("failed to find a numerically stable DV hidden-birth scale in the sweep")
                metrics["max_dv_hidden_normalized_residual"] = max(
                    metrics["max_dv_hidden_normalized_residual"],
                    min(dv_residuals),
                )

            H1 = _random_spd(rng, n)
            H2 = _random_spd(rng, n)
            C1 = np.eye(n)
            if n == 1:
                D = np.eye(1)
            else:
                D = _random_surjective(rng, max(1, n - 1), n)
            contraction = gaussian_data_processing_contraction(H1, H2, C1, D)
            metrics["max_divergence_contraction_violation"] = max(
                metrics["max_divergence_contraction_violation"],
                max(0.0, contraction.forward_kl_coarse - contraction.forward_kl_fine),
                max(0.0, contraction.reverse_kl_coarse - contraction.reverse_kl_fine),
                max(0.0, contraction.bhattacharyya_coarse - contraction.bhattacharyya_fine),
                max(0.0, contraction.hellinger_sq_coarse - contraction.hellinger_sq_fine),
            )

    return SweepSummary(
        seed=seed,
        dims=dimensions,
        draws_per_dim=draws_per_dim,
        max_lift_projector_residual=metrics["max_lift_projector_residual"],
        max_h_projector_symmetry_residual=metrics["max_h_projector_symmetry_residual"],
        max_energy_split_residual=metrics["max_energy_split_residual"],
        max_tower_law_residual=metrics["max_tower_law_residual"],
        max_curvature_split_residual=metrics["max_curvature_split_residual"],
        max_local_t3_normalised_residual=metrics["max_local_t3_normalised_residual"],
        max_hidden_rank_mismatch=int(metrics["max_hidden_rank_mismatch"]),
        max_hidden_clock_residual=metrics["max_hidden_clock_residual"],
        max_hidden_loewner_violation=metrics["max_hidden_loewner_violation"],
        max_inverse_lambda_roundtrip_rel_error=metrics["max_inverse_lambda_roundtrip_rel_error"],
        max_inverse_x_roundtrip_rel_error=metrics["max_inverse_x_roundtrip_rel_error"],
        max_inverse_order_violation=metrics["max_inverse_order_violation"],
        max_transport_two_step_residual=metrics["max_transport_two_step_residual"],
        max_downstairs_associativity_residual=metrics["max_downstairs_associativity_residual"],
        max_long_chain_transport_residual=metrics["max_long_chain_transport_residual"],
        max_long_chain_clock_residual=metrics["max_long_chain_clock_residual"],
        max_local_hidden_birth_residual=metrics["max_local_hidden_birth_residual"],
        max_dv_visible_residual=metrics["max_dv_visible_residual"],
        max_dv_hidden_normalized_residual=metrics["max_dv_hidden_normalized_residual"],
        max_divergence_contraction_violation=metrics["max_divergence_contraction_violation"],
    )


def main() -> None:
    summary = run_validation_sweep()
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
