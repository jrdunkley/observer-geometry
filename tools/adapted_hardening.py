from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from nomogeo import closure_adapted_observer, whitened_perturbation
from nomogeo.exceptions import InputValidationError

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tools" / "outputs"


@dataclass(frozen=True)
class AdaptedHardeningSummary:
    seed: int
    commuting_trials: int
    bridge_trials: int
    max_commuting_eta: float
    max_commuting_visible_score_gap: float
    max_bridge_eta: float
    max_bridge_scaled_quartic_residual: float
    rotated_condition_acceptance: dict[str, int]
    real_bundle_commutator_norms: dict[str, float]
    real_bundle_noncommuting_rejections: dict[str, bool]
    real_bundle_exact_success_eta: dict[str, float]
    real_bundle_exact_success_score: dict[str, float]


def _sqrt_spd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T


def _orthogonal(rng: np.random.Generator, n: int) -> np.ndarray:
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return q


def _random_spd(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a + n * np.eye(n, dtype=float)


def _fro(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, ord="fro"))


def _read_numeric_csv(path: Path, skip_first: bool = True) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        labels: list[str] = []
        rows: list[list[float]] = []
        for row in reader:
            labels.append(row[0])
            values = [float(x) for x in row[1 if skip_first else 0 :]]
            rows.append(values)
    return labels, np.asarray(rows, dtype=float)


def _empirical_reference(data: np.ndarray, ridge: float = 1e-1) -> np.ndarray:
    centered = data - np.mean(data, axis=0, keepdims=True)
    covariance = np.cov(centered, rowvar=False, bias=False)
    covariance = 0.5 * (covariance + covariance.T)
    return covariance + ridge * np.eye(covariance.shape[0], dtype=float)


def _empirical_noncommuting_family_iris() -> tuple[np.ndarray, list[np.ndarray]]:
    labels, data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "iris_protocol_mismatch" / "raw" / "iris_slice.csv")
    H = _empirical_reference(data)
    family = []
    unique = sorted(set(labels))
    for label in unique:
        subset = data[[idx for idx, item in enumerate(labels) if item == label]]
        centered = subset - np.mean(subset, axis=0, keepdims=True)
        cov = np.cov(centered, rowvar=False, bias=False)
        family.append(0.5 * (cov + cov.T) - H)
    return H, family


def _empirical_noncommuting_family_leaderboard() -> tuple[np.ndarray, list[np.ndarray]]:
    _labels, data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "leaderboard_benchmark_slice" / "raw" / "score_slice.csv")
    H = _empirical_reference(data)
    centered = data - np.mean(data, axis=0, keepdims=True)
    family = []
    for row in centered:
        family.append(np.outer(row, row) - H)
    return H, family


def _empirical_noncommuting_family_bell() -> tuple[np.ndarray, list[np.ndarray]]:
    path = ROOT / "evidence" / "micro_real_bundles" / "bell_counts_bundle" / "raw" / "counts_by_context.csv"
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                [
                    float(row["n_pp"]),
                    float(row["n_pm"]),
                    float(row["n_mp"]),
                    float(row["n_mm"]),
                ]
            )
    counts = np.asarray(rows, dtype=float)
    probs = counts / np.sum(counts, axis=1, keepdims=True)
    H = _empirical_reference(probs, ridge=5e-2)
    centered = probs - np.mean(probs, axis=0, keepdims=True)
    family = [np.outer(row, row) - H for row in centered]
    return H, family


def _empirical_exact_family_iris() -> tuple[np.ndarray, list[np.ndarray]]:
    labels, data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "iris_protocol_mismatch" / "raw" / "iris_slice.csv")
    H = _empirical_reference(data)
    H_half = _sqrt_spd(H)
    U = np.linalg.eigh(H)[1]
    overall = np.mean(data, axis=0)
    rows = []
    for label in sorted(set(labels)):
        subset = data[[idx for idx, item in enumerate(labels) if item == label]]
        delta = np.mean(subset, axis=0) - overall
        rows.append(U.T @ delta)
    family = [H_half @ U @ np.diag(row) @ U.T @ H_half for row in rows]
    return H, family


def _empirical_exact_family_leaderboard() -> tuple[np.ndarray, list[np.ndarray]]:
    _labels, data = _read_numeric_csv(ROOT / "evidence" / "micro_real_bundles" / "leaderboard_benchmark_slice" / "raw" / "score_slice.csv")
    H = _empirical_reference(data)
    H_half = _sqrt_spd(H)
    U = np.linalg.eigh(H)[1]
    centered = data - np.mean(data, axis=0, keepdims=True)
    family = [H_half @ U @ np.diag(U.T @ row) @ U.T @ H_half for row in centered]
    return H, family


def _max_commutator_norm(H: np.ndarray, family: list[np.ndarray]) -> float:
    whitened = [whitened_perturbation(H, Delta) for Delta in family]
    best = 0.0
    for i, left in enumerate(whitened):
        for right in whitened[i + 1 :]:
            best = max(best, _fro(left @ right - right @ left))
    return best


def run_adapted_hardening(seed: int = 20260409, commuting_trials: int = 40, bridge_trials: int = 20) -> AdaptedHardeningSummary:
    rng = np.random.default_rng(seed)

    max_commuting_eta = 0.0
    max_commuting_visible_score_gap = 0.0
    for _ in range(commuting_trials):
        n = int(rng.integers(5, 9))
        m = int(rng.integers(1, min(4, n) + 1))
        task_count = int(rng.integers(2, 5))
        H = _random_spd(rng, n)
        H_half = _sqrt_spd(H)
        U = _orthogonal(rng, n)
        lambdas = rng.normal(size=(task_count, n))
        lambdas[:, m + 1 :] = 0.0
        family = [H_half @ U @ np.diag(row) @ U.T @ H_half for row in lambdas]
        result = closure_adapted_observer(H, family, m)
        mu = np.sum(lambdas**2, axis=0)
        expected_score = float(np.sum(np.sort(mu)[::-1][:m]))
        max_commuting_eta = max(max_commuting_eta, result.scores.eta)
        max_commuting_visible_score_gap = max(max_commuting_visible_score_gap, abs(result.scores.visible_score - expected_score))

    max_bridge_eta = 0.0
    max_bridge_scaled_quartic_residual = 0.0
    from nomogeo import dv_bridge, local_visible_calculus, visible_precision  # local import keeps the script narrow

    for _ in range(bridge_trials):
        H0 = _random_spd(rng, 7)
        J = rng.normal(size=(7, 7))
        J = J - J.T
        bridge = dv_bridge(H0, J)
        result = closure_adapted_observer(H0, [bridge.delta_dv], 2)
        local = local_visible_calculus(H0, result.C, bridge.delta_dv)
        eps = 0.2
        phi_eps = visible_precision(H0 + (eps * eps) * bridge.delta_dv, result.C)
        residual = _fro(phi_eps - np.eye(2) - (eps * eps) * local.V) / (eps**4)
        max_bridge_eta = max(max_bridge_eta, result.scores.eta)
        max_bridge_scaled_quartic_residual = max(max_bridge_scaled_quartic_residual, residual)

    rotated_condition_acceptance: dict[str, int] = {}
    for exponent in (4, 6, 8, 9):
        accepted = 0
        for _ in range(12):
            n = 7
            m = 3
            vals = np.geomspace(1.0, 10.0**exponent, n)
            Q = _orthogonal(rng, n)
            H = Q @ np.diag(vals) @ Q.T
            U = _orthogonal(rng, n)
            lambdas = np.array(
                [
                    [5.0, 4.99999999, 4.1, 1.0, 0.0, 0.0, 0.0],
                    [4.2, 4.19999999, 3.9, 0.5, 0.0, 0.0, 0.0],
                ],
                dtype=float,
            )
            H_half = Q @ np.diag(np.sqrt(vals)) @ Q.T
            family = [H_half @ U @ np.diag(row) @ U.T @ H_half for row in lambdas]
            try:
                closure_adapted_observer(H, family, m)
                accepted += 1
            except InputValidationError:
                pass
        rotated_condition_acceptance[f"cond_1e{exponent}"] = accepted

    real_bundle_commutator_norms: dict[str, float] = {}
    real_bundle_noncommuting_rejections: dict[str, bool] = {}
    for name, builder, rank in (
        ("iris_covariance_family", _empirical_noncommuting_family_iris, 2),
        ("leaderboard_outerproduct_family", _empirical_noncommuting_family_leaderboard, 2),
        ("bell_probability_family", _empirical_noncommuting_family_bell, 2),
    ):
        H, family = builder()
        real_bundle_commutator_norms[name] = _max_commutator_norm(H, family)
        try:
            closure_adapted_observer(H, family, rank)
            real_bundle_noncommuting_rejections[name] = False
        except InputValidationError:
            real_bundle_noncommuting_rejections[name] = True

    real_bundle_exact_success_eta: dict[str, float] = {}
    real_bundle_exact_success_score: dict[str, float] = {}
    for name, builder, rank in (
        ("iris_mean_spectral_family", _empirical_exact_family_iris, 2),
        ("leaderboard_model_spectral_family", _empirical_exact_family_leaderboard, 2),
    ):
        H, family = builder()
        result = closure_adapted_observer(H, family, rank)
        real_bundle_exact_success_eta[name] = result.scores.eta
        real_bundle_exact_success_score[name] = result.scores.visible_score

    return AdaptedHardeningSummary(
        seed=seed,
        commuting_trials=commuting_trials,
        bridge_trials=bridge_trials,
        max_commuting_eta=max_commuting_eta,
        max_commuting_visible_score_gap=max_commuting_visible_score_gap,
        max_bridge_eta=max_bridge_eta,
        max_bridge_scaled_quartic_residual=max_bridge_scaled_quartic_residual,
        rotated_condition_acceptance=rotated_condition_acceptance,
        real_bundle_commutator_norms=real_bundle_commutator_norms,
        real_bundle_noncommuting_rejections=real_bundle_noncommuting_rejections,
        real_bundle_exact_success_eta=real_bundle_exact_success_eta,
        real_bundle_exact_success_score=real_bundle_exact_success_score,
    )


def write_adapted_hardening_outputs(summary: AdaptedHardeningSummary) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "adapted_hardening_summary.json"
    path.write_text(json.dumps(asdict(summary), indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    summary = run_adapted_hardening()
    write_adapted_hardening_outputs(summary)
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
