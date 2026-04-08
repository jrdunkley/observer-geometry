from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from nomogeo import visible_precision

from examples.common import heatmap_svg, line_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
SIGNS = np.array([[1.0, 1.0], [1.0, -1.0]])
CLASS_LABELS = {
    0: "compatible",
    1: "variance-only obstruction",
    2: "correlator-only obstruction",
    3: "both obstructions",
}


def pair_selectors() -> dict[tuple[int, int], np.ndarray]:
    return {
        (0, 0): np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        (0, 1): np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
        (1, 0): np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        (1, 1): np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    }


def family_covariances(delta: float, rho: float) -> dict[tuple[int, int], np.ndarray]:
    variances_a = {(0, 0): 1.0, (0, 1): 1.0 + delta, (1, 0): 1.0, (1, 1): 1.0}
    variances_b = {(0, 0): 1.0, (0, 1): 1.0, (1, 0): 1.0, (1, 1): 1.0}
    family: dict[tuple[int, int], np.ndarray] = {}
    for x in range(2):
        for y in range(2):
            va = variances_a[(x, y)]
            vb = variances_b[(x, y)]
            covariance = SIGNS[x, y] * rho * np.sqrt(va * vb)
            family[(x, y)] = np.array([[va, covariance], [covariance, vb]], dtype=float)
    return family


def normalized_correlators(family: dict[tuple[int, int], np.ndarray]) -> np.ndarray:
    matrix = np.zeros((2, 2), dtype=float)
    for x in range(2):
        for y in range(2):
            cov = family[(x, y)]
            matrix[x, y] = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    return matrix


def arcsine_chsh_value(k: np.ndarray) -> float:
    return float(abs(np.arcsin(k[0, 0]) + np.arcsin(k[0, 1]) + np.arcsin(k[1, 0]) - np.arcsin(k[1, 1])))


def arcsine_margin(family: dict[tuple[int, int], np.ndarray]) -> float:
    return float(np.pi - arcsine_chsh_value(normalized_correlators(family)))


def variance_gap(family: dict[tuple[int, int], np.ndarray]) -> float:
    gaps = [
        abs(float(family[(0, 0)][0, 0] - family[(0, 1)][0, 0])),
        abs(float(family[(1, 0)][0, 0] - family[(1, 1)][0, 0])),
        abs(float(family[(0, 0)][1, 1] - family[(1, 0)][1, 1])),
        abs(float(family[(0, 1)][1, 1] - family[(1, 1)][1, 1])),
    ]
    return max(gaps)


def classify_family(delta: float, rho: float) -> int:
    family = family_covariances(delta, rho)
    has_variance_obstruction = variance_gap(family) > 1e-12
    has_correlator_obstruction = arcsine_margin(family) < -1e-12
    if not has_variance_obstruction and not has_correlator_obstruction:
        return 0
    if has_variance_obstruction and not has_correlator_obstruction:
        return 1
    if not has_variance_obstruction and has_correlator_obstruction:
        return 2
    return 3


def linear_completion_residual(family: dict[tuple[int, int], np.ndarray]) -> float:
    selectors = pair_selectors()
    index_pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    matrix = np.zeros((16, len(index_pairs)), dtype=float)
    rhs = np.zeros(16, dtype=float)
    row = 0
    for key in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        selector = selectors[key]
        target = family[key]
        for a in range(2):
            for b in range(2):
                left = selector[a]
                right = selector[b]
                for col, (i, j) in enumerate(index_pairs):
                    coefficient = left[i] * right[j]
                    if i != j:
                        coefficient += left[j] * right[i]
                    matrix[row, col] = coefficient
                rhs[row] = target[a, b]
                row += 1
    solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    residual = matrix @ solution - rhs
    return float(np.linalg.norm(residual, ord=np.inf))


def consistent_completion_matrix(rho: float, alpha: float, beta: float) -> np.ndarray:
    return np.array(
        [
            [1.0, alpha, rho, rho],
            [alpha, 1.0, rho, -rho],
            [rho, rho, 1.0, beta],
            [rho, -rho, beta, 1.0],
        ],
        dtype=float,
    )


def best_psd_completion(rho: float, resolution: int = 241) -> tuple[np.ndarray, float, float, float]:
    alphas = np.linspace(-1.0, 1.0, resolution)
    betas = np.linspace(-1.0, 1.0, resolution)
    best_gap = -np.inf
    best_alpha = 0.0
    best_beta = 0.0
    best_sigma = consistent_completion_matrix(rho, 0.0, 0.0)
    for alpha in alphas:
        for beta in betas:
            sigma = consistent_completion_matrix(rho, float(alpha), float(beta))
            gap = float(np.min(np.linalg.eigvalsh(sigma)))
            if gap > best_gap:
                best_gap = gap
                best_alpha = float(alpha)
                best_beta = float(beta)
                best_sigma = sigma
    return best_sigma, best_alpha, best_beta, best_gap


def pair_precision_residuals(sigma: np.ndarray) -> dict[str, float]:
    selectors = pair_selectors()
    family = family_covariances(delta=0.0, rho=float(sigma[0, 2]))
    residuals: dict[str, float] = {}
    h = np.linalg.inv(sigma)
    for key, selector in selectors.items():
        observed_precision = visible_precision(h, selector)
        target_precision = np.linalg.inv(family[key])
        residuals[f"{key[0]}{key[1]}"] = float(np.linalg.norm(observed_precision - target_precision, ord=np.inf))
    return residuals


def _write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    deltas = np.linspace(0.0, 0.35, 71)
    rhos = np.linspace(0.0, 0.95, 96)
    phase_rows: list[dict[str, float | int | str]] = []
    phase_grid = np.zeros((len(deltas), len(rhos)), dtype=int)

    for iy, delta in enumerate(deltas):
        for ix, rho in enumerate(rhos):
            family = family_covariances(float(delta), float(rho))
            cls = classify_family(float(delta), float(rho))
            phase_grid[iy, ix] = cls
            phase_rows.append(
                {
                    "delta": float(delta),
                    "rho": float(rho),
                    "class_id": cls,
                    "class_label": CLASS_LABELS[cls],
                    "variance_gap": variance_gap(family),
                    "arcsine_chsh_margin": arcsine_margin(family),
                    "linear_completion_residual": linear_completion_residual(family),
                }
            )

    _write_csv(OUT / "phase_diagram.csv", phase_rows)
    heatmap_svg(
        OUT / "phase_diagram.svg",
        rhos,
        deltas,
        phase_grid,
        palette={0: "#d8f3dc", 1: "#f9c74f", 2: "#f9844a", 3: "#c1121f"},
        labels=CLASS_LABELS,
        title="Bell Square: Correlator vs Full-Law Compatibility",
        x_label="rho",
        y_label="delta",
    )

    delta_slice = np.linspace(0.0, 0.35, 36)
    rho_slice = 0.6
    slice_margin = np.array([arcsine_margin(family_covariances(float(delta), rho_slice)) for delta in delta_slice])
    slice_var = np.array([variance_gap(family_covariances(float(delta), rho_slice)) for delta in delta_slice])
    line_chart_svg(
        OUT / "variance_slice.svg",
        delta_slice,
        [
            ("variance gap", slice_var, "#bc6c25"),
            ("arcsine margin", slice_margin, "#355070"),
        ],
        title="Open Separation Slice at rho = 0.6",
        x_label="delta",
        y_label="gap / margin",
    )

    compatible_sigma, alpha, beta, best_gap = best_psd_completion(rho=0.6)
    pair_residuals = pair_precision_residuals(compatible_sigma)
    sample_rows = [
        {
            "sample": "compatible",
            "delta": 0.0,
            "rho": 0.6,
            "class_label": CLASS_LABELS[classify_family(0.0, 0.6)],
            "correlator_matches_compatible": True,
            "variance_gap": variance_gap(family_covariances(0.0, 0.6)),
            "arcsine_chsh_margin": arcsine_margin(family_covariances(0.0, 0.6)),
            "linear_completion_residual": linear_completion_residual(family_covariances(0.0, 0.6)),
            "best_psd_gap": best_gap,
        },
        {
            "sample": "variance_only",
            "delta": 0.22,
            "rho": 0.6,
            "class_label": CLASS_LABELS[classify_family(0.22, 0.6)],
            "correlator_matches_compatible": True,
            "variance_gap": variance_gap(family_covariances(0.22, 0.6)),
            "arcsine_chsh_margin": arcsine_margin(family_covariances(0.22, 0.6)),
            "linear_completion_residual": linear_completion_residual(family_covariances(0.22, 0.6)),
            "best_psd_gap": float("nan"),
        },
        {
            "sample": "correlator_only",
            "delta": 0.0,
            "rho": 0.82,
            "class_label": CLASS_LABELS[classify_family(0.0, 0.82)],
            "correlator_matches_compatible": False,
            "variance_gap": variance_gap(family_covariances(0.0, 0.82)),
            "arcsine_chsh_margin": arcsine_margin(family_covariances(0.0, 0.82)),
            "linear_completion_residual": linear_completion_residual(family_covariances(0.0, 0.82)),
            "best_psd_gap": best_psd_completion(rho=0.82)[3],
        },
    ]
    _write_csv(OUT / "sample_points.csv", sample_rows)
    _write_csv(OUT / "comparison_table.csv", sample_rows)

    summary = {
        "claim": (
            "There exists an explicit Gaussian family where the normalized correlator summary is unchanged, "
            "the correlator-level Bell obstruction is absent, and common gluing still fails at the full law level."
        ),
        "parameter_ranges": {
            "delta": [float(deltas[0]), float(deltas[-1])],
            "rho": [float(rhos[0]), float(rhos[-1])],
        },
        "compatible_fraction": float(np.mean(phase_grid == 0)),
        "variance_only_fraction": float(np.mean(phase_grid == 1)),
        "correlator_only_fraction": float(np.mean(phase_grid == 2)),
        "both_fraction": float(np.mean(phase_grid == 3)),
        "compatible_completion_alpha": alpha,
        "compatible_completion_beta": beta,
        "compatible_completion_gap": best_gap,
        "compatible_pair_precision_residual_max": max(pair_residuals.values()),
        "compatible_pair_precision_residuals": pair_residuals,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
