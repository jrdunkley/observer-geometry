from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.fft import dct

from nomogeo import visible_precision
from nomodescent import ObserverSpec, classify_relation, staged_descent_check

from worked_examples.common import bar_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def free_gaussian_precision(n: int = 8, mass: float = 1.3, coupling: float = 0.7) -> np.ndarray:
    lap = 2.0 * np.eye(n)
    for i in range(n - 1):
        lap[i, i + 1] = -1.0
        lap[i + 1, i] = -1.0
    return mass * np.eye(n) + coupling * lap


def dct_basis(n: int) -> np.ndarray:
    return dct(np.eye(n), norm="ortho", axis=0)


def observers() -> tuple[ObserverSpec, ObserverSpec, ObserverSpec]:
    basis = dct_basis(8)
    mode4 = ObserverSpec(name="mode4", matrix=basis[:4, :])
    mode2 = ObserverSpec(name="mode2", matrix=basis[:2, :])
    block4 = ObserverSpec(
        name="block4",
        matrix=(
            np.array(
                [
                    [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
                ],
                dtype=float,
            )
            / np.sqrt(2.0)
        ),
    )
    return mode4, mode2, block4


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    H = free_gaussian_precision()
    mode4, mode2, block4 = observers()

    tower = staged_descent_check(H, [mode4, ObserverSpec(name="keep2of4", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])])
    relation_mode = classify_relation(mode2, mode4)
    relation_block = classify_relation(block4, mode4)

    phi_mode4 = visible_precision(H, mode4.matrix)
    step = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    phi_mode2 = visible_precision(phi_mode4, step)
    phi_block4 = visible_precision(H, block4.matrix)

    rows = [
        {"relation": "mode2_through_mode4", "classification": relation_mode.classification, "residual": relation_mode.residuals["a_through_b_residual"]},
        {"relation": "block4_vs_mode4", "classification": relation_block.classification, "residual": relation_block.residuals["a_through_b_residual"]},
        {"relation": "tower_mode4_to_mode2", "classification": tower.classification, "residual": tower.residuals["tower_residual"]},
        {"relation": "mode4_to_mode2_logdet_gap", "classification": "exact_visible_gap", "residual": float(np.linalg.slogdet(phi_mode4)[1] - np.linalg.slogdet(phi_mode2)[1])},
    ]
    _write_csv(OUT / "relation_summary.csv", rows)

    bar_chart_svg(
        OUT / "visible_logdet_comparison.svg",
        ["mode4", "mode2", "block4"],
        np.array(
            [
                np.linalg.slogdet(phi_mode4)[1],
                np.linalg.slogdet(phi_mode2)[1],
                np.linalg.slogdet(phi_block4)[1],
            ]
        ),
        title="Free Gaussian Quotient Descent",
        y_label="log det visible precision",
        color="#355070",
    )

    summary = {
        "tower_residual": tower.residuals["tower_residual"],
        "mode_relation": relation_mode.classification,
        "block_relation": relation_block.classification,
        "mode4_to_mode2_logdet_gap": float(np.linalg.slogdet(phi_mode4)[1] - np.linalg.slogdet(phi_mode2)[1]),
        "boundary": "interacting RG is deliberately deferred; this example is exact only in the free Gaussian regime",
    }
    audit = {
        "theorem_layer": "exact quotient tower law and observer relation tests in the free Gaussian regime",
        "residuals": {
            "tower_residual": tower.residuals["tower_residual"],
            "mode_factorisation_residual": relation_mode.residuals["a_through_b_residual"],
            "block_factorisation_residual": relation_block.residuals["a_through_b_residual"],
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
