from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from nomogeo import hidden_load, visible_precision

from examples.common import bar_chart_svg, line_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def latent_precision() -> np.ndarray:
    return np.array(
        [
            [4.0, 1.8, 1.4],
            [1.8, 3.0, 0.1],
            [1.4, 0.1, 2.2],
        ],
        dtype=float,
    )


def observers() -> dict[str, np.ndarray]:
    return {
        "plurality": np.array([[1.0, 0.0, 0.0]]),
        "approval": np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        "ranked_choice": np.eye(3, dtype=float),
    }


def observer_ceiling(name: str, h: np.ndarray) -> np.ndarray:
    if name == "plurality":
        return h[:1, :1]
    if name == "approval":
        return h[:2, :2]
    if name == "ranked_choice":
        return h
    raise KeyError(name)


def feature_labels(name: str) -> list[str]:
    if name == "plurality":
        return ["first_place"]
    if name == "approval":
        return ["first_place", "approval_support"]
    return ["first_place", "approval_support", "transfer_structure"]


def _write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h = latent_precision()

    summary_rows: list[dict[str, float | int | str]] = []
    spectra_rows: list[dict[str, float | int | str]] = []
    observer_rows: list[dict[str, str]] = []
    clocks: list[float] = []
    names: list[str] = []

    for name, observer in observers().items():
        observer_rows.append({"observer": name, "matrix": np.array2string(observer, separator=", ")})

    for name, observer in observers().items():
        phi = visible_precision(h, observer)
        ceiling = observer_ceiling(name, h)
        load = hidden_load(ceiling, phi, support_mode="ambient")
        eigs = np.linalg.eigvalsh(load.reduced_lambda) if load.reduced_lambda.size else np.array([0.0])
        logdet_gap = float(np.linalg.slogdet(ceiling)[1] - np.linalg.slogdet(phi)[1])
        summary_rows.append(
            {
                "observer": name,
                "visible_dim": observer.shape[0],
                "latent_dim": h.shape[0],
                "observer_rank_deficiency": h.shape[0] - observer.shape[0],
                "clock": float(load.clock),
                "hidden_rank": int(load.rank),
                "direct_logdet_gap": logdet_gap,
                "clock_minus_direct_gap_abs": abs(float(load.clock - logdet_gap)),
                "preserved_features": ",".join(feature_labels(name)),
            }
        )
        for idx, eig in enumerate(eigs, start=1):
            spectra_rows.append({"observer": name, "eigen_index": idx, "eigenvalue": float(eig)})
        names.append(name)
        clocks.append(float(load.clock))

    _write_csv(OUT / "observer_summary.csv", summary_rows)
    _write_csv(OUT / "hidden_spectra.csv", spectra_rows)
    _write_csv(OUT / "observer_maps.csv", observer_rows)

    spectra_by_name: dict[str, np.ndarray] = {}
    for name in names:
        values = [row["eigenvalue"] for row in spectra_rows if row["observer"] == name]
        padded = np.zeros(3, dtype=float)
        padded[: len(values)] = np.array(values, dtype=float)
        spectra_by_name[name] = padded

    bar_chart_svg(
        OUT / "clock_comparison.svg",
        names,
        np.array(clocks, dtype=float),
        title="Clock By Voting Observer",
        y_label="clock",
        color="#1d3557",
    )
    line_chart_svg(
        OUT / "hidden_spectra.svg",
        np.array([1.0, 2.0, 3.0], dtype=float),
        [(name, spectra_by_name[name], color) for name, color in zip(names, ["#c1121f", "#2a9d8f", "#577590"], strict=True)],
        title="Hidden-Load Spectra Across Voting Observers",
        x_label="eigenvalue index",
        y_label="eigenvalue",
    )

    phi_ranked = visible_precision(h, observers()["ranked_choice"])
    phi_approval = visible_precision(h, observers()["approval"])
    phi_plurality = visible_precision(h, observers()["plurality"])
    step_ranked_to_approval = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    step_approval_to_plurality = np.array([[1.0, 0.0]])
    composed_phi = visible_precision(visible_precision(phi_ranked, step_ranked_to_approval), step_approval_to_plurality)

    summary = {
        "claim": (
            "In this synthetic latent preference geometry, plurality, approval, and ranked choice act as nested observers "
            "with different hidden-load spectra and determinant-gap clocks."
        ),
        "clock_ordering_plurality_ge_approval_ge_ranked": bool(clocks[0] >= clocks[1] >= clocks[2]),
        "composition_residual_abs": float(np.linalg.norm(composed_phi - phi_plurality, ord=np.inf)),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
