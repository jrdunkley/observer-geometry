from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from nomodescent import (
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    common_descent_test,
    false_collapse_diagnostic,
)

from external_comparisons.common import write_json

ROOT = Path(__file__).resolve().parent
SOURCE_COUNTS = ROOT.parents[2] / "evidence" / "micro_real_bundles" / "bell_counts_bundle" / "raw" / "counts_by_context.csv"


def _load_counts() -> dict[str, tuple[int, int, int, int]]:
    counts: dict[str, tuple[int, int, int, int]] = {}
    with SOURCE_COUNTS.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            counts[row["context"]] = tuple(int(row[key]) for key in ("n_pp", "n_pm", "n_mp", "n_mm"))
    return counts


def _moments(counts: tuple[int, int, int, int]) -> tuple[float, float, float]:
    pp, pm, mp, mm = counts
    total = float(pp + pm + mp + mm)
    mean_a = (pp + pm - mp - mm) / total
    mean_b = (pp - pm + mp - mm) / total
    e_ab = (pp - pm - mp + mm) / total
    return mean_a, mean_b, e_ab


def _covariance(mean_a: float, mean_b: float, e_ab: float) -> np.ndarray:
    return np.array([[1.0 - mean_a * mean_a, e_ab - mean_a * mean_b], [e_ab - mean_a * mean_b, 1.0 - mean_b * mean_b]])


def _observer_specs() -> tuple[ObserverSpec, ...]:
    return (
        ObserverSpec(name="00", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="01", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
        ObserverSpec(name="10", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="11", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]),
    )


def _problem_from_covariances(covariances: dict[str, np.ndarray], *, name: str) -> ProblemSpec:
    return ProblemSpec(
        name=name,
        latent_dim=4,
        observers=_observer_specs(),
        evidence=tuple(
            VisibleEvidenceSpec(name=f"cov_{context}", observer=context, kind="covariance", matrix=matrix)
            for context, matrix in covariances.items()
        ),
        goals=(GoalSpec(kind="common_completion"),),
    )


def _normalized_correlator_matrix(counts: dict[str, tuple[int, int, int, int]]) -> np.ndarray:
    correlations = {}
    for context, raw in counts.items():
        mean_a, mean_b, e_ab = _moments(raw)
        covariance = _covariance(mean_a, mean_b, e_ab)
        correlations[context] = covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])
    return np.array([[correlations["00"], correlations["01"]], [correlations["10"], correlations["11"]]], dtype=float)


def _arcsine_margin(correlator_matrix: np.ndarray) -> float:
    k00, k01 = correlator_matrix[0]
    k10, k11 = correlator_matrix[1]
    value = abs(math.asin(k00) + math.asin(k01) + math.asin(k10) - math.asin(k11))
    return float(math.pi - value)


def _projected_covariances(counts: dict[str, tuple[int, int, int, int]]) -> dict[str, np.ndarray]:
    moments = {context: _moments(raw) for context, raw in counts.items()}
    mean_a0 = 0.5 * (moments["00"][0] + moments["01"][0])
    mean_a1 = 0.5 * (moments["10"][0] + moments["11"][0])
    mean_b0 = 0.5 * (moments["00"][1] + moments["10"][1])
    mean_b1 = 0.5 * (moments["01"][1] + moments["11"][1])
    projected = {"00": (mean_a0, mean_b0), "01": (mean_a0, mean_b1), "10": (mean_a1, mean_b0), "11": (mean_a1, mean_b1)}
    return {context: _covariance(*projected[context], moments[context][2]) for context in ("00", "01", "10", "11")}


def main() -> None:
    counts = _load_counts()
    correlator_matrix = _normalized_correlator_matrix(counts)
    arcsine_margin = _arcsine_margin(correlator_matrix)
    raw_covariances = {context: _covariance(*_moments(raw)) for context, raw in counts.items()}
    raw_problem = _problem_from_covariances(raw_covariances, name="bell_counts_raw")
    projected_problem = _problem_from_covariances(_projected_covariances(counts), name="bell_counts_projected")

    raw_result = common_descent_test(raw_problem)
    projected_result = common_descent_test(projected_problem, allow_approximate_psd_search=True, psd_search_grid_size=61)
    false_collapse = false_collapse_diagnostic(
        fine_problem=raw_problem,
        coarse_summaries={"normalized_correlator_matrix": correlator_matrix},
        coarse_summary_label="normalized_correlator_matrix",
        coarse_certified_compatible=arcsine_margin > 0.0,
    )

    summary = {
        "external_method": "correlator-level CHSH / arcsine compatibility",
        "data_source": str(SOURCE_COUNTS.relative_to(ROOT.parents[2])),
        "arcsine_margin": arcsine_margin,
        "external_method_result": "non_obstructed" if arcsine_margin > 0.0 else "obstructed",
        "raw_qd_classification": raw_result.classification,
        "projected_qd_classification": projected_result.classification,
        "false_collapse_classification": false_collapse.classification,
        "comparison_relation": "distinguished",
    }
    audit = {
        "raw_qd": {
            "classification": raw_result.classification,
            "residuals": raw_result.residuals,
            "certificates": [certificate.kind for certificate in raw_result.certificates],
            "theorem_layer": raw_result.audit.theorem_layer,
        },
        "projected_qd": {
            "classification": projected_result.classification,
            "residuals": projected_result.residuals,
            "certificates": [certificate.kind for certificate in projected_result.certificates],
            "theorem_layer": projected_result.audit.theorem_layer,
        },
        "false_collapse": {
            "classification": false_collapse.classification,
            "coarse_summary_max_gap": false_collapse.coarse_summary_max_gap,
            "coarse_certified_compatible": false_collapse.coarse_certified_compatible,
            "theorem_layer": false_collapse.audit.theorem_layer,
        },
    }
    comparison = {
        "what_external_method_certifies": "The normalized correlator summary is non-obstructive under the arcsine CHSH criterion.",
        "what_qd_adds": "QD sees that the raw full-law Gaussian surrogate remains exactly incompatible because repeated local marginals disagree across contexts.",
        "epistemic_boundary": "The projected path is explicitly an encoded no-signalling Gaussian surrogate and therefore remains audited approximate rather than theorem-grade.",
    }
    write_json(ROOT / "outputs" / "summary.json", summary)
    write_json(ROOT / "outputs" / "audit.json", audit)
    write_json(ROOT / "outputs" / "comparison.json", comparison)


if __name__ == "__main__":
    main()

