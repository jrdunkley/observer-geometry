from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from nomodescent import common_descent_test

from evidence import (
    ExtractionNote,
    ExtractionRecord,
    MatrixExcerpt,
    ObserverHypothesis,
    ProtocolExcerpt,
    QuotedTextExcerpt,
    SourceRef,
    assemble_problem_spec,
    build_evidence_bundle_from_excerpts,
)

from micro_real_bundles.common import assembly_result_to_dict, descent_result_to_dict, ensure_dir, write_json, write_rows_csv

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
EXTRACTED = ROOT / "extracted"


def _load_counts() -> dict[str, dict[str, object]]:
    counts: dict[str, dict[str, object]] = {}
    with (RAW / "counts_by_context.csv").open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            counts[row["context"]] = {
                "alpha": row["alpha"],
                "beta": row["beta"],
                "counts": tuple(int(row[key]) for key in ("n_pp", "n_pm", "n_mp", "n_mm")),
            }
    return counts


def _pair_probabilities(counts: tuple[int, int, int, int]) -> np.ndarray:
    return np.asarray(counts, dtype=float) / float(sum(counts))


def _moments(probabilities: np.ndarray) -> tuple[float, float, float]:
    values = np.array([[1.0, 1.0], [1.0, -1.0], [-1.0, 1.0], [-1.0, -1.0]], dtype=float)
    mean = (probabilities[:, None] * values).sum(axis=0)
    e_ab = float(np.dot(probabilities, np.array([1.0, -1.0, -1.0, 1.0], dtype=float)))
    return float(mean[0]), float(mean[1]), e_ab


def _covariance_from_moments(mean_a: float, mean_b: float, e_ab: float) -> np.ndarray:
    return np.array(
        [
            [1.0 - mean_a * mean_a, e_ab - mean_a * mean_b],
            [e_ab - mean_a * mean_b, 1.0 - mean_b * mean_b],
        ],
        dtype=float,
    )


def _observer_hypotheses() -> tuple[ObserverHypothesis, ...]:
    source = SourceRef(source="bell_counts_bundle", location="context selector table")
    extraction = ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True)
    return (
        ObserverHypothesis(name="00", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]], protocol_name="00", features=("A0", "B0"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="01", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]], protocol_name="01", features=("A0", "B1"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="10", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]], protocol_name="10", features=("A1", "B0"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="11", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]], protocol_name="11", features=("A1", "B1"), source_ref=source, extraction=extraction),
    )


def _protocol_excerpts(counts: dict[str, dict[str, object]]) -> tuple[ProtocolExcerpt, ...]:
    extraction = ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True)
    return tuple(
        ProtocolExcerpt(
            name=context,
            facts=(f"A{context[0]}", f"B{context[1]}"),
            candidate_family=(context,),
            source_ref=SourceRef(source="SciPostPhys.10.1.017", location=f"Table A2 context {context}"),
            extraction=extraction,
            selection_provenance="manual context extraction from Table A2",
            notes=f"alpha={payload['alpha']}, beta={payload['beta']}",
        )
        for context, payload in counts.items()
    )


def _counts_excerpt(counts: dict[str, dict[str, object]]) -> QuotedTextExcerpt:
    lines = [
        f"{context}: alpha={payload['alpha']}, beta={payload['beta']}, counts={payload['counts']}"
        for context, payload in counts.items()
    ]
    return QuotedTextExcerpt(
        name="bell_count_summary",
        quote=" | ".join(lines),
        source_ref=SourceRef(source="SciPostPhys.10.1.017", location="Appendix C / Table A2"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        selection_provenance="manual Bell count table transcription",
    )


def _raw_covariance_excerpts(counts: dict[str, dict[str, object]]) -> tuple[MatrixExcerpt, ...]:
    excerpts: list[MatrixExcerpt] = []
    for context, payload in counts.items():
        probabilities = _pair_probabilities(payload["counts"])
        mean_a, mean_b, e_ab = _moments(probabilities)
        covariance = _covariance_from_moments(mean_a, mean_b, e_ab)
        excerpts.append(
            MatrixExcerpt(
                name=f"cov_{context}",
                matrix_role="visible_object",
                matrix_kind="covariance",
                matrix=tuple(tuple(float(entry) for entry in row) for row in covariance),
                observer_name=context,
                source_ref=SourceRef(source="SciPostPhys.10.1.017", location=f"derived from counts {context}"),
                extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
                selection_provenance="exact covariance of +/-1 outcomes derived from explicit counts",
                notes="binary-law covariance used as Gaussian surrogate visible object",
            )
        )
    return tuple(excerpts)


def _projected_covariance_excerpts(counts: dict[str, dict[str, object]]) -> tuple[MatrixExcerpt, ...]:
    raw_moments = {context: _moments(_pair_probabilities(payload["counts"])) for context, payload in counts.items()}
    mean_a0 = 0.5 * (raw_moments["00"][0] + raw_moments["01"][0])
    mean_a1 = 0.5 * (raw_moments["10"][0] + raw_moments["11"][0])
    mean_b0 = 0.5 * (raw_moments["00"][1] + raw_moments["10"][1])
    mean_b1 = 0.5 * (raw_moments["01"][1] + raw_moments["11"][1])
    projected_means = {"00": (mean_a0, mean_b0), "01": (mean_a0, mean_b1), "10": (mean_a1, mean_b0), "11": (mean_a1, mean_b1)}

    excerpts: list[MatrixExcerpt] = []
    for context, (_raw_a, _raw_b, e_ab) in raw_moments.items():
        mean_a, mean_b = projected_means[context]
        covariance = _covariance_from_moments(mean_a, mean_b, e_ab)
        excerpts.append(
            MatrixExcerpt(
                name=f"cov_{context}",
                matrix_role="visible_object",
                matrix_kind="covariance",
                matrix=tuple(tuple(float(entry) for entry in row) for row in covariance),
                observer_name=context,
                source_ref=SourceRef(source="bell_counts_bundle", location=f"projected covariance {context}"),
                extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
                selection_provenance="no-signalling marginal projection after exact count extraction",
                notes="local means averaged across matching settings before Gaussian common-descent test",
            )
        )
    return tuple(excerpts)


def _bundle_notes(projected: bool) -> tuple[ExtractionNote, ...]:
    notes = [
        ExtractionNote(
            name="gaussian_surrogate",
            statement="Binary pair covariances are treated as Gaussian surrogate visible objects downstream.",
            load_bearing=True,
            source_ref=SourceRef(source="bell_counts_bundle", location="bundle design"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred"),
        )
    ]
    if projected:
        notes.append(
            ExtractionNote(
                name="no_signalling_projection",
                statement="Repeated local means are averaged across matching settings to remove finite-sample signalling before Gaussian common-descent testing.",
                load_bearing=True,
                source_ref=SourceRef(source="bell_counts_bundle", location="projection layer"),
                extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred"),
            )
        )
    return tuple(notes)


def build_bundles() -> tuple[object, object]:
    counts = _load_counts()
    base_excerpts = (_counts_excerpt(counts), *_protocol_excerpts(counts))
    raw_bundle = build_evidence_bundle_from_excerpts(
        name="bell_counts_raw_gaussian_surrogate",
        latent_dim=4,
        excerpts=(*base_excerpts, *_raw_covariance_excerpts(counts)),
        observer_hypotheses=_observer_hypotheses(),
        notes=_bundle_notes(projected=False),
        description="Exact Bell counts with raw finite-sample Gaussian-surrogate pair covariances",
        provenance="SciPost Phys. 10, 017 (2021), Table A2 count transcription",
        tags=("micro-real", "bell"),
    )
    projected_bundle = build_evidence_bundle_from_excerpts(
        name="bell_counts_projected_gaussian_surrogate",
        latent_dim=4,
        excerpts=(*base_excerpts, *_projected_covariance_excerpts(counts)),
        observer_hypotheses=_observer_hypotheses(),
        notes=_bundle_notes(projected=True),
        description="Bell counts with explicit no-signalling marginal projection before Gaussian-surrogate descent",
        provenance="SciPost Phys. 10, 017 (2021), Table A2 transcription plus encoded no-signalling projection",
        tags=("micro-real", "bell"),
    )
    return raw_bundle, projected_bundle


def main() -> None:
    ensure_dir(EXTRACTED)
    counts = _load_counts()
    raw_bundle, projected_bundle = build_bundles()
    raw_bundle.bundle.to_json(ROOT / "bundle_raw.json")
    projected_bundle.bundle.to_json(ROOT / "bundle.json")

    raw_assembly = assemble_problem_spec(raw_bundle.bundle)
    projected_assembly = assemble_problem_spec(projected_bundle.bundle)
    if raw_assembly.problem_spec is None or projected_assembly.problem_spec is None:
        raise RuntimeError("Bell bundle assemblies should succeed")

    raw_result = common_descent_test(raw_assembly.problem_spec, allow_approximate_psd_search=True, psd_search_grid_size=31)
    projected_result = common_descent_test(projected_assembly.problem_spec, allow_approximate_psd_search=True, psd_search_grid_size=61)

    raw_assembly.problem_spec.to_json(ROOT / "problem_spec_raw.json")
    projected_assembly.problem_spec.to_json(ROOT / "problem_spec.json")

    probability_rows: list[dict[str, object]] = []
    raw_covariances: dict[str, object] = {}
    projected_covariances: dict[str, object] = {}
    raw_matrix_map = {item.name: item.to_dict() for item in raw_bundle.bundle.matrix_observations}
    projected_matrix_map = {item.name: item.to_dict() for item in projected_bundle.bundle.matrix_observations}
    for context, payload in counts.items():
        probabilities = _pair_probabilities(payload["counts"])
        probability_rows.append({"context": context, "p_pp": float(probabilities[0]), "p_pm": float(probabilities[1]), "p_mp": float(probabilities[2]), "p_mm": float(probabilities[3])})
        raw_covariances[context] = raw_matrix_map[f"cov_{context}"]
        projected_covariances[context] = projected_matrix_map[f"cov_{context}"]
    write_rows_csv(EXTRACTED / "pair_probabilities.csv", probability_rows)
    write_json(EXTRACTED / "raw_covariances.json", raw_covariances)
    write_json(EXTRACTED / "projected_covariances.json", projected_covariances)

    summary = {
        "raw_downstream_classification": raw_result.classification,
        "projected_downstream_classification": projected_result.classification,
        "raw_linear_residual": raw_result.residuals.get("linear_residual", 0.0),
        "projected_best_psd_margin": projected_result.residuals.get("best_psd_margin", 0.0),
    }
    audit = {
        "raw_ingestion": raw_bundle.audit,
        "projected_ingestion": projected_bundle.audit,
        "raw_assembly": assembly_result_to_dict(raw_assembly),
        "projected_assembly": assembly_result_to_dict(projected_assembly),
        "raw_downstream": descent_result_to_dict(raw_result),
        "projected_downstream": descent_result_to_dict(projected_result),
    }
    write_json(ROOT / "summary.json", summary)
    write_json(ROOT / "audit.json", audit)


if __name__ == "__main__":
    main()

