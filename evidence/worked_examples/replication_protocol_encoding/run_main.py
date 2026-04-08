from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from nomogeo import gaussian_bhattacharyya_distance, observed_covariance
from nomodescent import ObserverSpec, classify_relation, minimal_refinement_search

from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    assemble_problem_spec,
    encode_matrix_observation,
    infer_observer_candidates,
)
from worked_examples.common.simple_svg import bar_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def latent_precision() -> np.ndarray:
    return np.array(
        [
            [3.3, 0.9, 0.6, 0.0],
            [0.9, 2.9, 0.4, 0.1],
            [0.6, 0.4, 2.2, 0.0],
            [0.0, 0.1, 0.0, 1.7],
        ],
        dtype=float,
    )


def build_bundle() -> EvidenceBundle:
    H = latent_precision()
    protocol_a_matrix = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]], dtype=float)
    protocol_b_matrix = np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0]], dtype=float)

    protocol_observations = (
        ProtocolObservation(
            name="protocol_a",
            facts=("measure_x1", "measure_x2"),
            candidate_family=("protocol_a_direct", "protocol_a_sum13"),
            source_ref=SourceRef(source="replication_curated_source", location="protocol A"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
        ProtocolObservation(
            name="protocol_b",
            facts=("measure_x2", "aggregate_x1_x3"),
            candidate_family=("protocol_b_sum13", "protocol_b_direct"),
            source_ref=SourceRef(source="replication_curated_source", location="protocol B"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
    )
    observer_hypotheses = (
        ObserverHypothesis(
            name="protocol_a_direct",
            matrix=protocol_a_matrix,
            protocol_name="protocol_a",
            features=("measure_x1", "measure_x2"),
            source_ref=SourceRef(source="replication_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol A is encoded as direct observation of x1 and x2",),
        ),
        ObserverHypothesis(
            name="protocol_a_sum13",
            matrix=protocol_b_matrix,
            protocol_name="protocol_a",
            features=("measure_x2", "aggregate_x1_x3"),
            source_ref=SourceRef(source="replication_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol A is encoded as an aggregate x1+x3 panel, not direct x1",),
        ),
        ObserverHypothesis(
            name="protocol_b_sum13",
            matrix=protocol_b_matrix,
            protocol_name="protocol_b",
            features=("measure_x2", "aggregate_x1_x3"),
            source_ref=SourceRef(source="replication_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol B is encoded as an aggregate x1+x3 panel with x2",),
        ),
        ObserverHypothesis(
            name="protocol_b_direct",
            matrix=protocol_a_matrix,
            protocol_name="protocol_b",
            features=("measure_x1", "measure_x2"),
            source_ref=SourceRef(source="replication_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol B is encoded as direct x1 and x2 only",),
        ),
    )
    matrices = (
        encode_matrix_observation(
            name="cov_protocol_a",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in observed_covariance(H, protocol_a_matrix)),
            observer_name="protocol_a",
            source="replication_curated_source",
            location="protocol A result table",
            quote="Reported covariance under protocol A.",
        ),
        encode_matrix_observation(
            name="cov_protocol_b",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in observed_covariance(H, protocol_b_matrix)),
            observer_name="protocol_b",
            source="replication_curated_source",
            location="protocol B result table",
            quote="Reported covariance under protocol B.",
        ),
    )
    notes = (
        ExtractionNote(
            name="observer_family_choice",
            statement="Protocol descriptions do not directly state observer matrices; they are encoded from a finite candidate family.",
            load_bearing=True,
            source_ref=SourceRef(source="evidence_design"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=False),
        ),
    )
    return EvidenceBundle(
        name="replication_protocol_encoding",
        latent_dim=4,
        matrix_observations=matrices,
        protocol_observations=protocol_observations,
        observer_hypotheses=observer_hypotheses,
        notes=notes,
        description="Protocol mismatch as finite observer-family encoding",
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    H = latent_precision()
    bundle = build_bundle()
    suggestion_a = infer_observer_candidates(bundle, "protocol_a")
    suggestion_b = infer_observer_candidates(bundle, "protocol_b")
    selection = {"protocol_a": suggestion_a.ranked_candidates[0], "protocol_b": suggestion_b.ranked_candidates[0]}
    assembly = assemble_problem_spec(bundle, observer_selection=selection)
    if assembly.problem_spec is None:
        raise RuntimeError("replication bundle should assemble after explicit observer selection")

    observer_map = assembly.problem_spec.observer_map()
    relation = classify_relation(observer_map["protocol_a_direct"], observer_map["protocol_b_sum13"])
    candidates = [
        observer_map["protocol_a_direct"],
        observer_map["protocol_b_sum13"],
    ]
    refinements = [
        ObserverSpec(name="panel_3", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="full_4", matrix=np.eye(4)),
    ]
    search = minimal_refinement_search(H, candidates, refinements)

    cov_a = assembly.problem_spec.evidence_map()["cov_protocol_a"].matrix
    cov_b = assembly.problem_spec.evidence_map()["cov_protocol_b"].matrix
    mismatch = gaussian_bhattacharyya_distance(cov_a, cov_b)

    rows = [
        {"target": "protocol_a", "winner": suggestion_a.ranked_candidates[0], "score": suggestion_a.scores[suggestion_a.ranked_candidates[0]]},
        {"target": "protocol_b", "winner": suggestion_b.ranked_candidates[0], "score": suggestion_b.scores[suggestion_b.ranked_candidates[0]]},
    ]
    _write_csv(OUT / "suggestion_summary.csv", rows)
    _write_csv(
        OUT / "candidate_scores.csv",
        [{"candidate": name, "score": score} for name, score in search.scores.items()],
    )
    bar_chart_svg(
        OUT / "suggestion_scores.svg",
        ["protocol_a", "protocol_b"],
        np.array([rows[0]["score"], rows[1]["score"]], dtype=float),
        title="Top Observer-Hypothesis Scores",
        y_label="finite-family score",
        color="#2a9d8f",
    )
    summary = {
        "assembly_classification": assembly.classification,
        "selected_observers": selection,
        "downstream_relation": relation.classification,
        "mismatch_score": float(mismatch),
        "refinement_winner": search.winner,
    }
    audit = {
        "assembly_audit": assembly.audit.__dict__,
        "suggestion_audit": {"protocol_a": suggestion_a.audit.__dict__, "protocol_b": suggestion_b.audit.__dict__},
        "refinement_falsification_route": list(search.audit.falsification_route),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

