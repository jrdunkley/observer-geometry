from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from nomodescent import classify_relation

from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    TableObservation,
    assemble_problem_spec,
    encode_matrix_observation,
    infer_observer_candidates,
)
from worked_examples.common.simple_svg import bar_chart_svg

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def build_bundle() -> EvidenceBundle:
    task_table = TableObservation(
        name="benchmark_scores",
        columns=("arith", "logic", "memory"),
        rows=("benchmark_A", "benchmark_B"),
        values=((0.91, 0.52, 0.48), (0.54, 0.88, 0.49)),
        source_ref=SourceRef(source="benchmark_curated_table", location="Table 1", quote="Benchmark task means."),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        units="normalized score",
    )
    protocol_observations = (
        ProtocolObservation(
            name="benchmark_A",
            facts=("arith", "logic"),
            candidate_family=("A_balanced", "A_math_only"),
            source_ref=SourceRef(source="benchmark_curated_table", location="task inclusion"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
        ProtocolObservation(
            name="benchmark_B",
            facts=("logic", "memory"),
            candidate_family=("B_balanced", "B_reasoning_only"),
            source_ref=SourceRef(source="benchmark_curated_table", location="task inclusion"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
    )
    hypotheses = (
        ObserverHypothesis(
            name="A_balanced",
            matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            protocol_name="benchmark_A",
            features=("arith", "logic"),
            source_ref=SourceRef(source="benchmark_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("benchmark A is modelled as a two-direction math-plus-logic observer",),
        ),
        ObserverHypothesis(
            name="A_math_only",
            matrix=[[1.0, 0.0, 0.0]],
            protocol_name="benchmark_A",
            features=("arith",),
            source_ref=SourceRef(source="benchmark_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("benchmark A is modelled as math-only despite mixed task inclusion",),
        ),
        ObserverHypothesis(
            name="B_balanced",
            matrix=[[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            protocol_name="benchmark_B",
            features=("logic", "memory"),
            source_ref=SourceRef(source="benchmark_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("benchmark B is modelled as logic-plus-memory observer",),
        ),
        ObserverHypothesis(
            name="B_reasoning_only",
            matrix=[[0.0, 1.0, 0.0]],
            protocol_name="benchmark_B",
            features=("logic",),
            source_ref=SourceRef(source="benchmark_model_choice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("benchmark B is modelled as reasoning-only despite memory tasks being present",),
        ),
    )
    matrices = (
        encode_matrix_observation(
            name="cov_benchmark_A",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=((1.0, 0.21), (0.21, 1.0)),
            observer_name="benchmark_A",
            source="benchmark_curated_table",
            location="covariance estimate A",
            quote="Synthetic benchmark covariance A.",
        ),
        encode_matrix_observation(
            name="cov_benchmark_B",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=((1.0, 0.17), (0.17, 1.0)),
            observer_name="benchmark_B",
            source="benchmark_curated_table",
            location="covariance estimate B",
            quote="Synthetic benchmark covariance B.",
        ),
    )
    notes = (
        ExtractionNote(
            name="benchmark_scope",
            statement="This is a synthetic benchmark-encoding example, not a claim about real LLM capability structure.",
            load_bearing=True,
            source_ref=SourceRef(source="evidence_design"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=False),
        ),
    )
    return EvidenceBundle(
        name="benchmark_blindness_encoding",
        latent_dim=3,
        table_observations=(task_table,),
        matrix_observations=matrices,
        protocol_observations=protocol_observations,
        observer_hypotheses=hypotheses,
        notes=notes,
        description="Synthetic benchmark evidence encoding with explicit observer-family ambiguity",
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bundle = build_bundle()
    unresolved = assemble_problem_spec(bundle)
    suggestion_a = infer_observer_candidates(bundle, "benchmark_A")
    suggestion_b = infer_observer_candidates(bundle, "benchmark_B")
    selection = {"benchmark_A": suggestion_a.ranked_candidates[0], "benchmark_B": suggestion_b.ranked_candidates[0]}
    assembled = assemble_problem_spec(bundle, observer_selection=selection)
    if assembled.problem_spec is None:
        raise RuntimeError("benchmark bundle should assemble after explicit selection")

    observer_map = assembled.problem_spec.observer_map()
    relation = classify_relation(observer_map["A_balanced"], observer_map["B_balanced"])
    rows = [
        {"stage": "unresolved", "classification": unresolved.classification, "required_decisions": "; ".join(unresolved.required_human_decisions)},
        {"stage": "resolved", "classification": relation.classification, "required_decisions": ""},
    ]
    _write_csv(OUT / "assembly_path.csv", rows)
    bar_chart_svg(
        OUT / "suggestion_gap.svg",
        ["benchmark_A", "benchmark_B"],
        np.array(
            [
                suggestion_a.audit.residuals["score_gap_top2"],
                suggestion_b.audit.residuals["score_gap_top2"],
            ],
            dtype=float,
        ),
        title="Observer-Hypothesis Separation",
        y_label="top-two score gap",
        color="#355070",
    )
    summary = {
        "initial_classification": unresolved.classification,
        "selected_observers": selection,
        "resolved_relation": relation.classification,
    }
    audit = {
        "unresolved_audit": unresolved.audit.__dict__,
        "suggestion_audit": {"benchmark_A": suggestion_a.audit.__dict__, "benchmark_B": suggestion_b.audit.__dict__},
        "resolved_audit": assembled.audit.__dict__,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

