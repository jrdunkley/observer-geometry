from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from nomodescent import classify_relation, common_descent_test

from evidence import (
    ExtractionNote,
    ExtractionRecord,
    ManualTableExcerpt,
    MatrixExcerpt,
    ObserverHypothesis,
    ProtocolExcerpt,
    SourceRef,
    assemble_problem_spec,
    build_evidence_bundle_from_excerpts,
    infer_observer_candidates,
)

from micro_real_bundles.common import (
    assembly_result_to_dict,
    descent_result_to_dict,
    ensure_dir,
    suggestion_result_to_dict,
    write_json,
    write_matrix_csv,
)

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
EXTRACTED = ROOT / "extracted"

TASKS = ("ifeval", "bbh", "math_lvl5", "mmlu_pro")


def _load_scores() -> tuple[list[str], np.ndarray]:
    models: list[str] = []
    rows: list[list[float]] = []
    with (RAW / "score_slice.csv").open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            models.append(row["model"])
            rows.append([float(row[key]) for key in TASKS])
    return models, np.asarray(rows, dtype=float)


def _build_bundle() -> tuple[object, np.ndarray]:
    models, scores = _load_scores()
    covariance = np.cov(scores, rowvar=False, bias=False)
    excerpts = (
        ManualTableExcerpt(
            name="leaderboard_score_slice",
            columns=TASKS,
            rows=tuple(models),
            values=tuple(tuple(float(entry) for entry in row) for row in scores),
            source_ref=SourceRef(source="open_llm_leaderboard_public_page", location="manual 3x4 score slice"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual transcription of three public leaderboard rows and four task columns",
            units="score",
        ),
        ProtocolExcerpt(
            name="chat_suite",
            facts=("ifeval", "bbh"),
            candidate_family=("chat_direct", "chat_blurred"),
            source_ref=SourceRef(source="open_llm_leaderboard_public_page", location="task inclusion"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual benchmark-suite encoding from task subset",
        ),
        ProtocolExcerpt(
            name="reasoning_suite",
            facts=("bbh", "mmlu_pro"),
            candidate_family=("reasoning_direct", "reasoning_blurred"),
            source_ref=SourceRef(source="open_llm_leaderboard_public_page", location="task inclusion"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual benchmark-suite encoding from task subset",
        ),
        MatrixExcerpt(
            name="cov_chat",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in covariance[np.ix_([0, 1], [0, 1])]),
            observer_name="chat_suite",
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="derived chat covariance"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            selection_provenance="sample covariance computed from the exact score slice over IFEval and BBH",
        ),
        MatrixExcerpt(
            name="cov_reasoning",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in covariance[np.ix_([1, 3], [1, 3])]),
            observer_name="reasoning_suite",
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="derived reasoning covariance"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            selection_provenance="sample covariance computed from the exact score slice over BBH and MMLU-PRO",
        ),
    )
    observers = (
        ObserverHypothesis(
            name="chat_direct",
            matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
            protocol_name="chat_suite",
            features=("ifeval", "bbh"),
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("chat suite is modelled as direct observation of IFEval and BBH directions",),
        ),
        ObserverHypothesis(
            name="chat_blurred",
            matrix=[[1.0, 0.0, 0.0, 0.0], [0.4, 0.6, 0.0, 0.0]],
            protocol_name="chat_suite",
            features=("ifeval", "bbh", "style"),
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("chat suite is modelled as an instruction axis plus a blurred chat-quality axis",),
        ),
        ObserverHypothesis(
            name="reasoning_direct",
            matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            protocol_name="reasoning_suite",
            features=("bbh", "mmlu_pro"),
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("reasoning suite is modelled as direct observation of BBH and MMLU-PRO directions",),
        ),
        ObserverHypothesis(
            name="reasoning_blurred",
            matrix=[[0.0, 0.5, 0.5, 0.0], [0.0, 0.0, 0.0, 1.0]],
            protocol_name="reasoning_suite",
            features=("bbh", "math_lvl5", "mmlu_pro"),
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("reasoning suite is modelled as a blurred reasoning-plus-math axis with an MMLU-PRO axis",),
        ),
    )
    notes = (
        ExtractionNote(
            name="latent_axis_choice",
            statement="The latent ambient space is identified with the four task-score directions of the reconstructed score slice.",
            load_bearing=True,
            source_ref=SourceRef(source="leaderboard_benchmark_slice", location="bundle design"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred"),
        ),
    )
    bundle = build_evidence_bundle_from_excerpts(
        name="leaderboard_benchmark_slice",
        latent_dim=4,
        excerpts=excerpts,
        observer_hypotheses=observers,
        notes=notes,
        description="Tiny public leaderboard score slice with explicit benchmark-observer ambiguity",
        provenance="Manual reconstruction from public leaderboard numeric facts",
        tags=("micro-real", "benchmark"),
    )
    return bundle, covariance


def main() -> None:
    ensure_dir(EXTRACTED)
    bundle_result, covariance = _build_bundle()
    bundle_result.bundle.to_json(ROOT / "bundle.json")

    unresolved = assemble_problem_spec(bundle_result.bundle)
    suggestion_chat = infer_observer_candidates(bundle_result.bundle, "chat_suite")
    suggestion_reasoning = infer_observer_candidates(bundle_result.bundle, "reasoning_suite")
    selection = {"chat_suite": suggestion_chat.ranked_candidates[0], "reasoning_suite": suggestion_reasoning.ranked_candidates[0]}
    assembled = assemble_problem_spec(bundle_result.bundle, observer_selection=selection)
    if assembled.problem_spec is None:
        raise RuntimeError("benchmark bundle should assemble after explicit observer selection")
    assembled.problem_spec.to_json(ROOT / "problem_spec.json")

    relation = classify_relation(
        assembled.problem_spec.observer_map()[selection["chat_suite"]],
        assembled.problem_spec.observer_map()[selection["reasoning_suite"]],
    )
    completion = common_descent_test(assembled.problem_spec)

    write_matrix_csv(EXTRACTED / "covariance_full.csv", covariance, ["ifeval", "bbh", "math_lvl5", "mmlu_pro"], ["ifeval", "bbh", "math_lvl5", "mmlu_pro"])
    write_matrix_csv(EXTRACTED / "covariance_chat.csv", covariance[np.ix_([0, 1], [0, 1])], ["ifeval", "bbh"], ["ifeval", "bbh"])
    write_matrix_csv(EXTRACTED / "covariance_reasoning.csv", covariance[np.ix_([1, 3], [1, 3])], ["bbh", "mmlu_pro"], ["bbh", "mmlu_pro"])

    summary = {
        "initial_assembly_classification": unresolved.classification,
        "selected_observers": selection,
        "relation_classification": relation.classification,
        "completion_classification": completion.classification,
    }
    audit = {
        "ingestion": bundle_result.audit,
        "unresolved_assembly": assembly_result_to_dict(unresolved),
        "suggestion_chat": suggestion_result_to_dict(suggestion_chat),
        "suggestion_reasoning": suggestion_result_to_dict(suggestion_reasoning),
        "assembled": assembly_result_to_dict(assembled),
        "relation": descent_result_to_dict(relation),
        "completion": descent_result_to_dict(completion),
    }
    write_json(ROOT / "summary.json", summary)
    write_json(ROOT / "audit.json", audit)


if __name__ == "__main__":
    main()

