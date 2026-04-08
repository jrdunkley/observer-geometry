from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from nomodescent import ObserverSpec, classify_relation, common_descent_test, minimal_refinement_search

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
)

from micro_real_bundles.common import (
    assembly_result_to_dict,
    descent_result_to_dict,
    ensure_dir,
    refinement_result_to_dict,
    write_json,
    write_matrix_csv,
)

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
EXTRACTED = ROOT / "extracted"


def _load_rows() -> tuple[list[str], np.ndarray]:
    species: list[str] = []
    rows: list[list[float]] = []
    with (RAW / "iris_slice.csv").open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            species.append(row["species"])
            rows.append([float(row[key]) for key in ("sepal_length", "sepal_width", "petal_length", "petal_width")])
    return species, np.asarray(rows, dtype=float)


def _build_bundle() -> tuple[object, np.ndarray]:
    species, rows = _load_rows()
    full_covariance = np.cov(rows, rowvar=False, bias=False)
    excerpts = (
        ManualTableExcerpt(
            name="iris_slice",
            columns=("sepal_length", "sepal_width", "petal_length", "petal_width"),
            rows=tuple(f"{species_name}_{index}" for index, species_name in enumerate(species, start=1)),
            values=tuple(tuple(float(entry) for entry in row) for row in rows),
            source_ref=SourceRef(source="uci_iris", location="12-row curated slice"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual 12-row slice chosen for tiny redistributable protocol mismatch bundle",
            notes="raw flower measurements",
        ),
        ProtocolExcerpt(
            name="sepal_panel",
            facts=("sepal_length", "sepal_width"),
            candidate_family=("sepal_panel",),
            source_ref=SourceRef(source="uci_iris", location="measurement columns"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual panel encoding from explicit column names",
        ),
        ProtocolExcerpt(
            name="petal_panel",
            facts=("petal_length", "petal_width"),
            candidate_family=("petal_panel",),
            source_ref=SourceRef(source="uci_iris", location="measurement columns"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
            selection_provenance="manual panel encoding from explicit column names",
        ),
        MatrixExcerpt(
            name="cov_sepal",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in full_covariance[:2, :2]),
            observer_name="sepal_panel",
            source_ref=SourceRef(source="iris_protocol_mismatch", location="derived sepal covariance"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            selection_provenance="sample covariance computed exactly from extracted rows",
        ),
        MatrixExcerpt(
            name="cov_petal",
            matrix_role="visible_object",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in full_covariance[2:, 2:]),
            observer_name="petal_panel",
            source_ref=SourceRef(source="iris_protocol_mismatch", location="derived petal covariance"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            selection_provenance="sample covariance computed exactly from extracted rows",
        ),
        MatrixExcerpt(
            name="cov_full",
            matrix_role="latent_matrix",
            matrix_kind="covariance",
            matrix=tuple(tuple(float(entry) for entry in row) for row in full_covariance),
            observer_name=None,
            source_ref=SourceRef(source="iris_protocol_mismatch", location="derived full covariance"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=False),
            selection_provenance="full covariance retained only as an independent validation object",
        ),
    )
    observers = (
        ObserverHypothesis(
            name="sepal_panel",
            matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
            protocol_name="sepal_panel",
            features=("sepal_length", "sepal_width"),
            source_ref=SourceRef(source="uci_iris", location="measurement columns"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
        ObserverHypothesis(
            name="petal_panel",
            matrix=[[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            protocol_name="petal_panel",
            features=("petal_length", "petal_width"),
            source_ref=SourceRef(source="uci_iris", location="measurement columns"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
    )
    notes = (
        ExtractionNote(
            name="sample_covariance_note",
            statement="Visible covariances are exact sample covariances of the extracted rows, not matrices written explicitly in the source table.",
            load_bearing=True,
            source_ref=SourceRef(source="iris_protocol_mismatch", location="bundle design"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred"),
        ),
    )
    bundle = build_evidence_bundle_from_excerpts(
        name="iris_protocol_mismatch",
        latent_dim=4,
        excerpts=excerpts,
        observer_hypotheses=observers,
        notes=notes,
        description="Real measurement-panel mismatch bundle from a tiny Iris slice",
        provenance="UCI Iris 12-row curated slice",
        tags=("micro-real", "protocol-mismatch"),
    )
    return bundle, full_covariance


def main() -> None:
    ensure_dir(EXTRACTED)
    bundle_result, full_covariance = _build_bundle()
    bundle_result.bundle.to_json(ROOT / "bundle.json")
    assembly = assemble_problem_spec(bundle_result.bundle)
    if assembly.problem_spec is None:
        raise RuntimeError("Iris bundle should assemble directly")
    assembly.problem_spec.to_json(ROOT / "problem_spec.json")

    observer_map = assembly.problem_spec.observer_map()
    relation = classify_relation(observer_map["sepal_panel"], observer_map["petal_panel"])
    completion = common_descent_test(assembly.problem_spec)
    refinement = minimal_refinement_search(
        np.linalg.inv(full_covariance),
        [observer_map["sepal_panel"], observer_map["petal_panel"]],
        [
            ObserverSpec(name="sepal_panel", matrix=observer_map["sepal_panel"].matrix),
            ObserverSpec(name="petal_panel", matrix=observer_map["petal_panel"].matrix),
            ObserverSpec(name="full_iris", matrix=np.eye(4)),
        ],
    )

    write_matrix_csv(EXTRACTED / "cov_sepal.csv", full_covariance[:2, :2], ["sepal_length", "sepal_width"], ["sepal_length", "sepal_width"])
    write_matrix_csv(EXTRACTED / "cov_petal.csv", full_covariance[2:, 2:], ["petal_length", "petal_width"], ["petal_length", "petal_width"])
    write_matrix_csv(EXTRACTED / "cov_full.csv", full_covariance, ["sepal_length", "sepal_width", "petal_length", "petal_width"], ["sepal_length", "sepal_width", "petal_length", "petal_width"])

    summary = {
        "assembly_classification": assembly.classification,
        "relation_classification": relation.classification,
        "completion_classification": completion.classification,
        "refinement_winner": refinement.winner,
    }
    audit = {
        "ingestion": bundle_result.audit,
        "assembly": assembly_result_to_dict(assembly),
        "relation": descent_result_to_dict(relation),
        "completion": descent_result_to_dict(completion),
        "refinement": refinement_result_to_dict(refinement),
    }
    write_json(ROOT / "summary.json", summary)
    write_json(ROOT / "audit.json", audit)


if __name__ == "__main__":
    main()

