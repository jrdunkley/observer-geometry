from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from nomodescent import common_descent_test

from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    TextClaim,
    assemble_problem_spec,
    encode_matrix_observation,
    encode_table_observation,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def observers() -> tuple[ObserverHypothesis, ...]:
    source = SourceRef(source="bell_curated_source", location="observer section")
    extraction = ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True)
    return (
        ObserverHypothesis(name="00", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]], protocol_name="00", features=("A0", "B0"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="01", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]], protocol_name="01", features=("A0", "B1"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="10", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]], protocol_name="10", features=("A1", "B0"), source_ref=source, extraction=extraction),
        ObserverHypothesis(name="11", matrix=[[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]], protocol_name="11", features=("A1", "B1"), source_ref=source, extraction=extraction),
    )


def family_covariances(delta: float, rho: float) -> dict[str, np.ndarray]:
    return {
        "00": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
        "01": np.array([[1.0 + delta, rho * np.sqrt(1.0 + delta)], [rho * np.sqrt(1.0 + delta), 1.0]], dtype=float),
        "10": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
        "11": np.array([[1.0, -rho], [-rho, 1.0]], dtype=float),
    }


def correlator_matrix(family: dict[str, np.ndarray]) -> np.ndarray:
    return np.array(
        [
            [family["00"][0, 1] / np.sqrt(family["00"][0, 0] * family["00"][1, 1]), family["01"][0, 1] / np.sqrt(family["01"][0, 0] * family["01"][1, 1])],
            [family["10"][0, 1] / np.sqrt(family["10"][0, 0] * family["10"][1, 1]), family["11"][0, 1] / np.sqrt(family["11"][0, 0] * family["11"][1, 1])],
        ],
        dtype=float,
    )


def build_bundle(name: str, delta: float, rho: float, *, include_full_law: bool) -> EvidenceBundle:
    family = family_covariances(delta, rho)
    tables = (
        encode_table_observation(
            name=f"{name}_correlators",
            columns=("B0", "B1"),
            rows=("A0", "A1"),
            values=tuple(tuple(float(entry) for entry in row) for row in correlator_matrix(family)),
            source="bell_curated_source",
            location="correlator summary table",
            quote="Normalized correlator matrix.",
            authoritative=False,
        ),
    )
    matrices = ()
    if include_full_law:
        matrices = tuple(
            encode_matrix_observation(
                name=f"sigma_{name}_{key}",
                matrix_role="visible_object",
                matrix_kind="covariance",
                matrix=tuple(tuple(float(entry) for entry in row) for row in matrix),
                observer_name=key,
                source="bell_curated_source",
                location=f"pair law {key}",
                quote="Pairwise Gaussian covariance matrix.",
                authoritative=True,
            )
            for key, matrix in family.items()
        )
    protocols = tuple(
        ProtocolObservation(
            name=key,
            facts=(f"A{key[0]}", f"B{key[1]}"),
            candidate_family=(key,),
            source_ref=SourceRef(source="bell_curated_source", location=f"context {key}"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        )
        for key in ("00", "01", "10", "11")
    )
    notes = (
        ExtractionNote(
            name="coarse_not_full_law",
            statement="Normalized correlator summaries are preserved separately from full Gaussian law evidence.",
            load_bearing=True,
            source_ref=SourceRef(source="evidence_design"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=False),
        ),
    )
    text_claims = (
        TextClaim(
            name="latent_dim",
            statement="Common latent Gaussian object is encoded in four Bell variables.",
            source_ref=SourceRef(source="bell_curated_source", location="setup note"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
        ),
    )
    return EvidenceBundle(
        name=name,
        latent_dim=4,
        text_claims=text_claims,
        table_observations=tables,
        matrix_observations=matrices,
        protocol_observations=protocols,
        observer_hypotheses=observers(),
        notes=notes,
        description="Bell square evidence bundle with explicit separation between correlator summaries and full pair laws",
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    audits: dict[str, object] = {}
    sample_specs = {
        "compatible": (0.0, 0.6, True),
        "variance_only": (0.22, 0.6, False),
        "correlator_only": (0.0, 0.82, True),
    }
    for sample, (delta, rho, allow_approx) in sample_specs.items():
        full_bundle = build_bundle(sample, delta, rho, include_full_law=True)
        coarse_bundle = build_bundle(f"{sample}_coarse_only", delta, rho, include_full_law=False)

        full_assembly = assemble_problem_spec(full_bundle)
        coarse_assembly = assemble_problem_spec(coarse_bundle)
        if full_assembly.problem_spec is None:
            raise RuntimeError("full Bell bundle should assemble")
        downstream = common_descent_test(full_assembly.problem_spec, allow_approximate_psd_search=allow_approx, psd_search_grid_size=31)

        rows.append(
            {
                "sample": sample,
                "full_law_assembly": full_assembly.classification,
                "correlator_only_assembly": coarse_assembly.classification,
                "downstream_classification": downstream.classification,
                "downstream_exact": downstream.exact,
                "linear_residual": downstream.residuals.get("linear_residual", 0.0),
            }
        )
        audits[sample] = {
            "full_assembly_audit": full_assembly.audit.__dict__,
            "coarse_assembly_audit": coarse_assembly.audit.__dict__,
            "downstream_theorem_layer": downstream.audit.theorem_layer,
            "downstream_certificates": [cert.kind for cert in downstream.certificates],
        }

    _write_csv(OUT / "encoding_summary.csv", rows)
    summary = {
        "coarse_only_classification": {row["sample"]: row["correlator_only_assembly"] for row in rows},
        "full_law_downstream_classification": {row["sample"]: row["downstream_classification"] for row in rows},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "audit.json").write_text(json.dumps(audits, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

