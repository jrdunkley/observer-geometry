from __future__ import annotations

import json
from pathlib import Path

import pytest

from evidence import EvidenceBundle, ExtractionRecord, SourceRef, TableObservation
from evidence.exceptions import EvidenceInputError


def test_exact_extraction_round_trip_preserves_values_and_provenance() -> None:
    table = TableObservation(
        name="scores",
        columns=("x", "y"),
        rows=("r1", "r2"),
        values=((1.0, 2.0), (3.0, 4.0)),
        source_ref=SourceRef(source="paper", location="Table 1", quote="numbers"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
    )
    bundle = EvidenceBundle(name="bundle", latent_dim=2, table_observations=(table,))
    path = Path(__file__).resolve().parent / "_bundle_roundtrip.json"
    try:
        bundle.to_json(path)
        restored = EvidenceBundle.from_json(path)
        assert restored.table_observations[0].values == ((1.0, 2.0), (3.0, 4.0))
        assert restored.table_observations[0].source_ref.location == "Table 1"
        assert json.loads(path.read_text(encoding="utf-8"))["table_observations"][0]["source_ref"]["quote"] == "numbers"
    finally:
        if path.exists():
            path.unlink()


def test_mislabelled_exact_extraction_is_rejected() -> None:
    with pytest.raises(EvidenceInputError, match="exact_extraction"):
        ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="inferred", authoritative=True)


def test_table_missing_required_indices_is_rejected() -> None:
    with pytest.raises(EvidenceInputError, match="column count"):
        TableObservation(
            name="bad",
            columns=("x", "y"),
            rows=("r1",),
            values=((1.0,),),
            source_ref=SourceRef(source="paper"),
            extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact"),
        )

