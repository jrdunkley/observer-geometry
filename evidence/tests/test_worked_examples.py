from __future__ import annotations

import json
from pathlib import Path
import runpy
import sys


def _run_and_load(script_path: Path, summary_path: Path) -> dict[str, object]:
    sys.path.insert(0, str(script_path.parents[2]))
    try:
        namespace = runpy.run_path(str(script_path))
        namespace["main"]()
    finally:
        sys.path.pop(0)
    return json.loads(summary_path.read_text(encoding="utf-8"))


def test_bell_encoding_reproduces_downstream_split() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "worked_examples" / "bell_evidence_encoding" / "run_main.py",
        root / "worked_examples" / "bell_evidence_encoding" / "outputs" / "summary.json",
    )
    assert summary["coarse_only_classification"]["compatible"] == "insufficient_evidence"
    assert summary["full_law_downstream_classification"]["variance_only"] == "incompatible_by_linear_inconsistency"


def test_replication_encoding_reproduces_non_nested_result() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "worked_examples" / "replication_protocol_encoding" / "run_main.py",
        root / "worked_examples" / "replication_protocol_encoding" / "outputs" / "summary.json",
    )
    assert summary["assembly_classification"] == "assembled_problem_spec"
    assert summary["downstream_relation"] == "non_nested_observers"


def test_benchmark_encoding_preserves_initial_underdetermination() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "worked_examples" / "benchmark_blindness_encoding" / "run_main.py",
        root / "worked_examples" / "benchmark_blindness_encoding" / "outputs" / "summary.json",
    )
    assert summary["initial_classification"] == "underdetermined_evidence"
    assert summary["resolved_relation"] == "non_nested_observers"
