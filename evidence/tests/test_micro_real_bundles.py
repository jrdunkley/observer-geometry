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


def test_bell_micro_real_bundle_runs_cleanly() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "micro_real_bundles" / "bell_counts_bundle" / "run_main.py",
        root / "micro_real_bundles" / "bell_counts_bundle" / "summary.json",
    )
    assert summary["raw_downstream_classification"] == "incompatible_by_linear_inconsistency"
    assert summary["projected_downstream_classification"] == "approximate_common_descent"


def test_iris_micro_real_bundle_runs_cleanly() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "micro_real_bundles" / "iris_protocol_mismatch" / "run_main.py",
        root / "micro_real_bundles" / "iris_protocol_mismatch" / "summary.json",
    )
    assert summary["relation_classification"] == "non_nested_observers"
    assert summary["completion_classification"] == "underdetermined_affine_family"
    assert summary["refinement_winner"] == "full_iris"


def test_benchmark_micro_real_bundle_runs_cleanly() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "micro_real_bundles" / "leaderboard_benchmark_slice" / "run_main.py",
        root / "micro_real_bundles" / "leaderboard_benchmark_slice" / "summary.json",
    )
    assert summary["initial_assembly_classification"] == "underdetermined_evidence"
    assert summary["relation_classification"] == "non_nested_observers"
    assert summary["completion_classification"] == "underdetermined_affine_family"
