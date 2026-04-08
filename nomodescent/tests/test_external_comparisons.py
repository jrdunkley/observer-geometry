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


def test_bell_external_comparison_detects_false_collapse() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "external_comparisons" / "bell_chsh_comparison" / "run_main.py",
        root / "external_comparisons" / "bell_chsh_comparison" / "outputs" / "summary.json",
    )
    assert summary["external_method_result"] == "non_obstructed"
    assert summary["raw_qd_classification"] == "incompatible_by_linear_inconsistency"
    assert summary["false_collapse_classification"] == "false_collapse_detected"


def test_protocol_cca_external_comparison_distinguishes_relation_and_completion() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "external_comparisons" / "protocol_cca_comparison" / "run_main.py",
        root / "external_comparisons" / "protocol_cca_comparison" / "outputs" / "summary.json",
    )
    assert summary["qd_relation_classification"] == "non_nested_observers"
    assert summary["qd_completion_classification"] == "underdetermined_affine_family"


def test_free_gaussian_rg_external_comparison_reproduces_schur_elimination() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = _run_and_load(
        root / "external_comparisons" / "free_gaussian_rg_comparison" / "run_main.py",
        root / "external_comparisons" / "free_gaussian_rg_comparison" / "outputs" / "summary.json",
    )
    assert summary["schur_residual"] < 1e-12
    assert summary["tower_classification"] == "exact_tower_agreement"
