from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from nomogeo import BatchTaskError, batch_map

import nomogeo.batch as batch_module

from tools.stack_soak_support import execute_soak_task, fail_soak_task
from tools.validation_sweep import run_validation_sweep

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tools" / "outputs"


@dataclass(frozen=True)
class StackSoakSummary:
    seed: int
    worker_modes_exercised: tuple[str, ...]
    task_counts: dict[str, int]
    backend_runtime_seconds: dict[str, float]
    backend_runtime_by_suite_seconds: dict[str, dict[str, float]]
    backend_agreement: dict[str, dict[str, float | int | bool]]
    failure_propagation: dict[str, dict[str, object]]
    forced_process_fallback_matches_serial: bool
    long_run_kernel_dashboard: dict[str, float]
    long_run_descent_dashboard: dict[str, float]
    long_run_evidence_dashboard: dict[str, float]
    long_run_example_dashboard: dict[str, float]
    validation_sweep: dict[str, object]


def run_stack_soak(
    *,
    seed: int = 20260408,
    kernel_count: int = 96,
    descent_count: int = 48,
    evidence_count: int = 24,
    example_count: int = 24,
    max_workers: int = 2,
) -> StackSoakSummary:
    suites = _build_task_suites(seed, kernel_count, descent_count, evidence_count, example_count)
    backends = ("serial", "thread", "process")
    suite_results: dict[str, dict[str, list[dict[str, object]]]] = {}
    backend_runtime_seconds: dict[str, float] = {}
    backend_runtime_by_suite_seconds: dict[str, dict[str, float]] = {}

    for backend in backends:
        backend_total_start = time.perf_counter()
        suite_results[backend] = {}
        backend_runtime_by_suite_seconds[backend] = {}
        for suite_name, tasks in suites.items():
            start = time.perf_counter()
            suite_results[backend][suite_name] = batch_map(
                execute_soak_task,
                tasks,
                backend=backend,
                max_workers=max_workers,
            )
            backend_runtime_by_suite_seconds[backend][suite_name] = float(time.perf_counter() - start)
        backend_runtime_seconds[backend] = float(time.perf_counter() - backend_total_start)

    backend_agreement = {
        "serial_vs_thread": _compare_backend_results(suite_results["serial"], suite_results["thread"]),
        "serial_vs_process": _compare_backend_results(suite_results["serial"], suite_results["process"]),
    }
    failure_propagation = {backend: _failure_probe(backend, max_workers=max_workers) for backend in backends}
    forced_process_fallback_matches_serial = _forced_process_fallback_matches(suites["kernel"][: min(8, len(suites["kernel"]))], max_workers=max_workers)

    kernel_dashboard = _suite_dashboard(suite_results["serial"]["kernel"])
    descent_dashboard = _suite_dashboard(suite_results["serial"]["descent"])
    evidence_dashboard = _suite_dashboard(suite_results["serial"]["evidence"])
    example_dashboard = _suite_dashboard(suite_results["serial"]["example"])
    validation_sweep = asdict(run_validation_sweep(seed=seed, dims=list(range(2, 9)), draws_per_dim=6))

    summary = StackSoakSummary(
        seed=seed,
        worker_modes_exercised=backends,
        task_counts={suite: len(tasks) for suite, tasks in suites.items()},
        backend_runtime_seconds=backend_runtime_seconds,
        backend_runtime_by_suite_seconds=backend_runtime_by_suite_seconds,
        backend_agreement=backend_agreement,
        failure_propagation=failure_propagation,
        forced_process_fallback_matches_serial=forced_process_fallback_matches_serial,
        long_run_kernel_dashboard=kernel_dashboard,
        long_run_descent_dashboard=descent_dashboard,
        long_run_evidence_dashboard=evidence_dashboard,
        long_run_example_dashboard=example_dashboard,
        validation_sweep=validation_sweep,
    )
    return summary


def write_stack_soak_outputs(summary: StackSoakSummary) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "stack_soak_summary.json").write_text(json.dumps(asdict(summary), indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    summary = run_stack_soak()
    write_stack_soak_outputs(summary)
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


def _build_task_suites(seed: int, kernel_count: int, descent_count: int, evidence_count: int, example_count: int) -> dict[str, list[dict[str, object]]]:
    return {
        "kernel": [
            {"task_kind": "kernel", "seed": seed + idx, "variant": "general" if idx % 4 else "boundary"}
            for idx in range(kernel_count)
        ],
        "descent": [
            {
                "task_kind": "descent",
                "seed": seed + 10_000 + idx,
                "variant": ("factorisation", "tower", "completion_exact", "completion_incompatible", "refinement")[idx % 5],
            }
            for idx in range(descent_count)
        ],
        "evidence": [
            {
                "task_kind": "evidence",
                "seed": seed + 20_000 + idx,
                "variant": ("underdetermined", "assembled", "suggestion")[idx % 3],
            }
            for idx in range(evidence_count)
        ],
        "example": [
            {
                "task_kind": "example",
                "seed": seed + 30_000 + idx,
                "variant": ("entanglement", "bell", "arrow")[idx % 3],
            }
            for idx in range(example_count)
        ],
    }


def _compare_backend_results(
    reference: dict[str, list[dict[str, object]]],
    candidate: dict[str, list[dict[str, object]]],
) -> dict[str, float | int | bool]:
    max_numeric_diff = 0.0
    classification_mismatch_count = 0
    ordering_mismatch = False
    for suite_name, ref_rows in reference.items():
        cand_rows = candidate[suite_name]
        if len(ref_rows) != len(cand_rows):
            ordering_mismatch = True
            classification_mismatch_count += abs(len(ref_rows) - len(cand_rows))
            continue
        for ref_row, cand_row in zip(ref_rows, cand_rows, strict=True):
            if ref_row.get("task_kind") != cand_row.get("task_kind") or ref_row.get("variant") != cand_row.get("variant") or ref_row.get("seed") != cand_row.get("seed"):
                ordering_mismatch = True
            keys = set(ref_row) | set(cand_row)
            for key in keys:
                ref_value = ref_row.get(key)
                cand_value = cand_row.get(key)
                if isinstance(ref_value, bool) or isinstance(cand_value, bool):
                    if ref_value != cand_value:
                        classification_mismatch_count += 1
                    continue
                if isinstance(ref_value, (int, float)) and isinstance(cand_value, (int, float)):
                    max_numeric_diff = max(max_numeric_diff, abs(float(ref_value) - float(cand_value)))
                    continue
                if ref_value != cand_value:
                    classification_mismatch_count += 1
    return {
        "ordering_match": not ordering_mismatch,
        "classification_mismatch_count": classification_mismatch_count,
        "max_numeric_difference": max_numeric_diff,
        "all_match": (not ordering_mismatch) and classification_mismatch_count == 0,
    }


def _failure_probe(backend: str, *, max_workers: int) -> dict[str, object]:
    tasks = [
        {"label": "ok_0", "fail": False},
        {"label": "bad_1", "fail": True},
        {"label": "ok_2", "fail": False},
    ]
    try:
        batch_map(fail_soak_task, tasks, backend=backend, max_workers=max_workers)
    except BatchTaskError as exc:
        return {
            "index": exc.index,
            "task": exc.task,
            "message": str(exc),
        }
    raise RuntimeError(f"failure probe unexpectedly passed under backend '{backend}'")


def _forced_process_fallback_matches(tasks: list[dict[str, object]], *, max_workers: int) -> bool:
    serial = batch_map(execute_soak_task, tasks, backend="serial", max_workers=max_workers)
    original = batch_module.ProcessPoolExecutor

    class BlockedProcessPool:
        def __init__(self, *args: object, **kwargs: object) -> None:
            raise PermissionError("forced process blockage")

    batch_module.ProcessPoolExecutor = BlockedProcessPool
    try:
        fallback = batch_map(execute_soak_task, tasks, backend="process", max_workers=max_workers)
    finally:
        batch_module.ProcessPoolExecutor = original
    return serial == fallback


def _suite_dashboard(rows: list[dict[str, object]]) -> dict[str, float]:
    numeric_by_key: dict[str, list[float]] = {}
    for row in rows:
        for key, value in row.items():
            if key in {"seed"}:
                continue
            if isinstance(value, bool):
                numeric_by_key.setdefault(key, []).append(1.0 if value else 0.0)
            elif isinstance(value, (int, float)):
                numeric_by_key.setdefault(key, []).append(float(value))
    dashboard: dict[str, float] = {}
    for key, values in numeric_by_key.items():
        array = np.asarray(values, dtype=float)
        dashboard[f"max_{key}"] = float(np.max(array))
        dashboard[f"median_{key}"] = float(np.median(array))
    return dashboard


if __name__ == "__main__":
    main()
