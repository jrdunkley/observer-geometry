from __future__ import annotations

import csv
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

import numpy as np


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(to_jsonable(payload), indent=2, sort_keys=True), encoding="utf-8")


def write_rows_csv(path: Path, rows: list[dict[str, object]]) -> None:
    ensure_dir(path.parent)
    if not rows:
        raise ValueError("rows must be non-empty")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_matrix_csv(path: Path, matrix: np.ndarray, row_labels: list[str], column_labels: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row", *column_labels])
        array = np.asarray(matrix, dtype=float)
        for label, row in zip(row_labels, array, strict=True):
            writer.writerow([label, *[float(value) for value in row]])


def descent_result_to_dict(result: Any) -> dict[str, object]:
    return {
        "classification": result.classification,
        "exact": result.exact,
        "residuals": to_jsonable(result.residuals),
        "certificates": [to_jsonable(certificate) for certificate in result.certificates],
        "audit": to_jsonable(result.audit),
        "details": to_jsonable(result.details),
    }


def assembly_result_to_dict(result: Any) -> dict[str, object]:
    return {
        "classification": result.classification,
        "exact": result.exact,
        "selected_observers": list(result.selected_observers),
        "unresolved_observer_groups": list(result.unresolved_observer_groups),
        "required_human_decisions": list(result.required_human_decisions),
        "audit": to_jsonable(result.audit),
        "details": to_jsonable(result.details),
    }


def suggestion_result_to_dict(result: Any) -> dict[str, object]:
    return {
        "target": result.target,
        "ranked_candidates": list(result.ranked_candidates),
        "scores": to_jsonable(result.scores),
        "assumptions": list(result.assumptions),
        "audit": to_jsonable(result.audit),
        "details": to_jsonable(result.details),
    }


def refinement_result_to_dict(result: Any) -> dict[str, object]:
    return {
        "winner": result.winner,
        "objective": result.objective,
        "scores": to_jsonable(result.scores),
        "certificates": [to_jsonable(certificate) for certificate in result.certificates],
        "audit": to_jsonable(result.audit),
        "details": to_jsonable(result.details),
    }


def to_jsonable(value: Any) -> Any:
    if hasattr(value, "to_dict"):
        return to_jsonable(value.to_dict())
    if is_dataclass(value):
        return {key: to_jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, np.ndarray):
        return [[float(entry) for entry in row] for row in np.asarray(value, dtype=float)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value
