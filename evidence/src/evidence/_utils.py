from __future__ import annotations

from typing import Any

import numpy as np

from .exceptions import EvidenceInputError

EXTRACTION_MODES = {"exact_extraction", "encoded_inference", "open_ambiguity"}
EPISTEMIC_STATUSES = {"exact", "parsed", "inferred", "conjectural", "ambiguous"}


def require_choice(name: str, value: str, allowed: set[str]) -> str:
    if value not in allowed:
        raise EvidenceInputError(f"{name} must be one of {sorted(allowed)}")
    return value


def as_matrix(name: str, value: Any) -> np.ndarray:
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2:
        raise EvidenceInputError(f"{name} must be a 2D array")
    return matrix


def matrix_to_list(value: np.ndarray) -> list[list[float]]:
    array = np.asarray(value, dtype=float)
    return [[float(entry) for entry in row] for row in array]


def validate_epistemic_pair(mode: str, status: str, authoritative: bool) -> None:
    require_choice("extraction_mode", mode, EXTRACTION_MODES)
    require_choice("epistemic_status", status, EPISTEMIC_STATUSES)
    if mode == "exact_extraction" and status != "exact":
        raise EvidenceInputError("exact_extraction items must have epistemic_status='exact'")
    if mode == "open_ambiguity" and status != "ambiguous":
        raise EvidenceInputError("open_ambiguity items must have epistemic_status='ambiguous'")
    if mode == "open_ambiguity" and authoritative:
        raise EvidenceInputError("open_ambiguity items cannot be authoritative for downstream assembly")
