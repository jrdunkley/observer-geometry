from __future__ import annotations

from typing import Iterable

import numpy as np
import scipy.linalg as la


def as_array(value: np.ndarray | Iterable[Iterable[float]]) -> np.ndarray:
    return np.asarray(value, dtype=float)


def symmetric_basis(n: int) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        for j in range(i, n):
            matrix = np.zeros((n, n), dtype=float)
            matrix[i, j] = 1.0
            matrix[j, i] = 1.0 if i != j else 1.0
            basis.append(matrix)
    return basis


def vectorize_symmetric(matrix: np.ndarray) -> np.ndarray:
    n = matrix.shape[0]
    entries = []
    for i in range(n):
        for j in range(i, n):
            entries.append(float(matrix[i, j]))
    return np.array(entries, dtype=float)


def devectorize_symmetric(vector: np.ndarray, n: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i, n):
            matrix[i, j] = vector[idx]
            matrix[j, i] = vector[idx]
            idx += 1
    return matrix


def operator_constraint_rows(C: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    basis = symmetric_basis(n)
    rows = []
    coords = []
    for a in range(C.shape[0]):
        for b in range(C.shape[0]):
            row = np.array([float((C @ B @ C.T)[a, b]) for B in basis], dtype=float)
            rows.append(row)
            coords.append((a, b))
    return np.vstack(rows), np.array(coords, dtype=int)


def null_space(matrix: np.ndarray, rcond: float = 1e-12) -> np.ndarray:
    return la.null_space(matrix, rcond=rcond)


def min_eig(matrix: np.ndarray) -> float:
    if matrix.size == 0:
        return 0.0
    return float(np.min(np.linalg.eigvalsh(matrix)))
