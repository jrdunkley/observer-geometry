from __future__ import annotations

import numpy as np
import scipy.linalg as la


def random_spd(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a + n * np.eye(n, dtype=float)


def random_surjective(rng: np.random.Generator, m: int, n: int) -> np.ndarray:
    while True:
        c = rng.normal(size=(m, n))
        if np.linalg.matrix_rank(c) == m:
            return c


def random_skew(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a - a.T


def random_psd_on_subspace(rng: np.random.Generator, n: int, rank: int) -> tuple[np.ndarray, np.ndarray]:
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    basis = q[:, :rank]
    values = np.linspace(1.1, 2.3, rank)
    matrix = basis @ np.diag(values) @ basis.T
    return matrix, basis


def orthogonal_matrix(rng: np.random.Generator, n: int) -> np.ndarray:
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return q


def random_invertible(rng: np.random.Generator, n: int) -> np.ndarray:
    while True:
        matrix = rng.normal(size=(n, n))
        if abs(np.linalg.det(matrix)) > 1e-2:
            return matrix


def random_psd(rng: np.random.Generator, n: int, rank: int | None = None) -> np.ndarray:
    target_rank = n if rank is None else rank
    factor = rng.normal(size=(n, target_rank))
    return factor @ factor.T


def null_space_basis(matrix: np.ndarray, rcond: float = 1e-12) -> np.ndarray:
    return la.null_space(matrix, rcond=rcond)


def schur_complement(block_matrix: np.ndarray, visible_dim: int) -> np.ndarray:
    vv = block_matrix[:visible_dim, :visible_dim]
    vh = block_matrix[:visible_dim, visible_dim:]
    hv = block_matrix[visible_dim:, :visible_dim]
    hh = block_matrix[visible_dim:, visible_dim:]
    if hh.size == 0:
        return vv
    return vv - vh @ np.linalg.solve(hh, hv)


def square_value(x: int) -> int:
    return x * x


def add_pair(x: int, y: int) -> int:
    return x + y


def fail_on_negative(x: int) -> int:
    if x < 0:
        raise ValueError("negative input")
    return x * x


def visible_precision_trace_task(H: np.ndarray, C: np.ndarray) -> float:
    from nomogeo import visible_precision

    return float(np.trace(visible_precision(H, C)))
