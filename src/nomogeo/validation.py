from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as la
from numpy.typing import NDArray

from .exceptions import InputValidationError, SupportError

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Tolerances:
    atol: float = 1e-10
    rtol: float = 1e-8


@dataclass(frozen=True)
class SupportDecomposition:
    basis: Array
    complement: Array
    eigenvalues: Array
    rank: int
    cutoff: float


def resolve_tolerances(tolerances: Tolerances | None = None) -> Tolerances:
    return tolerances if tolerances is not None else Tolerances()


def to_float_array(name: str, matrix: Array) -> Array:
    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2:
        raise InputValidationError(f"{name} must be a 2D array")
    return array


def symmetrize(matrix: Array) -> Array:
    return 0.5 * (matrix + matrix.T)


def square_matrix(name: str, matrix: Array) -> Array:
    array = to_float_array(name, matrix)
    if array.shape[0] != array.shape[1]:
        raise InputValidationError(f"{name} must be square")
    return array


def rank_cutoff(values: Array, tolerances: Tolerances) -> float:
    scale = float(np.max(np.abs(values))) if values.size else 0.0
    return max(tolerances.atol, tolerances.rtol * max(1.0, scale))


def validate_symmetric_matrix(name: str, matrix: Array, tolerances: Tolerances) -> Array:
    array = square_matrix(name, matrix)
    if not np.allclose(array, array.T, atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError(f"{name} must be symmetric")
    return symmetrize(array)


def validate_spd_matrix(name: str, matrix: Array, tolerances: Tolerances) -> Array:
    array = validate_symmetric_matrix(name, matrix, tolerances)
    eigenvalues = np.linalg.eigvalsh(array)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    if float(np.min(eigenvalues)) <= cutoff:
        raise InputValidationError(f"{name} must be strictly positive definite")
    return array


def validate_psd_matrix(name: str, matrix: Array, tolerances: Tolerances) -> Array:
    array = validate_symmetric_matrix(name, matrix, tolerances)
    eigenvalues = np.linalg.eigvalsh(array)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    if float(np.min(eigenvalues)) < -cutoff:
        raise InputValidationError(f"{name} must be positive semidefinite")
    return array


def validate_surjective_map(C: Array, latent_dim: int, tolerances: Tolerances) -> Array:
    matrix = to_float_array("C", C)
    if matrix.shape[1] != latent_dim:
        raise InputValidationError("C has incompatible latent dimension")
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    cutoff = rank_cutoff(singular_values, tolerances)
    rank = int(np.sum(singular_values > cutoff))
    if rank != matrix.shape[0]:
        raise InputValidationError("C must be surjective")
    return matrix


def support_decomposition_psd(matrix: Array, tolerances: Tolerances) -> SupportDecomposition:
    psd = validate_psd_matrix("matrix", matrix, tolerances)
    eigenvalues, eigenvectors = la.eigh(psd)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    keep = eigenvalues > cutoff
    basis = eigenvectors[:, keep]
    complement = eigenvectors[:, ~keep]
    return SupportDecomposition(
        basis=basis,
        complement=complement,
        eigenvalues=eigenvalues[keep],
        rank=int(np.sum(keep)),
        cutoff=cutoff,
    )


def range_basis(matrix: Array, tolerances: Tolerances) -> Array:
    array = to_float_array("matrix", matrix)
    u, singular_values, _ = la.svd(array, full_matrices=True)
    cutoff = rank_cutoff(singular_values, tolerances)
    rank = int(np.sum(singular_values > cutoff))
    return u[:, :rank]


def restrict(matrix: Array, basis: Array) -> Array:
    if basis.size == 0:
        return np.zeros((0, 0), dtype=float)
    return symmetrize(basis.T @ matrix @ basis)


def embed(reduced: Array, basis: Array, ambient_dim: int) -> Array:
    if reduced.size == 0:
        return np.zeros((ambient_dim, ambient_dim), dtype=float)
    return symmetrize(basis @ reduced @ basis.T)


def solve_spd(matrix: Array, rhs: Array) -> Array:
    factor = la.cho_factor(matrix, lower=True, check_finite=True)
    return np.asarray(la.cho_solve(factor, rhs, check_finite=True), dtype=float)


def inverse_spd(matrix: Array) -> Array:
    return symmetrize(solve_spd(matrix, np.eye(matrix.shape[0], dtype=float)))


def logdet_spd(matrix: Array) -> float:
    sign, value = np.linalg.slogdet(matrix)
    if sign <= 0:
        raise InputValidationError("matrix must be positive definite to compute logdet")
    return float(value)


def sqrt_psd(matrix: Array, tolerances: Tolerances) -> Array:
    psd = validate_psd_matrix("matrix", matrix, tolerances)
    eigenvalues, eigenvectors = la.eigh(psd)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    eigenvalues = np.where(eigenvalues > cutoff, eigenvalues, 0.0)
    return symmetrize((eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T)


def inv_sqrt_spd(matrix: Array, tolerances: Tolerances) -> Array:
    spd = validate_spd_matrix("matrix", matrix, tolerances)
    eigenvalues, eigenvectors = la.eigh(spd)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    if np.min(eigenvalues) <= cutoff:
        raise InputValidationError("matrix must be positive definite")
    return symmetrize((eigenvectors * (1.0 / np.sqrt(eigenvalues))) @ eigenvectors.T)


def spectral_factor_psd(matrix: Array, tolerances: Tolerances) -> Array:
    psd = validate_psd_matrix("matrix", matrix, tolerances)
    eigenvalues, eigenvectors = la.eigh(psd)
    cutoff = rank_cutoff(eigenvalues, tolerances)
    keep = eigenvalues > cutoff
    if not np.any(keep):
        return np.zeros((matrix.shape[0], 0), dtype=float)
    return eigenvectors[:, keep] * np.sqrt(eigenvalues[keep])


def assert_support_preserving(
    X: Array,
    support_basis: Array,
    complement: Array,
    tolerances: Tolerances,
) -> None:
    if complement.size == 0:
        return
    cross = support_basis.T @ X @ complement
    null_block = complement.T @ X @ complement
    residual = max(float(np.linalg.norm(cross)), float(np.linalg.norm(null_block)))
    if residual > 10.0 * max(tolerances.atol, tolerances.rtol):
        raise SupportError("X does not preserve the support of the ceiling")

