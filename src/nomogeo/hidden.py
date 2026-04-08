from __future__ import annotations

import numpy as np

from .exceptions import InputValidationError, SupportError
from .types import HiddenLoadResult, LinearAlgebraMetadata, RealisationResult
from .validation import (
    Tolerances,
    assert_support_preserving,
    embed,
    inv_sqrt_spd,
    logdet_spd,
    rank_cutoff,
    resolve_tolerances,
    restrict,
    solve_spd,
    spectral_factor_psd,
    sqrt_psd,
    support_decomposition_psd,
    symmetrize,
    validate_psd_matrix,
    validate_symmetric_matrix,
)


def _resolve_support(
    T: np.ndarray,
    X: np.ndarray,
    support_mode: str,
    tolerances: Tolerances,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool]:
    ceiling = validate_psd_matrix("T", T, tolerances)
    x = validate_symmetric_matrix("X", X, tolerances)
    if ceiling.shape != x.shape:
        raise InputValidationError("T and X must have the same shape")

    decomposition = support_decomposition_psd(ceiling, tolerances)
    if support_mode not in {"auto", "support", "ambient"}:
        raise InputValidationError("support_mode must be 'auto', 'support', or 'ambient'")
    if support_mode == "ambient":
        if decomposition.rank != ceiling.shape[0]:
            raise SupportError("ambient mode requires a strictly positive definite ceiling")
        basis = np.eye(ceiling.shape[0], dtype=float)
        complement = np.zeros((ceiling.shape[0], 0), dtype=float)
        support_restricted = False
    else:
        basis = decomposition.basis
        complement = decomposition.complement
        support_restricted = decomposition.rank != ceiling.shape[0]

    assert_support_preserving(x, basis, complement, tolerances)
    return ceiling, x, basis, complement, support_restricted


def hidden_load(
    T: np.ndarray,
    X: np.ndarray,
    support_mode: str = "auto",
    tolerances: Tolerances | None = None,
) -> HiddenLoadResult:
    """Compute the hidden load Lambda beneath the ceiling T."""
    tol = resolve_tolerances(tolerances)
    ceiling, x, basis, _complement, support_restricted = _resolve_support(T, X, support_mode, tol)
    t_s = restrict(ceiling, basis)
    x_s = restrict(x, basis)

    if basis.shape[1] == 0:
        reduced_lambda = np.zeros((0, 0), dtype=float)
        lambda_ambient = np.zeros_like(ceiling)
        metadata = LinearAlgebraMetadata(
            atol=tol.atol,
            rtol=tol.rtol,
            rank_tol=max(tol.atol, tol.rtol),
            method="eigh",
            ambient_dim=ceiling.shape[0],
            support_rank=0,
            support_restricted=True,
            notes=("zero-support ceiling",),
        )
        return HiddenLoadResult(
            lambda_=lambda_ambient,
            reduced_lambda=reduced_lambda,
            X=np.zeros_like(x),
            ceiling=ceiling,
            support_basis=basis,
            rank=0,
            clock=0.0,
            metadata=metadata,
        )

    x_s = validate_symmetric_matrix("X|S", x_s, tol)
    min_x = float(np.min(np.linalg.eigvalsh(x_s)))
    cutoff = rank_cutoff(np.linalg.eigvalsh(x_s), tol)
    if min_x <= cutoff:
        raise SupportError("X|S must be strictly positive definite on the active support")

    gap = symmetrize(t_s - x_s)
    gap_eigs = np.linalg.eigvalsh(gap)
    if float(np.min(gap_eigs)) < -rank_cutoff(gap_eigs, tol):
        raise SupportError("X must lie below T on the active support")

    sqrt_t = sqrt_psd(t_s, tol)
    reduced_lambda = symmetrize(sqrt_t @ solve_spd(x_s, sqrt_t) - np.eye(t_s.shape[0], dtype=float))
    lambda_eigs = np.linalg.eigvalsh(reduced_lambda)
    lambda_cutoff = rank_cutoff(lambda_eigs, tol)
    if float(np.min(lambda_eigs)) < -lambda_cutoff:
        raise SupportError("computed hidden load is not positive semidefinite")

    lambda_ambient = embed(reduced_lambda, basis, ceiling.shape[0])
    positive_rank = int(np.sum(lambda_eigs > lambda_cutoff))
    clock_value = logdet_spd(np.eye(reduced_lambda.shape[0], dtype=float) + reduced_lambda)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=lambda_cutoff,
        method="eigh",
        ambient_dim=ceiling.shape[0],
        support_rank=basis.shape[1],
        support_restricted=support_restricted,
    )
    return HiddenLoadResult(
        lambda_=lambda_ambient,
        reduced_lambda=reduced_lambda,
        X=x,
        ceiling=ceiling,
        support_basis=basis,
        rank=positive_rank,
        clock=clock_value,
        metadata=metadata,
    )


def visible_from_hidden_load(
    T: np.ndarray,
    Lambda: np.ndarray,
    support_mode: str = "auto",
    lambda_representation: str | None = None,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Reconstruct X = T^{1/2} (I + Lambda)^{-1} T^{1/2} on the active support."""
    tol = resolve_tolerances(tolerances)
    ceiling = validate_psd_matrix("T", T, tol)
    decomposition = support_decomposition_psd(ceiling, tol)
    if support_mode not in {"auto", "support", "ambient"}:
        raise InputValidationError("support_mode must be 'auto', 'support', or 'ambient'")
    if support_mode == "ambient" and decomposition.rank != ceiling.shape[0]:
        raise SupportError("ambient mode requires a strictly positive definite ceiling")
    basis = np.eye(ceiling.shape[0], dtype=float) if support_mode == "ambient" else decomposition.basis
    t_s = restrict(ceiling, basis)
    support_dim = basis.shape[1]

    if support_dim == 0:
        lambda_array = validate_symmetric_matrix("Lambda", Lambda, tol)
        if lambda_array.shape not in {(0, 0), ceiling.shape}:
            raise InputValidationError("Lambda has incompatible shape for a zero-support ceiling")
        if lambda_array.size != 0 and np.linalg.norm(lambda_array) > max(tol.atol, tol.rtol):
            raise InputValidationError("zero-support ceilings only admit the zero hidden load")
        return np.zeros_like(ceiling)

    lambda_array = validate_symmetric_matrix("Lambda", Lambda, tol)
    reduced_shape = (support_dim, support_dim)
    ambient_shape = ceiling.shape
    ambiguous = reduced_shape == ambient_shape

    if lambda_representation is not None and lambda_representation not in {"reduced", "ambient"}:
        raise InputValidationError("lambda_representation must be 'reduced' or 'ambient'")
    if lambda_representation is None:
        if ambiguous and lambda_array.shape == ambient_shape:
            raise InputValidationError(
                "Lambda representation is ambiguous for full-rank ceilings; "
                "set lambda_representation='reduced' or 'ambient'"
            )
        if lambda_array.shape == ambient_shape and not ambiguous:
            lambda_representation = "ambient"
        elif lambda_array.shape == reduced_shape:
            lambda_representation = "reduced"
        else:
            raise InputValidationError("Lambda has incompatible shape")

    if lambda_representation == "ambient":
        if lambda_array.shape != ambient_shape:
            raise InputValidationError("ambient Lambda must have the same shape as T")
        lambda_s = restrict(lambda_array, basis)
    else:
        if lambda_array.shape != reduced_shape:
            raise InputValidationError("reduced Lambda must match the active-support dimension")
        lambda_s = lambda_array

    if lambda_s.size == 0:
        return np.zeros_like(ceiling)

    lambda_eigs = np.linalg.eigvalsh(lambda_s)
    if float(np.min(lambda_eigs)) < -rank_cutoff(lambda_eigs, tol):
        raise InputValidationError("Lambda must be positive semidefinite on the active support")

    if lambda_representation == "ambient" and decomposition.complement.size != 0:
        cross = basis.T @ lambda_array @ decomposition.complement
        null_block = decomposition.complement.T @ lambda_array @ decomposition.complement
        residual = max(float(np.linalg.norm(cross)), float(np.linalg.norm(null_block)))
        if residual > 10.0 * max(tol.atol, tol.rtol):
            raise InputValidationError("ambient Lambda must preserve the active support")

    sqrt_t = sqrt_psd(t_s, tol)
    pi = solve_spd(np.eye(lambda_s.shape[0], dtype=float) + lambda_s, np.eye(lambda_s.shape[0], dtype=float))
    x_s = symmetrize(sqrt_t @ pi @ sqrt_t)
    return embed(x_s, basis, ceiling.shape[0])


def inverse_visible_class(
    T: np.ndarray,
    Lambda: np.ndarray,
    support_mode: str = "auto",
    lambda_representation: str | None = None,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Theorem-facing alias for the fixed-ceiling inverse map Lambda -> X beneath T."""
    return visible_from_hidden_load(
        T,
        Lambda,
        support_mode=support_mode,
        lambda_representation=lambda_representation,
        tolerances=tolerances,
    )


def hidden_contraction(Lambda: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Return the canonical contraction factor K = (I + Lambda)^(-1/2)."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_symmetric_matrix("Lambda", Lambda, tol)
    if lambda_array.shape == (0, 0):
        return np.zeros((0, 0), dtype=float)
    lambda_eigs = np.linalg.eigvalsh(lambda_array)
    lambda_cutoff = rank_cutoff(lambda_eigs, tol)
    if float(np.min(lambda_eigs)) < -lambda_cutoff:
        raise InputValidationError("Lambda must be positive semidefinite on the active support")
    shifted = np.eye(lambda_array.shape[0], dtype=float) + lambda_array
    return inv_sqrt_spd(shifted, tol)


def load_from_hidden_contraction(K: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Recover Lambda from a contraction factor K with K^T K = (I + Lambda)^(-1)."""
    tol = resolve_tolerances(tolerances)
    factor = np.asarray(K, dtype=float)
    if factor.ndim != 2 or factor.shape[0] != factor.shape[1]:
        raise InputValidationError("K must be a square 2D array")
    if factor.shape == (0, 0):
        return np.zeros((0, 0), dtype=float)

    pi = symmetrize(factor.T @ factor)
    pi_eigs = np.linalg.eigvalsh(pi)
    pi_cutoff = rank_cutoff(pi_eigs, tol)
    if float(np.min(pi_eigs)) <= pi_cutoff:
        raise InputValidationError("K must have full rank on the active support")
    if float(np.max(pi_eigs)) > 1.0 + 10.0 * pi_cutoff:
        raise InputValidationError("K must be a contraction on the active support")

    lambda_array = symmetrize(solve_spd(pi, np.eye(pi.shape[0], dtype=float)) - np.eye(pi.shape[0], dtype=float))
    lambda_eigs = np.linalg.eigvalsh(lambda_array)
    lambda_cutoff = rank_cutoff(lambda_eigs, tol)
    if float(np.min(lambda_eigs)) < -lambda_cutoff:
        raise InputValidationError("K does not define a positive hidden load")
    return lambda_array


def canonical_hidden_realisation(
    T: np.ndarray,
    X: np.ndarray,
    support_mode: str = "auto",
    tolerances: Tolerances | None = None,
) -> RealisationResult:
    """Construct the canonical hidden realisation on the active support."""
    load = hidden_load(T, X, support_mode=support_mode, tolerances=tolerances)
    tol = resolve_tolerances(tolerances)
    t_s = restrict(load.ceiling, load.support_basis)
    sqrt_t = sqrt_psd(t_s, tol)
    sqrt_lambda = sqrt_psd(load.reduced_lambda, tol)

    upper = np.hstack([t_s, sqrt_t @ sqrt_lambda])
    lower = np.hstack([sqrt_lambda @ sqrt_t, np.eye(load.reduced_lambda.shape[0], dtype=float) + load.reduced_lambda])
    matrix = np.vstack([upper, lower])
    return RealisationResult(
        matrix=symmetrize(matrix),
        factor=sqrt_lambda,
        support_basis=load.support_basis,
        rank=load.rank,
        metadata=load.metadata,
    )


def minimal_hidden_realisation(
    T: np.ndarray,
    X: np.ndarray,
    support_mode: str = "auto",
    tolerances: Tolerances | None = None,
) -> RealisationResult:
    """Construct a minimal hidden realisation using a Gram factor of Lambda."""
    load = hidden_load(T, X, support_mode=support_mode, tolerances=tolerances)
    tol = resolve_tolerances(tolerances)
    t_s = restrict(load.ceiling, load.support_basis)
    sqrt_t = sqrt_psd(t_s, tol)
    factor = spectral_factor_psd(load.reduced_lambda, tol)

    upper = np.hstack([t_s, sqrt_t @ factor])
    lower = np.hstack([factor.T @ sqrt_t, np.eye(factor.shape[1], dtype=float) + factor.T @ factor])
    matrix = np.vstack([upper, lower])
    return RealisationResult(
        matrix=symmetrize(matrix),
        factor=factor,
        support_basis=load.support_basis,
        rank=load.rank,
        metadata=load.metadata,
    )


def transport_hidden_load(Lambda: np.ndarray, M: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Transport hidden load under sequential composition."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_symmetric_matrix("Lambda", Lambda, tol)
    m_array = validate_symmetric_matrix("M", M, tol)
    if lambda_array.shape != m_array.shape:
        raise InputValidationError("Lambda and M must have the same shape")
    if lambda_array.shape == (0, 0):
        return np.zeros((0, 0), dtype=float)

    lambda_eigs = np.linalg.eigvalsh(lambda_array)
    lambda_cutoff = rank_cutoff(lambda_eigs, tol)
    if float(np.min(lambda_eigs)) < -lambda_cutoff:
        raise InputValidationError("Lambda must be positive semidefinite on the active support")

    m_eigs = np.linalg.eigvalsh(m_array)
    m_cutoff = rank_cutoff(m_eigs, tol)
    if float(np.min(m_eigs)) < -m_cutoff:
        raise InputValidationError("M must be positive semidefinite on the active support")

    root = sqrt_psd(np.eye(lambda_array.shape[0], dtype=float) + lambda_array, tol)
    total = lambda_array + root @ m_array @ root
    return symmetrize(total)


def clock(Lambda: np.ndarray, tolerances: Tolerances | None = None) -> float:
    """Compute tau(Lambda) = log det(I + Lambda)."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_symmetric_matrix("Lambda", Lambda, tol)
    if lambda_array.shape == (0, 0):
        return 0.0
    lambda_eigs = np.linalg.eigvalsh(lambda_array)
    lambda_cutoff = rank_cutoff(lambda_eigs, tol)
    if float(np.min(lambda_eigs)) < -lambda_cutoff:
        raise InputValidationError("Lambda must be positive semidefinite on the active support")
    shifted = np.eye(lambda_array.shape[0], dtype=float) + lambda_array
    return logdet_spd(shifted)
