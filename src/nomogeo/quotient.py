from __future__ import annotations

import numpy as np

from .core import visible_precision
from .types import GaussianContractionResult, LinearAlgebraMetadata
from .validation import (
    Tolerances,
    logdet_spd,
    resolve_tolerances,
    solve_spd,
    symmetrize,
    validate_spd_matrix,
    validate_surjective_map,
)
from .exceptions import InputValidationError


def observed_covariance(H: np.ndarray, C: np.ndarray, tolerances: Tolerances | None = None) -> np.ndarray:
    """Compute the observed covariance C H^{-1} C^T."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    c = validate_surjective_map(C, h.shape[0], tol)
    return symmetrize(c @ solve_spd(h, c.T))


def gaussian_forward_kl(cov_p: np.ndarray, cov_q: np.ndarray, tolerances: Tolerances | None = None) -> float:
    tol = resolve_tolerances(tolerances)
    p = validate_spd_matrix("cov_p", cov_p, tol)
    q = validate_spd_matrix("cov_q", cov_q, tol)
    if p.shape != q.shape:
        raise InputValidationError("covariance matrices must have the same shape")
    m = p.shape[0]
    return float(0.5 * (np.trace(solve_spd(q, p)) - m + logdet_spd(q) - logdet_spd(p)))


def gaussian_reverse_kl(cov_p: np.ndarray, cov_q: np.ndarray, tolerances: Tolerances | None = None) -> float:
    return gaussian_forward_kl(cov_q, cov_p, tolerances=tolerances)


def gaussian_bhattacharyya_distance(
    cov_p: np.ndarray,
    cov_q: np.ndarray,
    tolerances: Tolerances | None = None,
) -> float:
    tol = resolve_tolerances(tolerances)
    p = validate_spd_matrix("cov_p", cov_p, tol)
    q = validate_spd_matrix("cov_q", cov_q, tol)
    mean_cov = symmetrize(0.5 * (p + q))
    return float(0.5 * (logdet_spd(mean_cov) - 0.5 * (logdet_spd(p) + logdet_spd(q))))


def gaussian_hellinger_squared(
    cov_p: np.ndarray,
    cov_q: np.ndarray,
    tolerances: Tolerances | None = None,
) -> float:
    distance = gaussian_bhattacharyya_distance(cov_p, cov_q, tolerances=tolerances)
    value = 1.0 - np.exp(-distance)
    return float(max(0.0, min(1.0, value)))


def observer_collapse_descends(
    H1: np.ndarray,
    H2: np.ndarray,
    C1: np.ndarray,
    D: np.ndarray,
    tolerances: Tolerances | None = None,
) -> bool:
    """Check that equality for a richer observer descends to a coarsening."""
    tol = resolve_tolerances(tolerances)
    phi1 = visible_precision(H1, C1, tolerances=tol)
    phi2 = visible_precision(H2, C1, tolerances=tol)
    if not np.allclose(phi1, phi2, atol=tol.atol, rtol=tol.rtol):
        return False
    c0 = np.asarray(D, dtype=float) @ np.asarray(C1, dtype=float)
    coarse1 = visible_precision(H1, c0, tolerances=tol)
    coarse2 = visible_precision(H2, c0, tolerances=tol)
    return bool(np.allclose(coarse1, coarse2, atol=tol.atol, rtol=tol.rtol))


def gaussian_data_processing_contraction(
    H1: np.ndarray,
    H2: np.ndarray,
    C1: np.ndarray,
    D: np.ndarray,
    tolerances: Tolerances | None = None,
) -> GaussianContractionResult:
    """Compute fine/coarse Gaussian divergences under a coarsening C0 = D C1."""
    tol = resolve_tolerances(tolerances)
    c0 = np.asarray(D, dtype=float) @ np.asarray(C1, dtype=float)
    cov1_fine = observed_covariance(H1, C1, tolerances=tol)
    cov2_fine = observed_covariance(H2, C1, tolerances=tol)
    cov1_coarse = observed_covariance(H1, c0, tolerances=tol)
    cov2_coarse = observed_covariance(H2, c0, tolerances=tol)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="cholesky",
        ambient_dim=np.asarray(H1).shape[0],
        support_rank=np.asarray(C1).shape[0],
        visible_dim=np.asarray(C1).shape[0],
    )
    return GaussianContractionResult(
        forward_kl_fine=gaussian_forward_kl(cov1_fine, cov2_fine, tolerances=tol),
        forward_kl_coarse=gaussian_forward_kl(cov1_coarse, cov2_coarse, tolerances=tol),
        reverse_kl_fine=gaussian_reverse_kl(cov1_fine, cov2_fine, tolerances=tol),
        reverse_kl_coarse=gaussian_reverse_kl(cov1_coarse, cov2_coarse, tolerances=tol),
        hellinger_sq_fine=gaussian_hellinger_squared(cov1_fine, cov2_fine, tolerances=tol),
        hellinger_sq_coarse=gaussian_hellinger_squared(cov1_coarse, cov2_coarse, tolerances=tol),
        bhattacharyya_fine=gaussian_bhattacharyya_distance(cov1_fine, cov2_fine, tolerances=tol),
        bhattacharyya_coarse=gaussian_bhattacharyya_distance(cov1_coarse, cov2_coarse, tolerances=tol),
        metadata=metadata,
    )

