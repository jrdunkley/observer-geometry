from __future__ import annotations

import numpy as np

from .types import DVBridgeResult, LinearAlgebraMetadata
from .validation import (
    Tolerances,
    inv_sqrt_spd,
    rank_cutoff,
    resolve_tolerances,
    symmetrize,
    validate_spd_matrix,
)
from .exceptions import InputValidationError


def dv_bridge(H0: np.ndarray, Jhat: np.ndarray, tolerances: Tolerances | None = None) -> DVBridgeResult:
    """Compute H_DV = H0 + 1/4 Jhat H0^{-1} Jhat^T and its Gram factor."""
    tol = resolve_tolerances(tolerances)
    h0 = validate_spd_matrix("H0", H0, tol)
    jhat = np.asarray(Jhat, dtype=float)
    if jhat.shape != h0.shape:
        raise InputValidationError("Jhat must match the shape of H0")
    if not np.allclose(jhat + jhat.T, 0.0, atol=tol.atol, rtol=tol.rtol):
        raise InputValidationError("Jhat must be antisymmetric")

    inv_sqrt = inv_sqrt_spd(h0, tol)
    gram_factor = 0.5 * jhat @ inv_sqrt
    delta = symmetrize(gram_factor @ gram_factor.T)
    h_dv = symmetrize(h0 + delta)
    eigs = np.linalg.eigvalsh(delta)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=rank_cutoff(eigs, tol),
        method="spectral",
        ambient_dim=h0.shape[0],
        support_rank=h0.shape[0],
    )
    return DVBridgeResult(h_dv=h_dv, delta_dv=delta, gram_factor=gram_factor, metadata=metadata)

