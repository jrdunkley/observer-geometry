from __future__ import annotations

import numpy as np

from .exceptions import InputValidationError
from .types import (
    LinearAlgebraMetadata,
    RankKCovariancePerturbationResult,
    RankOneCovariancePerturbationResult,
    ResidualMarginResult,
)
from .validation import (
    Tolerances,
    inverse_spd,
    rank_cutoff,
    resolve_tolerances,
    symmetrize,
    validate_spd_matrix,
)


def rank_one_covariance_perturbation(
    covariance_base: np.ndarray,
    signal: np.ndarray,
    visible_dim: int,
    epsilon: float,
    tolerances: Tolerances | None = None,
) -> RankOneCovariancePerturbationResult:
    """Diagnose a rank-one covariance/Fisher perturbation under a coordinate split.

    The theorem domain is ``Sigma_epsilon = Sigma_0 + epsilon f f^T`` with a
    visible coordinate block of size ``visible_dim``. In a white/aligned
    background the hidden-gap increment is one-channel. In a generic coloured
    background it is the difference of two rank-one terms and has rank at most
    two. This is a covariance/Fisher diagnostic, not a generic non-Gaussian
    full-law selector.
    """
    if epsilon < 0.0:
        raise InputValidationError("epsilon must be nonnegative")
    factor = np.sqrt(float(epsilon)) * np.asarray(signal, dtype=float)
    rank_k = rank_k_covariance_perturbation(covariance_base, factor[:, None], visible_dim, tolerances)
    tol = resolve_tolerances(tolerances)
    sigma0 = rank_k.covariance_base
    f = np.asarray(signal, dtype=float)
    if f.ndim != 1 or f.shape[0] != sigma0.shape[0]:
        raise InputValidationError("signal must be a 1D vector with length matching covariance_base")
    n = sigma0.shape[0]
    if visible_dim <= 0 or visible_dim >= n:
        raise InputValidationError("visible_dim must be strictly between 0 and the ambient dimension")

    m = int(visible_dim)
    phi0 = inverse_spd(validate_spd_matrix("covariance_base visible block", sigma0[:m, :m], tol))
    f_v = f[:m]
    u_v = rank_k.precision_base[:m, :] @ f
    w_v = phi0 @ f_v
    beta_denom = 1.0 + epsilon * float(f @ rank_k.precision_base @ f)
    gamma_denom = 1.0 + epsilon * float(f_v @ phi0 @ f_v)
    if beta_denom <= 0.0 or gamma_denom <= 0.0:
        raise InputValidationError("rank-one covariance update left theorem domain")

    denom = float(np.linalg.norm(u_v) * np.linalg.norm(w_v))
    if denom <= max(tol.atol, tol.rtol):
        direction_alignment = 1.0
    else:
        direction_alignment = abs(float(u_v @ w_v)) / denom
        direction_alignment = min(1.0, max(0.0, direction_alignment))
    one_channel = bool(rank_k.update_rank <= 1)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=rank_k.metadata.rank_tol,
        method="rank-one-covariance-sherman-morrison",
        ambient_dim=n,
        support_rank=rank_k.update_rank,
        visible_dim=m,
        condition_number=float(np.linalg.cond(sigma0)),
        notes=(
            "coordinate visible split",
            "covariance/Fisher diagnostic, not a generic non-Gaussian full-law selector",
        ),
    )
    return RankOneCovariancePerturbationResult(
        covariance_base=sigma0,
        covariance_perturbed=rank_k.covariance_perturbed,
        precision_base=rank_k.precision_base,
        precision_perturbed=rank_k.precision_perturbed,
        hidden_gap_increment=rank_k.hidden_gap_increment,
        formula_increment=rank_k.formula_increment,
        formula_residual=rank_k.formula_residual,
        full_precision_visible_direction=u_v,
        visible_precision_direction=w_v,
        direction_alignment=direction_alignment,
        singular_values=rank_k.singular_values,
        update_rank=rank_k.update_rank,
        one_channel=one_channel,
        metadata=metadata,
    )


def rank_k_covariance_perturbation(
    covariance_base: np.ndarray,
    perturbation_factor: np.ndarray,
    visible_dim: int,
    tolerances: Tolerances | None = None,
) -> RankKCovariancePerturbationResult:
    """Diagnose a rank-k covariance/Fisher perturbation under a coordinate split.

    The theorem domain is ``Sigma_1 = Sigma_0 + F F^T`` with a coordinate
    visible block of size ``visible_dim``. The hidden-gap increment is a
    difference of two rank-k positive semidefinite terms, hence has rank at
    most ``2k``. This is finite covariance/Fisher algebra, not a generic
    non-Gaussian law selector.
    """
    tol = resolve_tolerances(tolerances)
    sigma0 = validate_spd_matrix("covariance_base", covariance_base, tol)
    factor = np.asarray(perturbation_factor, dtype=float)
    if factor.ndim == 1:
        factor = factor[:, None]
    if factor.ndim != 2 or factor.shape[0] != sigma0.shape[0] or factor.shape[1] < 1:
        raise InputValidationError("perturbation_factor must have shape (ambient_dim, k)")
    if not np.all(np.isfinite(factor)):
        raise InputValidationError("perturbation_factor must be finite")
    n = sigma0.shape[0]
    if visible_dim <= 0 or visible_dim >= n:
        raise InputValidationError("visible_dim must be strictly between 0 and the ambient dimension")

    m = int(visible_dim)
    k = factor.shape[1]
    sigma1 = symmetrize(sigma0 + factor @ factor.T)
    h0 = inverse_spd(sigma0)
    h1 = inverse_spd(sigma1)
    sigma0_v = validate_spd_matrix("covariance_base visible block", sigma0[:m, :m], tol)
    sigma1_v = validate_spd_matrix("covariance_perturbed visible block", sigma1[:m, :m], tol)
    phi0 = inverse_spd(sigma0_v)
    phi1 = inverse_spd(sigma1_v)

    gap0 = symmetrize(h0[:m, :m] - phi0)
    gap1 = symmetrize(h1[:m, :m] - phi1)
    direct = symmetrize(gap1 - gap0)

    f_v = factor[:m, :]
    u_v = h0[:m, :] @ factor
    w_v = phi0 @ f_v
    beta = inverse_spd(np.eye(k) + factor.T @ h0 @ factor)
    gamma = inverse_spd(np.eye(k) + f_v.T @ phi0 @ f_v)
    full_term = symmetrize(u_v @ beta @ u_v.T)
    visible_term = symmetrize(w_v @ gamma @ w_v.T)
    formula = symmetrize(visible_term - full_term)

    singular_values = np.linalg.svd(formula, compute_uv=False)
    cutoff = rank_cutoff(singular_values, tol)
    update_rank = int(np.sum(singular_values > cutoff))
    residual = float(np.max(np.abs(direct - formula))) if direct.size else 0.0
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="rank-k-covariance-woodbury",
        ambient_dim=n,
        support_rank=update_rank,
        visible_dim=m,
        condition_number=float(np.linalg.cond(sigma0)),
        notes=(
            "coordinate visible split",
            "hidden-gap increment is a difference of two rank-k PSD terms",
            "covariance/Fisher diagnostic, not a generic non-Gaussian full-law selector",
        ),
    )
    return RankKCovariancePerturbationResult(
        covariance_base=sigma0,
        covariance_perturbed=sigma1,
        precision_base=h0,
        precision_perturbed=h1,
        perturbation_factor=factor.copy(),
        hidden_gap_increment=direct,
        formula_increment=formula,
        formula_residual=residual,
        full_precision_visible_factor=u_v,
        visible_precision_factor=w_v,
        full_precision_term=full_term,
        visible_precision_term=visible_term,
        singular_values=singular_values,
        update_rank=update_rank,
        rank_bound=2 * k,
        metadata=metadata,
    )


def residual_margin_ordering(
    quadratic_gap: float,
    residual_bound: float,
    tolerances: Tolerances | None = None,
) -> ResidualMarginResult:
    """Certify whether a positive quadratic branch gap survives residuals.

    If two residual terms are bounded by ``residual_bound`` in absolute value,
    a quadratic ordering gap ``quadratic_gap`` is forced only when
    ``quadratic_gap > 2 * residual_bound``. Equality is not a strict certificate.
    """
    tol = resolve_tolerances(tolerances)
    gap = float(quadratic_gap)
    bound = float(residual_bound)
    if gap < 0.0:
        raise InputValidationError("quadratic_gap must be nonnegative")
    if bound < 0.0:
        raise InputValidationError("residual_bound must be nonnegative")
    required_gap = 2.0 * bound
    margin = gap - required_gap
    robust = bool(margin > max(tol.atol, tol.rtol * max(1.0, gap, required_gap)))
    worst_case_gap = margin
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="residual-margin-ordering",
        ambient_dim=2,
        support_rank=0,
        visible_dim=None,
        notes=("strict certificate requires quadratic_gap > 2 * residual_bound",),
    )
    return ResidualMarginResult(
        quadratic_gap=gap,
        residual_bound=bound,
        required_gap=required_gap,
        margin=margin,
        worst_case_gap=worst_case_gap,
        robust=robust,
        adversarial_reversal_possible=not robust,
        metadata=metadata,
    )
