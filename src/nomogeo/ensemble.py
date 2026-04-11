from __future__ import annotations

import numpy as np

from .core import visible_precision
from .exceptions import InputValidationError
from .hidden import hidden_load
from .types import (
    CeilingMediatedLocalQuadraticEnsembleResult,
    CoordinateLocalQuadraticEnsembleResult,
    IntrinsicLocalQuadraticEnsembleResult,
    LinearAlgebraMetadata,
)
from .validation import (
    Tolerances,
    logdet_spd,
    rank_cutoff,
    resolve_tolerances,
    validate_spd_matrix,
    validate_surjective_map,
)


def intrinsic_local_quadratic_ensemble(
    hessians: np.ndarray,
    C: np.ndarray,
    tolerances: Tolerances | None = None,
) -> IntrinsicLocalQuadraticEnsembleResult:
    """Evaluate the intrinsic quotient ensemble ``Phi_C(H_i)``.

    This mode is observer-intrinsic: it returns only visible precision samples
    and functorial log-determinant summaries. Hidden load and clock require a
    separately declared ceiling family and live in
    ``ceiling_mediated_local_quadratic_ensemble``.
    """
    tol = resolve_tolerances(tolerances)
    samples = _validate_hessian_stack(hessians, tol)
    sample_count, n, _ = samples.shape
    c = validate_surjective_map(C, n, tol)
    m = c.shape[0]

    phis = np.zeros((sample_count, m, m), dtype=float)
    logdets = np.zeros(sample_count, dtype=float)
    condition_numbers = np.zeros(sample_count, dtype=float)
    for idx in range(sample_count):
        h = validate_spd_matrix(f"hessians[{idx}]", samples[idx], tol)
        phi = visible_precision(h, c, tolerances=tol)
        phis[idx] = phi
        logdets[idx] = logdet_spd(phi)
        condition_numbers[idx] = float(np.linalg.cond(h))

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=rank_cutoff(logdets, tol),
        method="intrinsic-local-quadratic-ensemble",
        ambient_dim=n,
        support_rank=sample_count,
        visible_dim=m,
        condition_number=float(np.max(condition_numbers)),
        support_restricted=False,
        notes=(
            "intrinsic quotient ensemble: Phi_C(H_i) only",
            "hidden load and clock require explicit ceilings",
        ),
    )
    return IntrinsicLocalQuadraticEnsembleResult(
        phis=phis,
        logdet_phis=logdets,
        mean_logdet_phi=float(np.mean(logdets)),
        std_logdet_phi=float(np.std(logdets)),
        min_logdet_phi=float(np.min(logdets)),
        max_logdet_phi=float(np.max(logdets)),
        sample_count=sample_count,
        metadata=metadata,
    )


def ceiling_mediated_local_quadratic_ensemble(
    hessians: np.ndarray,
    C: np.ndarray,
    ceilings: np.ndarray,
    tolerances: Tolerances | None = None,
) -> CeilingMediatedLocalQuadraticEnsembleResult:
    """Evaluate hidden load and clock for an ensemble with declared ceilings."""
    tol = resolve_tolerances(tolerances)
    intrinsic = intrinsic_local_quadratic_ensemble(hessians, C, tolerances=tol)
    ceiling_stack = _validate_ceiling_stack(ceilings, intrinsic.sample_count, intrinsic.phis.shape[1], tol)
    sample_count, m, _ = ceiling_stack.shape

    lambdas = np.zeros((sample_count, m, m), dtype=float)
    clocks = np.zeros(sample_count, dtype=float)
    ranks = np.zeros(sample_count, dtype=float)
    ceiling_conditions = np.zeros(sample_count, dtype=float)
    for idx in range(sample_count):
        ceiling = validate_spd_matrix(f"ceilings[{idx}]", ceiling_stack[idx], tol)
        load = hidden_load(ceiling, intrinsic.phis[idx], support_mode="ambient", tolerances=tol)
        lambdas[idx] = load.lambda_
        clocks[idx] = load.clock
        ranks[idx] = float(load.rank)
        ceiling_conditions[idx] = float(np.linalg.cond(ceiling))

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=rank_cutoff(clocks, tol),
        method="ceiling-mediated-local-quadratic-ensemble",
        ambient_dim=np.asarray(hessians).shape[1],
        support_rank=sample_count,
        visible_dim=m,
        condition_number=float(np.max(ceiling_conditions)),
        support_restricted=False,
        notes=(
            "ceiling-mediated ensemble with explicitly supplied ceilings",
            "clock summaries are descriptive ensemble statistics, not full-law cumulants",
        ),
    )
    return CeilingMediatedLocalQuadraticEnsembleResult(
        phis=intrinsic.phis,
        ceilings=ceiling_stack,
        lambdas=lambdas,
        clocks=clocks,
        hidden_ranks=ranks,
        mean_clock=float(np.mean(clocks)),
        std_clock=float(np.std(clocks)),
        min_clock=float(np.min(clocks)),
        max_clock=float(np.max(clocks)),
        sample_count=sample_count,
        metadata=metadata,
    )


def coordinate_local_quadratic_ensemble(
    hessians: np.ndarray,
    visible_dim: int,
    tolerances: Tolerances | None = None,
) -> CoordinateLocalQuadraticEnsembleResult:
    """Evaluate exact local quadratic observer geometry samplewise.

    This is a coordinate-split ensemble helper. Each sample is an SPD local
    Hessian/Fisher object; the visible observer is the first ``visible_dim``
    coordinates. The returned summaries are descriptive statistics over exact
    samplewise quadratic outputs, not full non-Gaussian law cumulants.
    """
    tol = resolve_tolerances(tolerances)
    samples = _validate_hessian_stack(hessians, tol)
    sample_count, n, _ = samples.shape
    if visible_dim <= 0 or visible_dim >= n:
        raise InputValidationError("visible_dim must be strictly between 0 and the ambient dimension")
    m = int(visible_dim)
    C = np.eye(m, n, dtype=float)
    ceilings = np.zeros((sample_count, m, m), dtype=float)
    condition_numbers = np.zeros(sample_count, dtype=float)

    for idx in range(sample_count):
        h = validate_spd_matrix(f"hessians[{idx}]", samples[idx], tol)
        ceilings[idx] = h[:m, :m]
        condition_numbers[idx] = float(np.linalg.cond(h))
    mediated = ceiling_mediated_local_quadratic_ensemble(samples, C, ceilings, tolerances=tol)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=mediated.metadata.rank_tol,
        method="coordinate-local-quadratic-ensemble",
        ambient_dim=n,
        support_rank=sample_count,
        visible_dim=m,
        condition_number=float(np.max(condition_numbers)),
        support_restricted=False,
        notes=(
            "samplewise exact local quadratic geometry",
            "coordinate visible split only",
            "clock summaries are descriptive ensemble statistics, not full-law cumulants",
        ),
    )
    return CoordinateLocalQuadraticEnsembleResult(
        phis=mediated.phis,
        ceilings=ceilings,
        lambdas=mediated.lambdas,
        clocks=mediated.clocks,
        hidden_ranks=mediated.hidden_ranks,
        mean_clock=mediated.mean_clock,
        std_clock=mediated.std_clock,
        min_clock=mediated.min_clock,
        max_clock=mediated.max_clock,
        sample_count=sample_count,
        metadata=metadata,
    )


def _validate_hessian_stack(hessians: np.ndarray, tolerances: Tolerances) -> np.ndarray:
    samples = np.asarray(hessians, dtype=float)
    if samples.ndim != 3 or samples.shape[1] != samples.shape[2]:
        raise InputValidationError("hessians must have shape (sample_count, n, n)")
    if samples.shape[0] < 1:
        raise InputValidationError("hessians must contain at least one sample")
    for idx in range(samples.shape[0]):
        validate_spd_matrix(f"hessians[{idx}]", samples[idx], tolerances)
    return samples


def _validate_ceiling_stack(
    ceilings: np.ndarray,
    sample_count: int,
    visible_dim: int,
    tolerances: Tolerances,
) -> np.ndarray:
    array = np.asarray(ceilings, dtype=float)
    if array.ndim == 2:
        if array.shape != (visible_dim, visible_dim):
            raise InputValidationError("2D ceilings must have shape (visible_dim, visible_dim)")
        stack = np.broadcast_to(array, (sample_count, visible_dim, visible_dim)).copy()
    elif array.ndim == 3:
        if array.shape != (sample_count, visible_dim, visible_dim):
            raise InputValidationError("3D ceilings must have shape (sample_count, visible_dim, visible_dim)")
        stack = array.copy()
    else:
        raise InputValidationError("ceilings must be a 2D fixed ceiling or 3D ceiling stack")
    for idx in range(sample_count):
        validate_spd_matrix(f"ceilings[{idx}]", stack[idx], tolerances)
    return stack
