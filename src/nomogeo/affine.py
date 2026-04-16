from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .exceptions import InputValidationError
from .regime_types import TowerEliminationRecord, TowerStageRecord
from .types import (
    AffineHiddenBranchReversalResult,
    AffineHiddenReductionResult,
    AffineHiddenStageResult,
    GuardedFibreDominanceResult,
    LinearAlgebraMetadata,
)
from .validation import (
    Tolerances,
    logdet_spd,
    rank_cutoff,
    resolve_tolerances,
    solve_spd,
    symmetrize,
    validate_spd_matrix,
)


def variable_precision_affine_hidden_reduction(
    A: np.ndarray | float,
    J: np.ndarray,
    D: np.ndarray,
    tolerances: Tolerances | None = None,
) -> AffineHiddenReductionResult:
    """Reduce an exact affine-hidden Gaussian fibre with visible-dependent precision.

    The exact law sector is

    ``p(v,h) proportional to exp(-A(v) - 1/2 h^T D(v) h - J(v)^T h)``.

    The returned visible action omits the global additive constant independent
    of ``v``:

    ``A(v) + 1/2 log det D(v) - 1/2 J(v)^T D(v)^(-1) J(v)``.
    """
    tol = resolve_tolerances(tolerances)
    action, coupling, precision, sample_shape, hidden_dim = _prepare_affine_inputs(A, J, D, tol)
    flat_count = int(np.prod(sample_shape, dtype=int)) if sample_shape else 1
    flat_action = action.reshape(flat_count)
    flat_coupling = coupling.reshape(flat_count, hidden_dim)
    flat_precision = precision.reshape(flat_count, hidden_dim, hidden_dim)

    hidden_mean = np.zeros((flat_count, hidden_dim), dtype=float)
    variational_action = np.zeros(flat_count, dtype=float)
    fibre_volume = np.zeros(flat_count, dtype=float)
    condition_numbers = np.zeros(flat_count, dtype=float)

    for idx in range(flat_count):
        d_i = validate_spd_matrix(f"D[{idx}]", flat_precision[idx], tol)
        j_i = flat_coupling[idx]
        solved = solve_spd(d_i, j_i)
        hidden_mean[idx] = -solved
        variational_action[idx] = float(flat_action[idx] - 0.5 * j_i.T @ solved)
        fibre_volume[idx] = 0.5 * logdet_spd(d_i)
        condition_numbers[idx] = float(np.linalg.cond(d_i))

    visible_action = variational_action + fibre_volume
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="variable-precision-affine-hidden-reduction",
        ambient_dim=hidden_dim,
        support_rank=flat_count,
        condition_number=float(np.max(condition_numbers)),
        support_restricted=False,
        notes=(
            "exact affine-hidden Gaussian fibre reduction",
            "visible action is defined up to an additive constant independent of visible coordinates",
            "not a generic non-Gaussian marginalisation engine",
        ),
    )
    return AffineHiddenReductionResult(
        action=action.copy(),
        coupling=coupling.copy(),
        hidden_precision=precision.copy(),
        hidden_mean=hidden_mean.reshape(sample_shape + (hidden_dim,)),
        variational_action=variational_action.reshape(sample_shape),
        fibre_volume=fibre_volume.reshape(sample_shape),
        visible_action=visible_action.reshape(sample_shape),
        sample_shape=sample_shape,
        metadata=metadata,
    )


def staged_affine_hidden_elimination(
    A: float,
    J: np.ndarray,
    D: np.ndarray,
    eliminate: Sequence[int],
    tolerances: Tolerances | None = None,
) -> AffineHiddenStageResult:
    """Eliminate a declared hidden block from one affine-hidden fibre."""
    tol = resolve_tolerances(tolerances)
    if not np.isscalar(A):
        raise InputValidationError("A must be a scalar for staged elimination")
    coupling = np.asarray(J, dtype=float)
    if coupling.ndim != 1:
        raise InputValidationError("J must be a 1D coupling vector for staged elimination")
    precision = validate_spd_matrix("D", D, tol)
    hidden_dim = precision.shape[0]
    if coupling.shape != (hidden_dim,):
        raise InputValidationError("J must have shape (hidden_dim,)")

    eliminated = _validate_eliminate_indices(eliminate, hidden_dim)
    keep = tuple(idx for idx in range(hidden_dim) if idx not in eliminated)
    d_ee = precision[np.ix_(eliminated, eliminated)]
    j_e = coupling[list(eliminated)]
    solved_j = solve_spd(d_ee, j_e)
    action_shift = 0.5 * logdet_spd(d_ee) - 0.5 * float(j_e.T @ solved_j)
    shifted_action = float(A) + action_shift
    condition_number = float(np.linalg.cond(d_ee))

    if keep:
        d_kk = precision[np.ix_(keep, keep)]
        d_ke = precision[np.ix_(keep, eliminated)]
        d_ek = precision[np.ix_(eliminated, keep)]
        j_k = coupling[list(keep)]
        schur = symmetrize(d_kk - d_ke @ solve_spd(d_ee, d_ek))
        shifted_coupling = j_k - d_ke @ solved_j
        visible_action = None
        output_precision = schur
        output_coupling = shifted_coupling
        support_rank = len(keep)
    else:
        visible_action = shifted_action
        output_precision = np.zeros((0, 0), dtype=float)
        output_coupling = np.zeros(0, dtype=float)
        support_rank = 0

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=rank_cutoff(np.linalg.eigvalsh(d_ee), tol),
        method="staged-affine-hidden-elimination",
        ambient_dim=hidden_dim,
        support_rank=support_rank,
        condition_number=condition_number,
        support_restricted=False,
        notes=(
            "one-stage exact Schur elimination in the affine-hidden sector",
            "visible_action is populated only when all hidden coordinates have been eliminated",
        ),
    )
    return AffineHiddenStageResult(
        action=shifted_action,
        coupling=output_coupling,
        hidden_precision=output_precision,
        eliminated_indices=eliminated,
        kept_indices=keep,
        action_shift=action_shift,
        visible_action=visible_action,
        metadata=metadata,
    )


def tower_affine_hidden_elimination(
    A: float,
    J: np.ndarray,
    D: np.ndarray,
    stages: Sequence[tuple[str, Sequence[int]]],
    tolerances: Tolerances | None = None,
) -> TowerEliminationRecord:
    """Sequential Schur-complement elimination of multiple hidden blocks.

    Chains exact affine-hidden elimination across a declared sequence of
    hidden blocks.  At each stage, the specified block is eliminated via
    Schur complement, accumulating the log-determinant contribution.

    This implements the sequential tower law from Quotient_Descent_MASTER:
    Schur complement composition is exact and order-invariant.

    Parameters
    ----------
    A : float
        Initial scalar action.
    J : ndarray, shape (hidden_dim,)
        Initial coupling vector.
    D : ndarray, shape (hidden_dim, hidden_dim)
        Initial precision matrix (SPD).
    stages : sequence of (block_name, eliminate_indices) pairs
        Each stage specifies a human-readable block name and a sequence
        of indices (relative to the *current* remaining hidden space)
        to eliminate at that stage.

        IMPORTANT: indices at each stage are relative to the remaining
        hidden coordinates after all previous stages, not to the
        original D.  This is the natural convention for layered models
        (eliminate layer 1, then layer 2 of the remainder, etc.).
    tolerances : Tolerances, optional

    Returns
    -------
    TowerEliminationRecord
        Per-stage records with block names, dimensions, log-det
        increments, and condition numbers.  Plus the final reduced
        action, coupling, and precision.

    Notes
    -----
    The tower does NOT construct a ReducedLocalDatum.  That is the
    caller's responsibility, because the caller knows the cone/orbit
    context of the visible problem.  The tower only handles the exact
    hidden-elimination layer.

    Contact surfaces: Kalman smoothing, state-space evidence, banded
    Gaussian latent models, hierarchical random effects.
    """
    tol = resolve_tolerances(tolerances)
    if not np.isscalar(A):
        raise InputValidationError("A must be a scalar for tower elimination")
    coupling = np.asarray(J, dtype=float).copy()
    if coupling.ndim != 1:
        raise InputValidationError("J must be a 1D coupling vector")
    precision = validate_spd_matrix("D", D, tol).copy()
    hidden_dim = precision.shape[0]
    if coupling.shape != (hidden_dim,):
        raise InputValidationError("J must have shape (hidden_dim,)")

    if not stages:
        raise InputValidationError("stages must contain at least one elimination stage")

    current_action = float(A)
    current_coupling = coupling
    current_precision = precision
    total_half_log_det = 0.0
    stage_records: list[TowerStageRecord] = []

    for stage_idx, (block_name, eliminate) in enumerate(stages):
        current_dim = current_precision.shape[0]
        if current_dim == 0:
            raise InputValidationError(
                f"Stage {stage_idx} ({block_name!r}): no hidden "
                f"coordinates remain to eliminate"
            )

        # Validate and resolve indices for this stage.
        eliminated = _validate_eliminate_indices(eliminate, current_dim)
        keep = tuple(idx for idx in range(current_dim) if idx not in eliminated)

        d_ee = current_precision[np.ix_(eliminated, eliminated)]
        j_e = current_coupling[list(eliminated)]
        solved_j = solve_spd(d_ee, j_e)

        half_log_det_inc = 0.5 * logdet_spd(d_ee)
        action_shift = half_log_det_inc - 0.5 * float(j_e.T @ solved_j)
        condition_number = float(np.linalg.cond(d_ee))

        current_action = current_action + action_shift
        total_half_log_det += half_log_det_inc

        stage_records.append(TowerStageRecord(
            block_name=str(block_name),
            hidden_dim=len(eliminated),
            half_log_det_increment=half_log_det_inc,
            schur_condition_number=condition_number,
        ))

        if keep:
            d_kk = current_precision[np.ix_(keep, keep)]
            d_ke = current_precision[np.ix_(keep, eliminated)]
            d_ek = current_precision[np.ix_(eliminated, keep)]
            j_k = current_coupling[list(keep)]
            current_precision = symmetrize(
                d_kk - d_ke @ solve_spd(d_ee, d_ek)
            )
            current_coupling = j_k - d_ke @ solved_j
        else:
            current_precision = np.zeros((0, 0), dtype=float)
            current_coupling = np.zeros(0, dtype=float)

    return TowerEliminationRecord(
        stages=tuple(stage_records),
        total_half_log_det=total_half_log_det,
        final_action=current_action,
        final_coupling=current_coupling,
        final_precision=current_precision,
    )


def affine_hidden_branch_reversal(
    variational_action: np.ndarray,
    fibre_volume: np.ndarray,
    tolerances: Tolerances | None = None,
) -> AffineHiddenBranchReversalResult:
    """Detect branch changes caused by affine-hidden fibre-volume terms.

    Lower action is better. This helper is exact only for finite branch data
    already reduced in the affine-hidden Gaussian-fibre sector.
    """
    tol = resolve_tolerances(tolerances)
    var = np.asarray(variational_action, dtype=float)
    fib = np.asarray(fibre_volume, dtype=float)
    if var.ndim != 1 or fib.ndim != 1 or var.shape != fib.shape:
        raise InputValidationError("variational_action and fibre_volume must be 1D arrays with matching shape")
    if var.size < 1:
        raise InputValidationError("at least one branch is required")
    if not np.all(np.isfinite(var)) or not np.all(np.isfinite(fib)):
        raise InputValidationError("branch actions must be finite")
    visible = var + fib
    cutoff = max(tol.atol, tol.rtol * max(1.0, float(np.max(np.abs(visible)))))
    var_min = float(np.min(var))
    vis_min = float(np.min(visible))
    var_winners = tuple(int(idx) for idx in np.flatnonzero(var <= var_min + cutoff))
    vis_winners = tuple(int(idx) for idx in np.flatnonzero(visible <= vis_min + cutoff))
    preserved = bool(set(var_winners) & set(vis_winners))

    count = var.size
    reversal_matrix = np.zeros((count, count), dtype=float)
    for a in range(count):
        for b in range(count):
            if a == b:
                continue
            var_gap = var[b] - var[a]
            fib_advantage = fib[a] - fib[b]
            if var_gap > cutoff and fib_advantage > var_gap + cutoff:
                reversal_matrix[a, b] = 1.0
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="affine-hidden-branch-reversal",
        ambient_dim=count,
        support_rank=len(vis_winners),
        visible_dim=None,
        notes=(
            "finite branch comparison in the exact affine-hidden Gaussian-fibre sector",
            "lower visible action wins",
            "not a generic marginalisation or non-Gaussian branch probability",
        ),
    )
    return AffineHiddenBranchReversalResult(
        variational_action=var.copy(),
        fibre_volume=fib.copy(),
        visible_action=visible,
        variational_winners=var_winners,
        visible_winners=vis_winners,
        preserved=preserved,
        reversal=not preserved,
        branch_reversal_matrix=reversal_matrix,
        metadata=metadata,
    )


def guarded_fibre_dominance(
    variational_action: np.ndarray,
    fibre_volume: np.ndarray,
    sample_weights: np.ndarray | None = None,
    denominator_floor: float = 1e-12,
    norm: str = "l2",
    tolerances: Tolerances | None = None,
) -> GuardedFibreDominanceResult:
    """Return guarded centered norms for fibre-volume dominance diagnostics."""
    tol = resolve_tolerances(tolerances)
    var = np.asarray(variational_action, dtype=float)
    fib = np.asarray(fibre_volume, dtype=float)
    if var.shape != fib.shape or var.size < 1:
        raise InputValidationError("variational_action and fibre_volume must have matching nonempty shape")
    if not np.all(np.isfinite(var)) or not np.all(np.isfinite(fib)):
        raise InputValidationError("diagnostic arrays must be finite")
    floor = float(denominator_floor)
    if not np.isfinite(floor) or floor < 0.0:
        raise InputValidationError("denominator_floor must be nonnegative")

    weights = None
    weight_sum = None
    if sample_weights is not None:
        weights = np.asarray(sample_weights, dtype=float)
        if weights.shape != var.shape:
            raise InputValidationError("sample_weights must match action shape")
        if not np.all(np.isfinite(weights)) or np.any(weights < 0.0) or float(np.sum(weights)) <= 0.0:
            raise InputValidationError("sample_weights must be finite, nonnegative, and nonzero")
        weight_sum = float(np.sum(weights))
        weights = weights / weight_sum

    fibre_norm = _centered_norm(fib, weights, norm)
    variational_norm = _centered_norm(var, weights, norm)
    ratio_defined = bool(variational_norm >= floor)
    ratio = None if not ratio_defined else float(fibre_norm / variational_norm)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="guarded-fibre-dominance",
        ambient_dim=int(var.size),
        support_rank=int(var.size),
        visible_dim=None,
        notes=(
            "diagnostic centered-norm tuple with declared norm and denominator floor",
            "ratio is undefined when the variational centered norm is below the floor",
            "not an invariant naked fibre-dominance scalar",
        ),
    )
    return GuardedFibreDominanceResult(
        fibre_centered_norm=float(fibre_norm),
        variational_centered_norm=float(variational_norm),
        ratio=ratio,
        ratio_defined=ratio_defined,
        denominator_floor=floor,
        norm=norm,
        sample_weight_sum=weight_sum,
        metadata=metadata,
    )


def _prepare_affine_inputs(
    A: np.ndarray | float,
    J: np.ndarray,
    D: np.ndarray,
    tolerances: Tolerances,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, ...], int]:
    precision = np.asarray(D, dtype=float)
    if precision.ndim < 2 or precision.shape[-1] != precision.shape[-2]:
        raise InputValidationError("D must have shape sample_shape + (hidden_dim, hidden_dim)")
    hidden_dim = precision.shape[-1]
    sample_shape = tuple(int(dim) for dim in precision.shape[:-2])
    if hidden_dim < 1:
        raise InputValidationError("hidden_dim must be positive")

    coupling = np.asarray(J, dtype=float)
    expected_coupling_shape = sample_shape + (hidden_dim,)
    if coupling.shape != expected_coupling_shape:
        raise InputValidationError("J must have shape sample_shape + (hidden_dim,)")

    action_array = np.asarray(A, dtype=float)
    if action_array.shape == ():
        action = np.broadcast_to(action_array, sample_shape).astype(float).copy()
    elif action_array.shape == sample_shape:
        action = action_array.astype(float, copy=True)
    else:
        raise InputValidationError("A must be scalar or have shape sample_shape")

    flat_precision = precision.reshape((-1, hidden_dim, hidden_dim)) if sample_shape else precision.reshape((1, hidden_dim, hidden_dim))
    for idx in range(flat_precision.shape[0]):
        validate_spd_matrix(f"D[{idx}]", flat_precision[idx], tolerances)

    return action, coupling.copy(), precision.copy(), sample_shape, hidden_dim


def _centered_norm(values: np.ndarray, weights: np.ndarray | None, norm: str) -> float:
    if weights is None:
        centered = values - float(np.mean(values))
        if norm == "l2":
            return float(np.linalg.norm(centered.ravel()) / np.sqrt(centered.size))
        if norm == "fro":
            return float(np.linalg.norm(centered.ravel()))
        if norm == "range":
            return float(np.max(centered) - np.min(centered))
    else:
        mean = float(np.sum(weights * values))
        centered = values - mean
        if norm == "l2":
            return float(np.sqrt(np.sum(weights * centered * centered)))
        if norm == "fro":
            return float(np.linalg.norm(centered.ravel()))
        if norm == "range":
            return float(np.max(centered) - np.min(centered))
    raise InputValidationError("norm must be 'l2', 'fro', or 'range'")


def _validate_eliminate_indices(eliminate: Sequence[int], hidden_dim: int) -> tuple[int, ...]:
    try:
        indices = tuple(int(idx) for idx in eliminate)
    except TypeError as exc:
        raise InputValidationError("eliminate must be a nonempty sequence of indices") from exc
    if not indices:
        raise InputValidationError("eliminate must contain at least one index")
    if len(set(indices)) != len(indices):
        raise InputValidationError("eliminate indices must be unique")
    if min(indices) < 0 or max(indices) >= hidden_dim:
        raise InputValidationError("eliminate indices are out of range")
    return indices
