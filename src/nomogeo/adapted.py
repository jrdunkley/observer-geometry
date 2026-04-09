from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .exceptions import InputValidationError
from .types import (
    ClosureAdaptedObserverResult,
    ClosureScoresResult,
    LeakageChannelsResult,
    LinearAlgebraMetadata,
    ObserverComparisonResult,
)
from .validation import (
    Tolerances,
    inv_sqrt_spd,
    rank_cutoff,
    resolve_tolerances,
    sqrt_psd,
    symmetrize,
    to_float_array,
    validate_spd_matrix,
    validate_symmetric_matrix,
)


def whitened_perturbation(
    H: np.ndarray,
    Delta: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return the H-whitened perturbation H^(-1/2) Delta H^(-1/2)."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    delta = validate_symmetric_matrix("Delta", Delta, tol)
    if delta.shape != h.shape:
        raise InputValidationError("Delta must have the same shape as H")
    h_inv_half = inv_sqrt_spd(h, tol)
    return symmetrize(h_inv_half @ delta @ h_inv_half)


def observer_from_subspace(
    H: np.ndarray,
    B: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Build the observer C = B^T H^(1/2) from an H-whitened visible plane."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    basis = _validate_stiefel_basis(B, ambient_dim=h.shape[0], tolerances=tol)
    h_half = sqrt_psd(h, tol)
    return basis.T @ h_half


def closure_scores(
    H: np.ndarray,
    family: Sequence[np.ndarray] | np.ndarray,
    B: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ClosureScoresResult:
    """Compute exact closure leakage L, retained score S, and hidden fraction eta."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    basis = _validate_stiefel_basis(B, ambient_dim=h.shape[0], tolerances=tol)
    whitened_family = _prepare_whitened_family(h, family, tol)
    projector = symmetrize(basis @ basis.T)

    leakage = 0.0
    visible_score = 0.0
    for whitened in whitened_family:
        commutator = whitened @ projector - projector @ whitened
        leakage += 0.5 * float(np.linalg.norm(commutator, ord="fro") ** 2)
        compressed = symmetrize(basis.T @ whitened @ basis)
        visible_score += float(np.linalg.norm(compressed, ord="fro") ** 2)

    total_curvature = leakage + visible_score
    eta_cutoff = max(tol.atol, tol.rtol * max(1.0, total_curvature))
    eta = 0.0 if total_curvature <= eta_cutoff else float(leakage / total_curvature)
    notes = _condition_notes(h, tol)
    if total_curvature <= eta_cutoff:
        notes = notes + ("eta set to zero because total curvature is numerically zero",)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="whitened-closure-scores",
        ambient_dim=h.shape[0],
        support_rank=basis.shape[1],
        visible_dim=basis.shape[1],
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=notes,
    )
    return ClosureScoresResult(
        leakage=leakage,
        visible_score=visible_score,
        eta=eta,
        total_curvature=total_curvature,
        projector=projector,
        metadata=metadata,
    )


def leakage_channels(
    H: np.ndarray,
    Delta: np.ndarray,
    B: np.ndarray,
    tolerances: Tolerances | None = None,
) -> LeakageChannelsResult:
    """Return the exact leakage channels of one perturbation relative to a visible plane."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    basis = _validate_stiefel_basis(B, ambient_dim=h.shape[0], tolerances=tol)
    whitened = whitened_perturbation(h, Delta, tolerances=tol)
    projector = symmetrize(basis @ basis.T)
    visible_operator = symmetrize(basis.T @ whitened @ basis)
    coupling = (np.eye(h.shape[0], dtype=float) - projector) @ whitened @ basis
    left, singular_values, right_t = np.linalg.svd(coupling, full_matrices=False)
    leakage_gram = symmetrize(coupling.T @ coupling)
    cutoff = 10.0 * max(tol.atol, tol.rtol * max(1.0, float(np.linalg.norm(coupling, ord="fro"))))
    keep = singular_values > cutoff
    visible_channels = right_t[keep].T if np.any(keep) else np.zeros((basis.shape[1], 0), dtype=float)
    hidden_channels = left[:, keep] if np.any(keep) else np.zeros((h.shape[0], 0), dtype=float)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="leakage-channel-svd",
        ambient_dim=h.shape[0],
        support_rank=int(np.sum(keep)),
        visible_dim=basis.shape[1],
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=_condition_notes(h, tol),
    )
    return LeakageChannelsResult(
        whitened_perturbation=whitened,
        visible_operator=visible_operator,
        leakage_gram=leakage_gram,
        singular_values=singular_values[keep],
        visible_channel_basis=visible_channels,
        hidden_channel_basis=hidden_channels,
        coupling=coupling,
        projector=projector,
        metadata=metadata,
    )


def compare_observers(
    H: np.ndarray,
    family: Sequence[np.ndarray] | np.ndarray,
    B_left: np.ndarray,
    B_right: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ObserverComparisonResult:
    """Compare two same-rank observers on the same task family."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    left_basis = _validate_stiefel_basis(B_left, ambient_dim=h.shape[0], tolerances=tol)
    right_basis = _validate_stiefel_basis(B_right, ambient_dim=h.shape[0], tolerances=tol)
    if left_basis.shape[1] != right_basis.shape[1]:
        raise InputValidationError("B_left and B_right must have the same visible rank")

    left_scores = closure_scores(h, family, left_basis, tolerances=tol)
    right_scores = closure_scores(h, family, right_basis, tolerances=tol)
    leakage_delta = float(left_scores.leakage - right_scores.leakage)
    visible_score_delta = float(left_scores.visible_score - right_scores.visible_score)
    eta_delta = float(left_scores.eta - right_scores.eta)
    total_curvature_delta = float(left_scores.total_curvature - right_scores.total_curvature)
    left_dominates = (
        left_scores.leakage <= right_scores.leakage + tol.atol
        and left_scores.visible_score >= right_scores.visible_score - tol.atol
        and (
            left_scores.leakage < right_scores.leakage - tol.atol
            or left_scores.visible_score > right_scores.visible_score + tol.atol
        )
    )
    right_dominates = (
        right_scores.leakage <= left_scores.leakage + tol.atol
        and right_scores.visible_score >= left_scores.visible_score - tol.atol
        and (
            right_scores.leakage < left_scores.leakage - tol.atol
            or right_scores.visible_score > left_scores.visible_score + tol.atol
        )
    )
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="closure-score-comparison",
        ambient_dim=h.shape[0],
        support_rank=left_basis.shape[1],
        visible_dim=left_basis.shape[1],
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=_condition_notes(h, tol),
    )
    return ObserverComparisonResult(
        left_scores=left_scores,
        right_scores=right_scores,
        leakage_delta=leakage_delta,
        visible_score_delta=visible_score_delta,
        eta_delta=eta_delta,
        total_curvature_delta=total_curvature_delta,
        left_dominates=left_dominates,
        right_dominates=right_dominates,
        metadata=metadata,
    )


def closure_adapted_observer(
    H: np.ndarray,
    family: Sequence[np.ndarray] | np.ndarray,
    rank: int,
    *,
    mode: str = "commuting_exact",
    tolerances: Tolerances | None = None,
) -> ClosureAdaptedObserverResult:
    """Synthesize an exact closure-adapted observer for the commuting regime."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    if isinstance(rank, bool) or not isinstance(rank, int):
        raise InputValidationError("rank must be an integer")
    if rank < 1 or rank > h.shape[0]:
        raise InputValidationError("rank must satisfy 1 <= rank <= ambient dimension")
    if mode != "commuting_exact":
        raise InputValidationError("only mode='commuting_exact' is currently implemented")

    whitened_family = _prepare_whitened_family(h, family, tol)
    _validate_pairwise_commuting(whitened_family, tol)
    common_basis = _common_eigenbasis(whitened_family, tol)

    spectral_energies = np.zeros(common_basis.shape[1], dtype=float)
    for idx in range(common_basis.shape[1]):
        vector = common_basis[:, idx]
        for whitened in whitened_family:
            spectral_energies[idx] += float(vector.T @ whitened @ whitened @ vector)
    order = np.argsort(spectral_energies)[::-1]
    selected = tuple(int(i) for i in order[:rank])
    basis = common_basis[:, list(selected)]
    observer = observer_from_subspace(h, basis, tolerances=tol)
    scores = closure_scores(h, family, basis, tolerances=tol)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="commuting-common-eigenbasis",
        ambient_dim=h.shape[0],
        support_rank=rank,
        visible_dim=rank,
        condition_number=_condition_number_spd(h),
        support_restricted=False,
        notes=("commuting exact solver",) + _condition_notes(h, tol),
    )
    return ClosureAdaptedObserverResult(
        B=basis,
        C=observer,
        projector=scores.projector,
        scores=scores,
        common_basis=common_basis,
        spectral_energies=spectral_energies,
        selected_indices=selected,
        metadata=metadata,
    )


def _validate_stiefel_basis(B: np.ndarray, ambient_dim: int, tolerances: Tolerances) -> np.ndarray:
    basis = to_float_array("B", B)
    if basis.shape[0] != ambient_dim:
        raise InputValidationError("B has incompatible ambient dimension")
    if basis.shape[1] < 1 or basis.shape[1] > ambient_dim:
        raise InputValidationError("B must have between 1 and ambient_dim columns")
    gram = symmetrize(basis.T @ basis)
    if not np.allclose(gram, np.eye(basis.shape[1]), atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError("B must have orthonormal columns")
    return basis


def _prepare_whitened_family(
    H: np.ndarray,
    family: Sequence[np.ndarray] | np.ndarray,
    tolerances: Tolerances,
) -> list[np.ndarray]:
    if isinstance(family, np.ndarray):
        members = [family]
    else:
        members = list(family)
    if not members:
        raise InputValidationError("family must contain at least one perturbation")
    whitened = []
    for idx, delta in enumerate(members):
        validated = validate_symmetric_matrix(f"family[{idx}]", delta, tolerances)
        if validated.shape != H.shape:
            raise InputValidationError(f"family[{idx}] must have the same shape as H")
        whitened.append(whitened_perturbation(H, validated, tolerances=tolerances))
    return whitened


def _validate_pairwise_commuting(family: Sequence[np.ndarray], tolerances: Tolerances) -> None:
    for i, left in enumerate(family):
        for j in range(i + 1, len(family)):
            right = family[j]
            commutator = left @ right - right @ left
            scale = max(
                1.0,
                float(np.linalg.norm(left, ord="fro")),
                float(np.linalg.norm(right, ord="fro")),
            )
            cutoff = 10.0 * max(tolerances.atol, tolerances.rtol * scale)
            if float(np.linalg.norm(commutator, ord="fro")) > cutoff:
                raise InputValidationError("mode='commuting_exact' requires a pairwise commuting whitened family")


def _condition_number_spd(H: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(H)
    return float(np.max(eigenvalues) / np.min(eigenvalues))


def _condition_notes(H: np.ndarray, tolerances: Tolerances) -> tuple[str, ...]:
    condition = _condition_number_spd(H)
    if condition >= 1e9:
        return (f"H is extremely ill-conditioned (cond={condition:.3e}); whitening may be numerically unreliable",)
    if condition >= 1e8:
        return (f"H is severely ill-conditioned (cond={condition:.3e}); interpret closure numerics with care",)
    if condition >= 1e6:
        return (f"H is highly conditioned (cond={condition:.3e}); whitening diagnostics should be checked",)
    return ()


def _group_consecutive(values: np.ndarray, cutoff: float) -> list[tuple[int, int]]:
    if values.size == 0:
        return []
    groups: list[tuple[int, int]] = []
    start = 0
    for idx in range(1, values.size):
        if abs(float(values[idx] - values[idx - 1])) > cutoff:
            groups.append((start, idx))
            start = idx
    groups.append((start, values.size))
    return groups


def _common_eigenbasis(family: Sequence[np.ndarray], tolerances: Tolerances) -> np.ndarray:
    n = family[0].shape[0]
    blocks = [np.eye(n, dtype=float)]
    for operator in family:
        refined_blocks: list[np.ndarray] = []
        for block in blocks:
            reduced = symmetrize(block.T @ operator @ block)
            eigenvalues, eigenvectors = np.linalg.eigh(reduced)
            cutoff = 10.0 * rank_cutoff(eigenvalues, tolerances)
            for start, stop in _group_consecutive(eigenvalues, cutoff):
                refined_blocks.append(block @ eigenvectors[:, start:stop])
        blocks = refined_blocks

    basis = np.column_stack(blocks)
    for operator in family:
        reduced = basis.T @ operator @ basis
        diagonal = np.diag(np.diag(reduced))
        residual = float(np.linalg.norm(reduced - diagonal, ord="fro"))
        scale = max(1.0, float(np.linalg.norm(operator, ord="fro")))
        cutoff = 50.0 * max(tolerances.atol, tolerances.rtol * scale)
        if residual > cutoff:
            raise InputValidationError("failed to construct a common eigenbasis for the commuting family")
    return basis
