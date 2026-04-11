from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import scipy.linalg as la

from .exceptions import InputValidationError
from .types import (
    DeclaredFrontierLocalCertificateResult,
    DeclaredLadderDimensionCostResult,
    ExactBranchHessianResult,
    GeneralGraphFrontierHessianResult,
    LinearAlgebraMetadata,
    WeightedFamilyFrontierResult,
)
from .validation import (
    Tolerances,
    rank_cutoff,
    resolve_tolerances,
    symmetrize,
    to_float_array,
    validate_symmetric_matrix,
)


def weighted_family_frontier_scores(
    family: Sequence[np.ndarray] | np.ndarray,
    B: np.ndarray,
    weights: np.ndarray | None = None,
    mu: float = 0.0,
    tolerances: Tolerances | None = None,
) -> WeightedFamilyFrontierResult:
    """Evaluate the finite weighted-family leakage/visibility frontier."""
    tol = resolve_tolerances(tolerances)
    operators = _validate_operator_family(family, tol)
    n = operators[0].shape[0]
    basis = _validate_stiefel_basis(B, n, tol, allow_full_rank=True)
    weights_array = _validate_weights(weights, len(operators))
    penalty = _validate_mu(mu)
    projector = symmetrize(basis @ basis.T)

    moment = np.zeros((n, n), dtype=float)
    leakage = 0.0
    visible_score = 0.0
    for weight, operator in zip(weights_array, operators):
        moment += weight * operator @ operator
        compressed = basis.T @ operator @ basis
        visible_score += float(weight * np.linalg.norm(compressed, ord="fro") ** 2)
        commutator = operator @ projector - projector @ operator
        leakage += float(weight * 0.5 * np.linalg.norm(commutator, ord="fro") ** 2)
    moment = symmetrize(moment)
    captured = float(np.trace(projector @ moment))
    residual = float(captured - visible_score - leakage)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="weighted-family-frontier-scores",
        ambient_dim=n,
        support_rank=basis.shape[1],
        visible_dim=basis.shape[1],
        support_restricted=False,
        notes=(
            "finite weighted-family quadratic frontier evaluator",
            "not a noncommuting optimiser or full-law selector",
        ),
    )
    return WeightedFamilyFrontierResult(
        leakage=leakage,
        visible_score=visible_score,
        captured_curvature=captured,
        energy_split_residual=residual,
        penalized_score=float(visible_score - penalty * leakage),
        moment_operator=moment,
        projector=projector,
        weights=weights_array,
        metadata=metadata,
    )


def exact_branch_hessian(
    family: Sequence[np.ndarray] | np.ndarray,
    B: np.ndarray,
    weights: np.ndarray | None = None,
    mu: float = 0.0,
    tolerances: Tolerances | None = None,
) -> ExactBranchHessianResult:
    """Return the exact-branch Hessian contract on ``Hom(U, U_perp)``.

    Positive eigenvalues of the returned contract operator correspond to a
    local maximum of the penalised frontier, because the second variation is
    ``-2 * hessian_contract``.
    """
    tol = resolve_tolerances(tolerances)
    operators = _validate_operator_family(family, tol)
    n = operators[0].shape[0]
    basis = _validate_stiefel_basis(B, n, tol, allow_full_rank=False)
    complement = np.asarray(la.null_space(basis.T, rcond=max(tol.atol, tol.rtol)), dtype=float)
    if complement.shape[1] == 0:
        raise InputValidationError("exact_branch_hessian requires a nontrivial complement")
    weights_array = _validate_weights(weights, len(operators))
    penalty = _validate_mu(mu)
    m = basis.shape[1]
    hdim = complement.shape[1]

    blocks_u: list[np.ndarray] = []
    blocks_w: list[np.ndarray] = []
    off_block_norm = 0.0
    moment_u = np.zeros((m, m), dtype=float)
    moment_w = np.zeros((hdim, hdim), dtype=float)
    for weight, operator in zip(weights_array, operators):
        a_u = symmetrize(basis.T @ operator @ basis)
        a_w = symmetrize(complement.T @ operator @ complement)
        off = complement.T @ operator @ basis
        off_block_norm = max(off_block_norm, float(np.linalg.norm(off, ord="fro")))
        blocks_u.append(a_u)
        blocks_w.append(a_w)
        moment_u += weight * a_u @ a_u
        moment_w += weight * a_w @ a_w

    family_scale = max(1.0, max(float(np.linalg.norm(op, ord="fro")) for op in operators))
    off_cutoff = 50.0 * max(tol.atol, tol.rtol * family_scale)
    if off_block_norm > off_cutoff:
        raise InputValidationError("B is not an exact branch for the supplied family")

    tangent_basis = _tangent_basis(hdim, m)
    dim = len(tangent_basis)
    contract = np.zeros((dim, dim), dtype=float)
    for col, x in enumerate(tangent_basis):
        gx = moment_w @ x - x @ moment_u
        for row, y in enumerate(tangent_basis):
            intertwiner = 0.0
            for weight, a_u, a_w in zip(weights_array, blocks_u, blocks_w):
                cx = a_w @ x - x @ a_u
                cy = a_w @ y - y @ a_u
                intertwiner += float(weight * np.sum(cy * cx))
            gap_term = float(np.sum(y * gx))
            contract[row, col] = (1.0 + penalty) * intertwiner - gap_term
    contract = symmetrize(contract)
    second_variation = -2.0 * contract
    eigenvalues = np.linalg.eigvalsh(contract)
    cutoff = rank_cutoff(eigenvalues, tol)
    nullity = int(np.sum(np.abs(eigenvalues) <= cutoff))
    min_eig = float(np.min(eigenvalues))
    max_eig = float(np.max(eigenvalues))
    if min_eig > cutoff:
        status = "strict_max"
    elif max_eig < -cutoff:
        status = "local_min"
    elif min_eig < -cutoff and max_eig > cutoff:
        status = "saddle"
    else:
        status = "degenerate"

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="exact-branch-hessian",
        ambient_dim=n,
        support_rank=basis.shape[1],
        visible_dim=basis.shape[1],
        support_restricted=False,
        notes=(
            "fixed-rank exact-branch Hessian diagnostic",
            "positive contract eigenvalues indicate a local maximum of the penalised frontier",
            "not a branch optimiser or probability-support classifier",
        ),
    )
    return ExactBranchHessianResult(
        hessian_contract=contract,
        second_variation_operator=second_variation,
        eigenvalues=eigenvalues,
        min_eigenvalue=min_eig,
        nullity=nullity,
        status=status,
        off_block_norm=off_block_norm,
        basis=basis,
        complement_basis=complement,
        weights=weights_array,
        metadata=metadata,
    )


def declared_ladder_dimension_cost_intervals(
    scores: np.ndarray,
    dimensions: np.ndarray,
    tolerances: Tolerances | None = None,
) -> DeclaredLadderDimensionCostResult:
    """Return exact cost intervals for a declared finite observer ladder.

    Candidate ``j`` is scored by ``scores[j] - c * dimensions[j]`` for
    ``c >= 0``. The result is a finite declared-ladder phase diagram, not a
    global observer optimizer.
    """
    tol = resolve_tolerances(tolerances)
    score_array = np.asarray(scores, dtype=float)
    dim_array = np.asarray(dimensions, dtype=float)
    if score_array.ndim != 1 or dim_array.ndim != 1 or score_array.shape != dim_array.shape:
        raise InputValidationError("scores and dimensions must be 1D arrays with matching shape")
    if score_array.size < 1:
        raise InputValidationError("at least one candidate is required")
    if not np.all(np.isfinite(score_array)) or not np.all(np.isfinite(dim_array)):
        raise InputValidationError("scores and dimensions must be finite")
    if np.any(dim_array < 0.0):
        raise InputValidationError("dimensions must be nonnegative")

    count = score_array.size
    crossings = np.full((count, count), np.nan, dtype=float)
    lower = np.zeros(count, dtype=float)
    upper = np.full(count, np.inf, dtype=float)
    nonempty = np.ones(count, dtype=bool)
    cutoff = max(tol.atol, tol.rtol * max(1.0, float(np.max(np.abs(score_array)))))

    for a in range(count):
        for b in range(count):
            if a == b:
                continue
            coeff = float(dim_array[b] - dim_array[a])
            rhs = float(score_array[b] - score_array[a])
            if abs(coeff) <= cutoff:
                if rhs > cutoff:
                    nonempty[a] = False
                continue
            crossing = rhs / coeff
            crossings[a, b] = crossing
            if coeff > 0.0:
                lower[a] = max(lower[a], crossing)
            else:
                upper[a] = min(upper[a], crossing)
        if upper[a] < 0.0:
            nonempty[a] = False
        lower[a] = max(lower[a], 0.0)
        if lower[a] > upper[a] + cutoff:
            nonempty[a] = False

    adjusted_zero = score_array.copy()
    max_zero = float(np.max(adjusted_zero))
    zero_winners = tuple(int(idx) for idx in np.flatnonzero(adjusted_zero >= max_zero - cutoff))
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="declared-ladder-dimension-cost-intervals",
        ambient_dim=count,
        support_rank=int(np.sum(nonempty)),
        visible_dim=None,
        notes=(
            "finite declared-ladder dimension-cost comparison",
            "not a global Grassmannian optimizer or observer discovery algorithm",
        ),
    )
    return DeclaredLadderDimensionCostResult(
        scores=score_array.copy(),
        dimensions=dim_array.copy(),
        pairwise_crossings=crossings,
        interval_lower=lower,
        interval_upper=upper,
        interval_nonempty=nonempty,
        winner_at_zero=zero_winners,
        metadata=metadata,
    )


def general_graph_frontier_hessian(
    family: Sequence[np.ndarray] | np.ndarray,
    B: np.ndarray,
    weights: np.ndarray | None = None,
    mu: float = 0.0,
    complement_basis: np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> GeneralGraphFrontierHessianResult:
    """Return the declared-observer graph-chart frontier gradient and Hessian.

    This is exact local quadratic geometry at the supplied observer. It
    contains the exact-branch Hessian as the off-block-zero sector, but it is
    not a branch optimizer and not a full-law non-Gaussian selector.
    """
    tol = resolve_tolerances(tolerances)
    operators = _validate_operator_family(family, tol)
    n = operators[0].shape[0]
    basis = _validate_stiefel_basis(B, n, tol, allow_full_rank=False)
    complement = _validate_or_make_complement(complement_basis, basis, tol)
    weights_array = _validate_weights(weights, len(operators))
    penalty = _validate_mu(mu)
    m = basis.shape[1]
    hdim = complement.shape[1]

    gradient = np.zeros((hdim, m), dtype=float)
    off_block_norm = 0.0
    blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for weight, operator in zip(weights_array, operators):
        u = symmetrize(basis.T @ operator @ basis)
        d = symmetrize(complement.T @ operator @ complement)
        e = complement.T @ operator @ basis
        k2 = operator @ operator
        k2_uu = symmetrize(basis.T @ k2 @ basis)
        k2_ww = symmetrize(complement.T @ k2 @ complement)
        off_block_norm = max(off_block_norm, float(np.linalg.norm(e, ord="fro")))
        gradient += 2.0 * weight * ((2.0 + penalty) * e @ u - penalty * d @ e)
        blocks.append((u, d, e, k2_uu, k2_ww))

    tangent_basis = _tangent_basis(hdim, m)
    dim = len(tangent_basis)
    hessian = np.zeros((dim, dim), dtype=float)
    for col, x in enumerate(tangent_basis):
        for row, y in enumerate(tangent_basis):
            hessian[row, col] = 0.25 * (
                _graph_second_variation_quadratic(x + y, blocks, weights_array, penalty)
                - _graph_second_variation_quadratic(x - y, blocks, weights_array, penalty)
            )
    hessian = symmetrize(hessian)
    eigenvalues = np.linalg.eigvalsh(hessian)
    cutoff = rank_cutoff(eigenvalues, tol)
    stationarity_residual = float(np.linalg.norm(gradient, ord="fro"))
    grad_cutoff = max(tol.atol, tol.rtol * max(1.0, stationarity_residual))
    nullity = int(np.sum(np.abs(eigenvalues) <= cutoff))
    min_eig = float(np.min(eigenvalues))
    max_eig = float(np.max(eigenvalues))
    if stationarity_residual > grad_cutoff:
        status = "nonstationary"
    elif max_eig < -cutoff:
        status = "strict_local_max"
    elif min_eig > cutoff:
        status = "strict_local_min"
    elif min_eig < -cutoff and max_eig > cutoff:
        status = "saddle"
    else:
        status = "stationary_degenerate"

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="general-graph-frontier-hessian",
        ambient_dim=n,
        support_rank=basis.shape[1],
        visible_dim=basis.shape[1],
        support_restricted=False,
        notes=(
            "declared-observer exact local quadratic graph-chart Hessian",
            "raw matrix entries depend on the declared complement basis",
            "not a branch optimizer, probability-support classifier, or full-law selector",
        ),
    )
    return GeneralGraphFrontierHessianResult(
        gradient=gradient,
        gradient_vector=gradient.ravel().copy(),
        hessian_operator=hessian,
        second_variation_operator=hessian.copy(),
        eigenvalues=eigenvalues,
        min_eigenvalue=min_eig,
        max_eigenvalue=max_eig,
        nullity=nullity,
        status=status,
        stationarity_residual=stationarity_residual,
        off_block_norm=off_block_norm,
        basis=basis,
        complement_basis=complement,
        weights=weights_array,
        metadata=metadata,
    )


def declared_frontier_local_certificate(
    family: Sequence[np.ndarray] | np.ndarray,
    B: np.ndarray,
    weights: np.ndarray | None = None,
    mu: float = 0.0,
    mode: str = "max",
    rho: float = 1.0,
    complement_basis: np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> DeclaredFrontierLocalCertificateResult:
    """Return a sufficient local certificate for a declared frontier observer."""
    tol = resolve_tolerances(tolerances)
    if mode not in {"max", "min"}:
        raise InputValidationError("mode must be 'max' or 'min'")
    radius = float(rho)
    if not np.isfinite(radius) or radius <= 0.0 or radius > 1.0:
        raise InputValidationError("rho must satisfy 0 < rho <= 1")
    graph = general_graph_frontier_hessian(
        family,
        B,
        weights=weights,
        mu=mu,
        complement_basis=complement_basis,
        tolerances=tol,
    )
    eps = graph.stationarity_residual
    if mode == "max":
        lambda_margin = max(0.0, -graph.max_eigenvalue)
    else:
        lambda_margin = max(0.0, graph.min_eigenvalue)

    operators = _validate_operator_family(family, tol)
    weights_array = _validate_weights(weights, len(operators))
    penalty = _validate_mu(mu)
    lipschitz = _frontier_lipschitz_bound(operators, weights_array, penalty, graph.basis.shape[1], operators[0].shape[0], radius)
    if lambda_margin > 0.0:
        r0 = min(radius, np.inf if lipschitz == 0.0 else lambda_margin / (2.0 * lipschitz))
        left = 4.0 * eps / lambda_margin
        displacement = 2.0 * eps / lambda_margin
        passes = bool(left < r0)
    else:
        r0 = 0.0
        left = None
        displacement = None
        passes = False
    kind = f"sufficient_local_{mode}" if passes else "vacuous"
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=graph.metadata.rank_tol,
        method="declared-frontier-local-certificate",
        ambient_dim=graph.metadata.ambient_dim,
        support_rank=graph.metadata.support_rank,
        visible_dim=graph.metadata.visible_dim,
        support_restricted=False,
        notes=(
            "sufficient graph-chart local frontier certificate",
            "can be vacuous when residuals or Lipschitz constants dominate the Hessian margin",
            "not a global Grassmannian optimizer or law-level branch probability",
        ),
    )
    return DeclaredFrontierLocalCertificateResult(
        graph_hessian=graph,
        mode=mode,
        eps=eps,
        lambda_margin=float(lambda_margin),
        lipschitz_bound=float(lipschitz),
        r0=float(r0),
        left_4eps_over_lambda=None if left is None else float(left),
        displacement_bound_2eps_over_lambda=None if displacement is None else float(displacement),
        certificate_passes=passes,
        certificate_kind=kind,
        chart_radius=radius,
        metadata=metadata,
    )


def _validate_operator_family(family: Sequence[np.ndarray] | np.ndarray, tolerances: Tolerances) -> list[np.ndarray]:
    if isinstance(family, np.ndarray):
        if family.ndim == 2:
            members = [family]
        elif family.ndim == 3:
            members = [family[idx] for idx in range(family.shape[0])]
        else:
            raise InputValidationError("family must be a symmetric matrix or stack of symmetric matrices")
    else:
        members = list(family)
    if not members:
        raise InputValidationError("family must contain at least one symmetric matrix")
    operators: list[np.ndarray] = []
    shape: tuple[int, int] | None = None
    for idx, member in enumerate(members):
        operator = validate_symmetric_matrix(f"family[{idx}]", member, tolerances)
        if shape is None:
            shape = operator.shape
        elif operator.shape != shape:
            raise InputValidationError("all family members must have the same shape")
        operators.append(operator)
    return operators


def _validate_stiefel_basis(
    B: np.ndarray,
    ambient_dim: int,
    tolerances: Tolerances,
    *,
    allow_full_rank: bool,
) -> np.ndarray:
    basis = to_float_array("B", B)
    if basis.shape[0] != ambient_dim:
        raise InputValidationError("B has incompatible ambient dimension")
    upper = ambient_dim if allow_full_rank else ambient_dim - 1
    if basis.shape[1] < 1 or basis.shape[1] > upper:
        raise InputValidationError("B has invalid visible rank")
    gram = symmetrize(basis.T @ basis)
    if not np.allclose(gram, np.eye(basis.shape[1]), atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError("B must have orthonormal columns")
    return basis


def _validate_or_make_complement(
    complement_basis: np.ndarray | None,
    basis: np.ndarray,
    tolerances: Tolerances,
) -> np.ndarray:
    if complement_basis is None:
        complement = np.asarray(la.null_space(basis.T, rcond=max(tolerances.atol, tolerances.rtol)), dtype=float)
    else:
        complement = to_float_array("complement_basis", complement_basis)
    if complement.shape != (basis.shape[0], basis.shape[0] - basis.shape[1]):
        raise InputValidationError("complement_basis has incompatible shape")
    gram = symmetrize(complement.T @ complement)
    if not np.allclose(gram, np.eye(complement.shape[1]), atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError("complement_basis must have orthonormal columns")
    cross = complement.T @ basis
    if not np.allclose(cross, np.zeros_like(cross), atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError("complement_basis must be orthogonal to B")
    return complement


def _validate_weights(weights: np.ndarray | None, count: int) -> np.ndarray:
    if weights is None:
        return np.ones(count, dtype=float)
    array = np.asarray(weights, dtype=float)
    if array.shape != (count,):
        raise InputValidationError("weights must have shape (family_count,)")
    if not np.all(np.isfinite(array)):
        raise InputValidationError("weights must be finite")
    if np.any(array < 0.0):
        raise InputValidationError("weights must be nonnegative")
    if float(np.sum(array)) <= 0.0:
        raise InputValidationError("at least one weight must be positive")
    return array.copy()


def _validate_mu(mu: float) -> float:
    if isinstance(mu, bool) or not np.isscalar(mu):
        raise InputValidationError("mu must be a nonnegative scalar")
    value = float(mu)
    if not np.isfinite(value) or value < 0.0:
        raise InputValidationError("mu must be a nonnegative scalar")
    return value


def _graph_second_variation_quadratic(
    x: np.ndarray,
    blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    weights: np.ndarray,
    mu: float,
) -> float:
    total = 0.0
    xtx = x.T @ x
    for weight, (u, d, e, k2_uu, k2_ww) in zip(weights, blocks):
        b = e.T @ x + x.T @ e
        s2 = (
            float(np.trace(b @ b))
            + 2.0 * float(np.trace(u @ x.T @ d @ x))
            - 2.0 * float(np.trace(xtx @ u @ u))
        )
        t2 = float(np.trace(x.T @ k2_ww @ x)) - float(np.trace(xtx @ k2_uu))
        total += 2.0 * weight * ((1.0 + mu) * s2 - mu * t2)
    return float(total)


def _frontier_lipschitz_bound(
    operators: list[np.ndarray],
    weights: np.ndarray,
    mu: float,
    visible_dim: int,
    ambient_dim: int,
    rho: float,
) -> float:
    y = float(np.sqrt(1.0 + rho * rho))
    c1 = 2.0 * y + 2.0 * rho * y * y
    c2 = 2.0 + 8.0 * rho * y + (8.0 * rho * rho + 2.0) * y * y
    c3 = (24.0 * rho + 48.0 * rho**3) * y * y + 6.0 * (8.0 * rho * rho + 2.0) * y + 12.0 * rho
    a2 = 0.0
    for weight, operator in zip(weights, operators):
        a2 += float(weight * np.linalg.norm(operator, ord=2) ** 2)
    g1 = (3.0 * (1.0 + mu) * np.sqrt(float(visible_dim)) + mu * np.sqrt(float(ambient_dim))) * a2
    g2 = 6.0 * (1.0 + mu) * a2
    g3 = 6.0 * (1.0 + mu) * a2
    return float(g3 * c1**3 + 3.0 * g2 * c2 * c1 + g1 * c3)


def _tangent_basis(hidden_dim: int, visible_dim: int) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for row in range(hidden_dim):
        for col in range(visible_dim):
            item = np.zeros((hidden_dim, visible_dim), dtype=float)
            item[row, col] = 1.0
            basis.append(item)
    return basis
