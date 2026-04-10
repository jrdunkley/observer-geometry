from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import scipy.linalg as la

from .exceptions import InputValidationError
from .types import (
    IntervalHessianResult,
    KernelJetResult,
    LinearAlgebraMetadata,
    LocalCoupledBirthResult,
    SampledIntervalLeakageResult,
    SemisimpleEventBlockResult,
    SupportRestartResult,
    SupportStratumTransportResult,
)
from .validation import (
    Tolerances,
    inv_sqrt_spd,
    inverse_spd,
    rank_cutoff,
    resolve_tolerances,
    solve_spd,
    sqrt_psd,
    support_decomposition_psd,
    symmetrize,
    to_float_array,
    validate_psd_matrix,
    validate_spd_matrix,
    validate_surjective_map,
    validate_symmetric_matrix,
)


def pi_from_hidden_load(
    Lambda: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return Pi = (I + Lambda)^(-1) for a reduced hidden load."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_psd_matrix("Lambda", Lambda, tol)
    if lambda_array.shape == (0, 0):
        return np.zeros((0, 0), dtype=float)
    return inverse_spd(np.eye(lambda_array.shape[0], dtype=float) + lambda_array)


def hidden_load_from_pi(
    Pi: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return Lambda = Pi^(-1) - I for a precision-side contraction variable."""
    tol = resolve_tolerances(tolerances)
    pi = validate_spd_matrix("Pi", Pi, tol)
    eigs = np.linalg.eigvalsh(pi)
    cutoff = 10.0 * rank_cutoff(eigs, tol)
    if float(np.max(eigs)) > 1.0 + cutoff:
        raise InputValidationError("Pi must satisfy Pi <= I on the active support")
    lambda_array = symmetrize(inverse_spd(pi) - np.eye(pi.shape[0], dtype=float))
    lambda_eigs = np.linalg.eigvalsh(lambda_array)
    if float(np.min(lambda_eigs)) < -rank_cutoff(lambda_eigs, tol):
        raise InputValidationError("Pi does not define a positive hidden load")
    return lambda_array


def pi_rhs(
    Pi: np.ndarray,
    A_cpl: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Evaluate dot(Pi) = -Pi^(1/2) A_cpl Pi^(1/2)."""
    tol = resolve_tolerances(tolerances)
    pi = validate_spd_matrix("Pi", Pi, tol)
    generator = validate_symmetric_matrix("A_cpl", A_cpl, tol)
    _require_same_shape("Pi", pi, "A_cpl", generator)
    root = sqrt_psd(pi, tol)
    return symmetrize(-root @ generator @ root)


def lambda_rhs(
    Lambda: np.ndarray,
    A_cpl: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Evaluate dot(Lambda) = (I + Lambda)^(1/2) A_cpl (I + Lambda)^(1/2)."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_psd_matrix("Lambda", Lambda, tol)
    generator = validate_symmetric_matrix("A_cpl", A_cpl, tol)
    _require_same_shape("Lambda", lambda_array, "A_cpl", generator)
    if lambda_array.shape == (0, 0):
        return np.zeros((0, 0), dtype=float)
    root = sqrt_psd(np.eye(lambda_array.shape[0], dtype=float) + lambda_array, tol)
    return symmetrize(root @ generator @ root)


def clock_rate(A_cpl: np.ndarray, tolerances: Tolerances | None = None) -> float:
    """Evaluate dot(tau) = Tr(A_cpl)."""
    tol = resolve_tolerances(tolerances)
    generator = validate_symmetric_matrix("A_cpl", A_cpl, tol)
    return float(np.trace(generator))


def support_stratum_transport(
    Lambda: np.ndarray,
    A_cpl: np.ndarray,
    tolerances: Tolerances | None = None,
    *,
    require_psd_generator: bool = False,
) -> SupportStratumTransportResult:
    """Evaluate the exact support-stable hidden-load transport diagnostics."""
    tol = resolve_tolerances(tolerances)
    lambda_array = validate_psd_matrix("Lambda", Lambda, tol)
    generator = validate_symmetric_matrix("A_cpl", A_cpl, tol)
    _require_same_shape("Lambda", lambda_array, "A_cpl", generator)
    if lambda_array.shape == (0, 0):
        eigs = np.zeros(0, dtype=float)
        generator_psd = True
        pi = np.zeros((0, 0), dtype=float)
        l_rhs = np.zeros((0, 0), dtype=float)
        p_rhs = np.zeros((0, 0), dtype=float)
    else:
        eigs = np.linalg.eigvalsh(generator)
        cutoff = rank_cutoff(eigs, tol)
        generator_psd = float(np.min(eigs)) >= -cutoff
        if require_psd_generator and not generator_psd:
            raise InputValidationError("A_cpl must be positive semidefinite for positive hidden-load transport")
        pi = pi_from_hidden_load(lambda_array, tolerances=tol)
        l_rhs = lambda_rhs(lambda_array, generator, tolerances=tol)
        p_rhs = pi_rhs(pi, generator, tolerances=tol)

    notes = () if generator_psd else ("A_cpl is not positive semidefinite; generator exists but positive hidden-load transport is not licensed",)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="support-stratum-transport",
        ambient_dim=lambda_array.shape[0],
        support_rank=lambda_array.shape[0],
        notes=notes,
    )
    return SupportStratumTransportResult(
        lambda_=lambda_array,
        pi=pi,
        generator=generator,
        lambda_rhs=l_rhs,
        pi_rhs=p_rhs,
        clock_rate=float(np.trace(generator)),
        a_min=float(np.min(eigs)) if eigs.size else None,
        a_max=float(np.max(eigs)) if eigs.size else None,
        generator_psd=generator_psd,
        metadata=metadata,
    )


def comparison_envelope_bounds(
    Pi0: np.ndarray,
    a_min_integral: float,
    a_max_integral: float,
    tolerances: Tolerances | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return lower and upper Pi envelopes from the stratumwise comparison theorem."""
    tol = resolve_tolerances(tolerances)
    pi0 = validate_spd_matrix("Pi0", Pi0, tol)
    if a_min_integral > a_max_integral:
        raise InputValidationError("a_min_integral must be <= a_max_integral")
    lower = np.exp(-float(a_max_integral)) * pi0
    upper = np.exp(-float(a_min_integral)) * pi0
    return symmetrize(lower), symmetrize(upper)


def restart_hidden_load_birth(
    lambda_before: np.ndarray,
    old_basis: np.ndarray,
    new_basis: np.ndarray,
    tolerances: Tolerances | None = None,
) -> SupportRestartResult:
    """Apply the forced finite positive birth restart in arbitrary nested bases."""
    tol = resolve_tolerances(tolerances)
    before = validate_psd_matrix("lambda_before", lambda_before, tol)
    old = _validate_orthonormal_basis("old_basis", old_basis, tol)
    new = _validate_orthonormal_basis("new_basis", new_basis, tol)
    if before.shape != (old.shape[1], old.shape[1]):
        raise InputValidationError("lambda_before must match old_basis dimension")
    if old.shape[0] != new.shape[0]:
        raise InputValidationError("old_basis and new_basis must have the same ambient dimension")
    if new.shape[1] < old.shape[1]:
        raise InputValidationError("birth restart requires new_basis to have at least old support rank")
    basis_map = new.T @ old
    _assert_nested_basis("old_basis", old, "new_basis", new, basis_map, tol)
    after = symmetrize(basis_map @ before @ basis_map.T)
    return _restart_result("birth", before, after, old, new, basis_map, tol)


def restart_hidden_load_death(
    lambda_before: np.ndarray,
    old_basis: np.ndarray,
    survivor_basis: np.ndarray,
    tolerances: Tolerances | None = None,
) -> SupportRestartResult:
    """Apply the forced finite positive death restart by survivor compression."""
    tol = resolve_tolerances(tolerances)
    before = validate_psd_matrix("lambda_before", lambda_before, tol)
    old = _validate_orthonormal_basis("old_basis", old_basis, tol)
    survivor = _validate_orthonormal_basis("survivor_basis", survivor_basis, tol)
    if before.shape != (old.shape[1], old.shape[1]):
        raise InputValidationError("lambda_before must match old_basis dimension")
    if old.shape[0] != survivor.shape[0]:
        raise InputValidationError("old_basis and survivor_basis must have the same ambient dimension")
    if survivor.shape[1] > old.shape[1]:
        raise InputValidationError("death restart requires survivor_basis rank <= old support rank")
    basis_map = old.T @ survivor
    _assert_nested_basis("survivor_basis", survivor, "old_basis", old, basis_map, tol)
    after = symmetrize(basis_map.T @ before @ basis_map)
    return _restart_result("death", before, after, old, survivor, basis_map, tol)


def kernel_schur_jet_from_coefficients(
    coefficients: Sequence[np.ndarray],
    max_order: int | None = None,
    tolerances: Tolerances | None = None,
) -> KernelJetResult:
    """Classify finite-order support events from Taylor coefficients of V(epsilon)."""
    tol = resolve_tolerances(tolerances)
    coeffs = _validate_coefficients(coefficients, tol)
    if max_order is None:
        max_order = len(coeffs) - 1
    if max_order < 1:
        raise InputValidationError("max_order must be at least 1")
    max_order = min(max_order, len(coeffs) - 1)

    v0 = coeffs[0]
    decomposition = support_decomposition_psd(v0, tol)
    gap_basis = decomposition.basis
    kernel_basis = decomposition.complement
    kernel_dim = kernel_basis.shape[1]
    gap_dim = gap_basis.shape[1]
    if kernel_dim == 0:
        return _kernel_jet_empty(coeffs, gap_basis, np.zeros((0, 0), dtype=float), tol, "none")

    gap_operator = symmetrize(gap_basis.T @ v0 @ gap_basis) if gap_dim else np.zeros((0, 0), dtype=float)
    a_coeffs = [np.zeros((kernel_dim, kernel_dim), dtype=float)]
    b_coeffs = [np.zeros((kernel_dim, gap_dim), dtype=float)]
    c_coeffs = [np.zeros((gap_dim, gap_dim), dtype=float)]
    for idx in range(1, max_order + 1):
        coeff = coeffs[idx]
        a_coeffs.append(symmetrize(kernel_basis.T @ coeff @ kernel_basis))
        b_coeffs.append(kernel_basis.T @ coeff @ gap_basis if gap_dim else np.zeros((kernel_dim, 0), dtype=float))
        c_coeffs.append(symmetrize(gap_basis.T @ coeff @ gap_basis) if gap_dim else np.zeros((0, 0), dtype=float))

    if gap_dim:
        p_inv = inverse_spd(gap_operator)
        d_coeffs = [p_inv]
        for n in range(1, max_order + 1):
            total = np.zeros_like(gap_operator)
            for k in range(1, n + 1):
                total = total + c_coeffs[k] @ d_coeffs[n - k]
            d_coeffs.append(-p_inv @ total)
    else:
        d_coeffs = [np.zeros((0, 0), dtype=float) for _ in range(max_order + 1)]

    effective: list[np.ndarray] = []
    for n in range(1, max_order + 1):
        e_n = a_coeffs[n].copy()
        if gap_dim:
            for p in range(1, n + 1):
                for r in range(1, n + 1):
                    q = n - p - r
                    if q >= 0:
                        e_n = e_n - b_coeffs[p] @ d_coeffs[q] @ b_coeffs[r].T
        effective.append(symmetrize(e_n))

    return _kernel_jet_from_effective(coeffs, kernel_basis, gap_basis, gap_operator, tuple(effective), tol)


def classify_support_event_from_jet(
    jet: KernelJetResult,
    side_sign: int = 1,
    tolerances: Tolerances | None = None,
) -> str:
    """Classify the event type on a chosen side from a computed kernel jet."""
    tol = resolve_tolerances(tolerances)
    _validate_side_sign(side_sign)
    if jet.order is None or jet.leading_eigenvalues.size == 0:
        return "none" if jet.event_kind == "none" else "degenerate"
    scaled = (float(side_sign) ** jet.order) * jet.leading_eigenvalues
    cutoff = rank_cutoff(scaled, tol)
    positive = int(np.sum(scaled > cutoff))
    negative = int(np.sum(scaled < -cutoff))
    zero = scaled.size - positive - negative
    if zero:
        return "degenerate"
    if positive and negative:
        return "mixed"
    if positive:
        if jet.order % 2 == 0:
            return "touch"
        return "birth" if side_sign > 0 else "death"
    if negative:
        if jet.order % 2 == 0:
            return "degenerate"
        return "death" if side_sign > 0 else "birth"
    return "degenerate"


def semisimple_event_block(
    jet: KernelJetResult,
    block_indices: Sequence[int] | None = None,
    side_sign: int = 1,
    tolerances: Tolerances | None = None,
) -> SemisimpleEventBlockResult:
    """Return universal semisimple pole and clock diagnostics for a kernel-jet block."""
    tol = resolve_tolerances(tolerances)
    _validate_side_sign(side_sign)
    if jet.order is None or jet.leading_effective.shape == (0, 0):
        raise InputValidationError("jet has no leading event block")
    eigvals, eigvecs = np.linalg.eigh(jet.leading_effective)
    if block_indices is None:
        cutoff = rank_cutoff(eigvals, tol)
        indices = [idx for idx, value in enumerate(eigvals) if abs(float(value)) > cutoff]
    else:
        indices = [int(idx) for idx in block_indices]
    if not indices:
        raise InputValidationError("event block must contain at least one nonzero direction")
    if min(indices) < 0 or max(indices) >= eigvals.size:
        raise InputValidationError("block_indices are out of range")

    block_basis = eigvecs[:, indices]
    leading = symmetrize(block_basis.T @ jet.leading_effective @ block_basis)
    active_leading = (float(side_sign) ** jet.order) * leading
    active_eigs = np.linalg.eigvalsh(active_leading)
    cutoff = rank_cutoff(active_eigs, tol)
    if float(np.min(active_eigs)) <= cutoff:
        raise InputValidationError("selected block is not active positive definite on the chosen side")

    pole = -float(side_sign) * float(jet.order) / 2.0
    death_like = pole > 0.0
    birth_like = pole < 0.0
    dimension = len(indices)
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=cutoff,
        method="semisimple-event-block",
        ambient_dim=jet.kernel_basis.shape[0],
        support_rank=dimension,
        notes=() if death_like else ("birth-like block fails positive hidden-load transport near the event",),
    )
    return SemisimpleEventBlockResult(
        order=int(jet.order),
        side_sign=int(side_sign),
        dimension=dimension,
        leading_matrix=active_leading,
        pole_coefficient=pole,
        death_like=death_like,
        birth_like=birth_like,
        clock_log_coefficient=(float(jet.order) * dimension / 2.0) if death_like else 0.0,
        desingularisation_power=(float(jet.order) / 2.0) if death_like else 0.0,
        metadata=metadata,
    )


def default_hidden_complement(
    C: np.ndarray,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return an orthonormal hidden complement basis for the kernel of C."""
    tol = resolve_tolerances(tolerances)
    c = to_float_array("C", C)
    if c.ndim != 2:
        raise InputValidationError("C must be a 2D array")
    return np.asarray(la.null_space(c, rcond=max(tol.atol, tol.rtol)), dtype=float)


def local_coupled_birth(
    H: np.ndarray,
    H_dot: np.ndarray,
    H_ddot: np.ndarray,
    C: np.ndarray,
    C_dot: np.ndarray,
    Z: np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> LocalCoupledBirthResult:
    """Extract W = D_t^alpha V and A_cpl from explicit path derivative data."""
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    h_dot = validate_symmetric_matrix("H_dot", H_dot, tol)
    h_ddot = validate_symmetric_matrix("H_ddot", H_ddot, tol)
    _require_same_shape("H", h, "H_dot", h_dot)
    _require_same_shape("H", h, "H_ddot", h_ddot)
    c = validate_surjective_map(C, h.shape[0], tol)
    c_dot = to_float_array("C_dot", C_dot)
    if c_dot.shape != c.shape:
        raise InputValidationError("C_dot must have the same shape as C")
    hidden_basis = default_hidden_complement(c, tolerances=tol) if Z is None else _validate_hidden_complement(c, Z, tol)

    phi = inverse_spd(c @ inverse_spd(h) @ c.T)
    lift = solve_spd(h, c.T @ phi)
    hidden_metric = symmetrize(hidden_basis.T @ h @ hidden_basis)
    if hidden_metric.shape != (0, 0):
        hidden_metric = validate_spd_matrix("R", hidden_metric, tol)
    V = symmetrize(lift.T @ h_dot @ lift)
    B = lift.T @ h_dot @ hidden_basis
    beta = -c_dot @ hidden_basis
    if hidden_metric.shape == (0, 0):
        Q = np.zeros_like(V)
        observer_tensor = np.zeros_like(V)
        mixed_left = np.zeros_like(V)
        mixed_right = np.zeros_like(V)
    else:
        r_inv_b_t = solve_spd(hidden_metric, B.T)
        r_inv_beta_t_phi = solve_spd(hidden_metric, beta.T @ phi)
        Q = symmetrize(B @ r_inv_b_t)
        observer_tensor = symmetrize(phi @ beta @ solve_spd(hidden_metric, beta.T @ phi))
        mixed_left = phi @ beta @ r_inv_b_t
        mixed_right = B @ r_inv_beta_t_phi
    W = symmetrize(lift.T @ h_ddot @ lift - 2.0 * Q - mixed_left - mixed_right)

    decomposition = support_decomposition_psd(V, tol)
    active = decomposition.basis
    if active.shape[1] == 0:
        a_cpl = np.zeros((0, 0), dtype=float)
    else:
        v_s = symmetrize(active.T @ V @ active)
        w_s = symmetrize(active.T @ W @ active)
        inv_root = inv_sqrt_spd(v_s, tol)
        a_cpl = symmetrize(-0.5 * inv_root @ w_s @ inv_root)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=decomposition.cutoff,
        method="local-coupled-birth-extractor",
        ambient_dim=h.shape[0],
        support_rank=active.shape[1],
        visible_dim=c.shape[0],
        condition_number=_condition_number_spd(h),
        notes=_condition_notes(h),
    )
    return LocalCoupledBirthResult(
        phi=phi,
        lift=lift,
        hidden_basis=hidden_basis,
        hidden_metric=hidden_metric,
        V=V,
        B=B,
        beta=beta,
        Q=Q,
        observer_tensor=observer_tensor,
        W=W,
        active_support_basis=active,
        a_cpl=a_cpl,
        metadata=metadata,
    )


def sampled_interval_leakage(
    family: Sequence[np.ndarray] | np.ndarray,
    projector_or_basis: np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> SampledIntervalLeakageResult:
    """Evaluate sampled interval leakage and visibility for a projector or basis."""
    tol = resolve_tolerances(tolerances)
    members = _validate_family(family, tol)
    projector, rank = _projector_from_basis_or_projector(projector_or_basis, members[0].shape[0], tol)
    sample_weights = _validate_weights(weights, len(members))
    leakage, visible, stationarity = _sampled_interval_components(members, projector, sample_weights)
    residual = float(np.linalg.norm(projector @ stationarity - stationarity @ projector, ord="fro"))
    cutoff = 100.0 * max(tol.atol, tol.rtol * max(1.0, leakage, visible))
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="sampled-interval-leakage",
        ambient_dim=members[0].shape[0],
        support_rank=rank,
        visible_dim=rank,
        notes=("sampled exact closure only; not a continuum certificate",),
    )
    return SampledIntervalLeakageResult(
        leakage=leakage,
        visible_score=visible,
        stationarity_residual=residual,
        projector=projector,
        sample_count=len(members),
        weights=sample_weights,
        sampled_exact_closure=leakage <= cutoff,
        metadata=metadata,
    )


def sampled_interval_stationarity(
    family: Sequence[np.ndarray] | np.ndarray,
    projector_or_basis: np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> np.ndarray:
    """Return [P, sum_k w_k [A_k, [A_k, P]]] for a sampled interval family."""
    tol = resolve_tolerances(tolerances)
    members = _validate_family(family, tol)
    projector, _rank = _projector_from_basis_or_projector(projector_or_basis, members[0].shape[0], tol)
    sample_weights = _validate_weights(weights, len(members))
    _leakage, _visible, stationarity = _sampled_interval_components(members, projector, sample_weights)
    commutator = projector @ stationarity - stationarity @ projector
    return 0.5 * (commutator - commutator.T)


def sampled_interval_closure_check(
    family: Sequence[np.ndarray] | np.ndarray,
    projector_or_basis: np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> bool:
    """Return whether the supplied samples commute with the projector within tolerance."""
    return sampled_interval_leakage(family, projector_or_basis, weights=weights, tolerances=tolerances).sampled_exact_closure


def interval_hessian_at_exact_family(
    family: Sequence[np.ndarray] | np.ndarray,
    basis: np.ndarray,
    weights: Sequence[float] | np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> IntervalHessianResult:
    """Return the sampled Hessian diagnostic at an exact interval observer basis."""
    tol = resolve_tolerances(tolerances)
    members = _validate_family(family, tol)
    visible_basis = _validate_orthonormal_basis("basis", basis, tol)
    if visible_basis.shape[0] != members[0].shape[0]:
        raise InputValidationError("basis has incompatible ambient dimension")
    sample_weights = _validate_weights(weights, len(members))
    hidden_basis = np.asarray(la.null_space(visible_basis.T, rcond=max(tol.atol, tol.rtol)), dtype=float)
    hidden_dim = hidden_basis.shape[1]
    visible_dim = visible_basis.shape[1]
    hessian = np.zeros((hidden_dim * visible_dim, hidden_dim * visible_dim), dtype=float)
    spectral_gap: float | None = None
    for weight, member in zip(sample_weights, members):
        a_u = symmetrize(visible_basis.T @ member @ visible_basis)
        a_perp = symmetrize(hidden_basis.T @ member @ hidden_basis)
        if hidden_dim:
            eig_u = np.linalg.eigvalsh(a_u)
            eig_perp = np.linalg.eigvalsh(a_perp)
            gap = float(np.min(np.abs(eig_perp[:, None] - eig_u[None, :])))
            spectral_gap = gap if spectral_gap is None else min(spectral_gap, gap)
            operator = np.kron(np.eye(visible_dim), a_perp) - np.kron(a_u.T, np.eye(hidden_dim))
            hessian = hessian + 2.0 * float(weight) * (operator.T @ operator)
    if hidden_dim == 0:
        spectral_gap = None
        lower_bound = None
        locally_rigid = True
    else:
        hessian_eigs = np.linalg.eigvalsh(symmetrize(hessian))
        cutoff = rank_cutoff(hessian_eigs, tol)
        locally_rigid = float(np.min(hessian_eigs)) > cutoff
        lower_bound = None if spectral_gap is None else 2.0 * spectral_gap * spectral_gap * float(np.sum(sample_weights))
    metadata = LinearAlgebraMetadata(
        atol=tol.atol,
        rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="sampled-interval-hessian",
        ambient_dim=members[0].shape[0],
        support_rank=visible_dim,
        visible_dim=visible_dim,
        notes=("sampled Hessian diagnostic; assumes supplied basis is exact for the sampled family",),
    )
    return IntervalHessianResult(
        hessian_operator=symmetrize(hessian),
        spectral_gap=spectral_gap,
        rigidity_lower_bound=lower_bound,
        locally_rigid=locally_rigid,
        metadata=metadata,
    )


def _require_same_shape(left_name: str, left: np.ndarray, right_name: str, right: np.ndarray) -> None:
    if left.shape != right.shape:
        raise InputValidationError(f"{left_name} and {right_name} must have the same shape")


def _validate_orthonormal_basis(name: str, basis: np.ndarray, tolerances: Tolerances) -> np.ndarray:
    array = to_float_array(name, basis)
    gram = array.T @ array
    if not np.allclose(gram, np.eye(array.shape[1], dtype=float), atol=tolerances.atol, rtol=tolerances.rtol):
        raise InputValidationError(f"{name} must have orthonormal columns")
    return array


def _validate_hidden_complement(C: np.ndarray, Z: np.ndarray, tolerances: Tolerances) -> np.ndarray:
    z = to_float_array("Z", Z)
    expected_dim = C.shape[1] - C.shape[0]
    if z.shape != (C.shape[1], expected_dim):
        raise InputValidationError("Z must have shape (ambient_dim, ambient_dim - visible_dim)")
    if z.shape[1] == 0:
        return z
    singular_values = np.linalg.svd(z, compute_uv=False)
    cutoff = rank_cutoff(singular_values, tolerances)
    if int(np.sum(singular_values > cutoff)) != z.shape[1]:
        raise InputValidationError("Z must have full column rank")
    if not np.allclose(C @ z, 0.0, atol=10.0 * tolerances.atol, rtol=10.0 * tolerances.rtol):
        raise InputValidationError("Z must be a hidden complement satisfying C @ Z = 0")
    return z


def _assert_nested_basis(
    nested_name: str,
    nested: np.ndarray,
    container_name: str,
    container: np.ndarray,
    basis_map: np.ndarray,
    tolerances: Tolerances,
) -> None:
    residual = float(np.linalg.norm(nested - container @ basis_map, ord="fro"))
    cutoff = 100.0 * max(tolerances.atol, tolerances.rtol)
    if residual > cutoff:
        raise InputValidationError(f"{nested_name} must lie in the span of {container_name}")


def _restart_result(
    event_kind: str,
    before: np.ndarray,
    after: np.ndarray,
    old_basis: np.ndarray,
    new_basis: np.ndarray,
    basis_map: np.ndarray,
    tolerances: Tolerances,
) -> SupportRestartResult:
    after = validate_psd_matrix("lambda_after", after, tolerances)
    metadata = LinearAlgebraMetadata(
        atol=tolerances.atol,
        rtol=tolerances.rtol,
        rank_tol=max(tolerances.atol, tolerances.rtol),
        method=f"{event_kind}-restart",
        ambient_dim=old_basis.shape[0],
        support_rank=new_basis.shape[1],
        notes=("forced finite positive restart",),
    )
    return SupportRestartResult(
        event_kind=event_kind,
        lambda_before=before,
        lambda_after=after,
        old_basis=old_basis,
        new_basis=new_basis,
        basis_map=basis_map,
        metadata=metadata,
    )


def _validate_coefficients(coefficients: Sequence[np.ndarray], tolerances: Tolerances) -> list[np.ndarray]:
    if isinstance(coefficients, np.ndarray):
        raise InputValidationError("coefficients must be a sequence of square matrices")
    coeffs = [validate_symmetric_matrix(f"coefficients[{idx}]", matrix, tolerances) for idx, matrix in enumerate(coefficients)]
    if len(coeffs) < 2:
        raise InputValidationError("coefficients must include V0 and at least one higher-order coefficient")
    shape = coeffs[0].shape
    for idx, coeff in enumerate(coeffs):
        if coeff.shape != shape:
            raise InputValidationError(f"coefficients[{idx}] must have the same shape as coefficients[0]")
    validate_psd_matrix("coefficients[0]", coeffs[0], tolerances)
    return coeffs


def _kernel_jet_empty(
    coeffs: Sequence[np.ndarray],
    gap_basis: np.ndarray,
    gap_operator: np.ndarray,
    tolerances: Tolerances,
    event_kind: str,
) -> KernelJetResult:
    metadata = LinearAlgebraMetadata(
        atol=tolerances.atol,
        rtol=tolerances.rtol,
        rank_tol=max(tolerances.atol, tolerances.rtol),
        method="kernel-schur-jet",
        ambient_dim=coeffs[0].shape[0],
        support_rank=0,
        notes=("no kernel event",) if event_kind == "none" else ("no resolved kernel event in supplied coefficients",),
    )
    return KernelJetResult(
        order=None,
        kernel_basis=np.zeros((coeffs[0].shape[0], 0), dtype=float),
        gap_basis=gap_basis,
        gap_operator=gap_operator,
        effective_coefficients=(),
        leading_effective=np.zeros((0, 0), dtype=float),
        leading_eigenvalues=np.zeros(0, dtype=float),
        birth_count_forward=0,
        death_count_forward=0,
        zero_count=0,
        event_kind=event_kind,
        semisimple=False,
        metadata=metadata,
    )


def _kernel_jet_from_effective(
    coeffs: Sequence[np.ndarray],
    kernel_basis: np.ndarray,
    gap_basis: np.ndarray,
    gap_operator: np.ndarray,
    effective: tuple[np.ndarray, ...],
    tolerances: Tolerances,
) -> KernelJetResult:
    order = None
    leading = np.zeros((kernel_basis.shape[1], kernel_basis.shape[1]), dtype=float)
    for idx, coeff in enumerate(effective, start=1):
        eigs = np.linalg.eigvalsh(coeff)
        if float(np.max(np.abs(eigs))) > rank_cutoff(eigs, tolerances):
            order = idx
            leading = coeff
            break
    if order is None:
        result = _kernel_jet_empty(coeffs, gap_basis, gap_operator, tolerances, "degenerate")
        return KernelJetResult(
            order=None,
            kernel_basis=kernel_basis,
            gap_basis=gap_basis,
            gap_operator=gap_operator,
            effective_coefficients=effective,
            leading_effective=result.leading_effective,
            leading_eigenvalues=result.leading_eigenvalues,
            birth_count_forward=0,
            death_count_forward=0,
            zero_count=kernel_basis.shape[1],
            event_kind="degenerate",
            semisimple=False,
            metadata=result.metadata,
        )

    eigs = np.linalg.eigvalsh(leading)
    cutoff = rank_cutoff(eigs, tolerances)
    birth = int(np.sum(eigs > cutoff))
    death = int(np.sum(eigs < -cutoff))
    zero = int(eigs.size - birth - death)
    if zero:
        event_kind = "degenerate"
    elif birth and death:
        event_kind = "mixed"
    elif birth:
        event_kind = "touch" if order % 2 == 0 else "birth"
    elif death:
        event_kind = "degenerate" if order % 2 == 0 else "death"
    else:
        event_kind = "degenerate"
    metadata = LinearAlgebraMetadata(
        atol=tolerances.atol,
        rtol=tolerances.rtol,
        rank_tol=cutoff,
        method="kernel-schur-jet",
        ambient_dim=coeffs[0].shape[0],
        support_rank=kernel_basis.shape[1],
        notes=("leading small-eigenvalue jets, not full-spectrum equality",),
    )
    return KernelJetResult(
        order=order,
        kernel_basis=kernel_basis,
        gap_basis=gap_basis,
        gap_operator=gap_operator,
        effective_coefficients=effective,
        leading_effective=leading,
        leading_eigenvalues=eigs,
        birth_count_forward=birth,
        death_count_forward=death,
        zero_count=zero,
        event_kind=event_kind,
        semisimple=zero == 0,
        metadata=metadata,
    )


def _validate_side_sign(side_sign: int) -> None:
    if side_sign not in {-1, 1}:
        raise InputValidationError("side_sign must be -1 or 1")


def _validate_family(family: Sequence[np.ndarray] | np.ndarray, tolerances: Tolerances) -> list[np.ndarray]:
    if isinstance(family, np.ndarray):
        members = [family]
    else:
        members = list(family)
    if not members:
        raise InputValidationError("family must contain at least one symmetric matrix")
    validated = [validate_symmetric_matrix(f"family[{idx}]", member, tolerances) for idx, member in enumerate(members)]
    shape = validated[0].shape
    for idx, member in enumerate(validated):
        if member.shape != shape:
            raise InputValidationError(f"family[{idx}] must have the same shape as family[0]")
    return validated


def _validate_weights(weights: Sequence[float] | np.ndarray | None, count: int) -> np.ndarray:
    if weights is None:
        return np.ones(count, dtype=float)
    sample_weights = np.asarray(weights, dtype=float)
    if sample_weights.ndim != 1 or sample_weights.shape[0] != count:
        raise InputValidationError("weights must be a one-dimensional array matching the family length")
    if np.any(sample_weights <= 0.0):
        raise InputValidationError("weights must be strictly positive")
    return sample_weights


def _projector_from_basis_or_projector(
    projector_or_basis: np.ndarray,
    ambient_dim: int,
    tolerances: Tolerances,
) -> tuple[np.ndarray, int]:
    array = to_float_array("projector_or_basis", projector_or_basis)
    if array.shape[0] != ambient_dim:
        raise InputValidationError("projector_or_basis has incompatible ambient dimension")
    if array.shape[1] == ambient_dim and np.allclose(array, array.T, atol=tolerances.atol, rtol=tolerances.rtol):
        projector = symmetrize(array)
        if not np.allclose(projector @ projector, projector, atol=10.0 * tolerances.atol, rtol=10.0 * tolerances.rtol):
            raise InputValidationError("projector_or_basis square input must be an orthogonal projector")
        eigs = np.linalg.eigvalsh(projector)
        rank = int(np.sum(eigs > 0.5))
        return projector, rank
    basis = _validate_orthonormal_basis("projector_or_basis", array, tolerances)
    return symmetrize(basis @ basis.T), basis.shape[1]


def _sampled_interval_components(
    family: Sequence[np.ndarray],
    projector: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    leakage = 0.0
    visible = 0.0
    stationarity = np.zeros_like(projector)
    for weight, member in zip(weights, family):
        commutator = member @ projector - projector @ member
        leakage += float(weight) * 0.5 * float(np.linalg.norm(commutator, ord="fro") ** 2)
        compressed = projector @ member @ projector
        visible += float(weight) * float(np.linalg.norm(compressed, ord="fro") ** 2)
        stationarity = stationarity + float(weight) * (member @ commutator - commutator @ member)
    return leakage, visible, symmetrize(stationarity)


def _condition_number_spd(H: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(H)
    return float(np.max(eigenvalues) / np.min(eigenvalues))


def _condition_notes(H: np.ndarray) -> tuple[str, ...]:
    condition = _condition_number_spd(H)
    if condition >= 1e9:
        return (f"H is extremely ill-conditioned (cond={condition:.3e}); local coupled extraction may be numerically unreliable",)
    if condition >= 1e8:
        return (f"H is severely ill-conditioned (cond={condition:.3e}); interpret local coupled extraction with care",)
    if condition >= 1e6:
        return (f"H is highly conditioned (cond={condition:.3e}); local coupled diagnostics should be checked",)
    return ()
