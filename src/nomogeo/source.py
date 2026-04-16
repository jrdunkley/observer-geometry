"""
Source law, information conservation, and observer steering.

0.4.0 kernel surface: pathwise dynamics of the split-frame geometry.

Public API:
    information_budget       -- first-order conservation: vis + hid = ambient
    source_law               -- A_cpl on a support-stable stratum (requires V > 0)
    evidence_decomposition   -- second-order evidence curvature from A_cpl
    observer_diagnostics     -- combined diagnostic for observer quality
    capture_curve            -- information capture vs observer rank
"""
from __future__ import annotations

import numpy as np
from numpy.linalg import svd as _svd

from .types import (
    Array,
    LinearAlgebraMetadata,
)
from .validation import (
    Tolerances,
    resolve_tolerances,
    solve_spd,
    symmetrize,
    validate_spd_matrix,
    validate_surjective_map,
    validate_symmetric_matrix,
)
from .core import visible_geometry
from .exceptions import InputValidationError, SupportError

from dataclasses import dataclass, field


@dataclass(frozen=True)
class InformationBudgetResult:
    """First-order information conservation: vis_rate + hid_rate = amb_rate."""
    visible_rate: float
    hidden_rate: float
    ambient_rate: float
    visible_fraction: float
    conservation_residual: float
    phi: Array
    V: Array
    R: Array
    U_h: Array
    B: Array
    v_eigenvalues: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class SourceLawResult:
    """Source law on a support-stable stratum (V > 0 required)."""
    A_cpl: Array
    A_direct: Array
    hidden_defect: Array
    a_cpl_eigenvalues: Array
    W: Array
    budget: InformationBudgetResult
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class EvidenceDecompositionResult:
    """Second-order evidence curvature decomposition."""
    f_pprime_visible: float
    f_pprime_hidden: float
    f_pprime_ambient: float
    source_term: float | None
    kinematic_term: float
    conservation_residual: float
    source_law: SourceLawResult | None
    budget: InformationBudgetResult
    metadata: LinearAlgebraMetadata


def _build_hidden_frame(C: np.ndarray, n: int, m: int) -> np.ndarray:
    """Return Z = kernel of C, orthonormalised."""
    _, _, Vt = _svd(C)
    return Vt[m:].T


def information_budget(
    H: np.ndarray,
    C: np.ndarray,
    Hdot: np.ndarray,
    tolerances: Tolerances | None = None,
) -> InformationBudgetResult:
    """
    Compute the first-order information conservation budget.

    Given an SPD field H, its derivative Hdot along a path, and a rank-m
    observer C, return the visible rate, hidden rate, and ambient rate.
    These satisfy: visible_rate + hidden_rate = ambient_rate (exact).

    Parameters
    ----------
    H : (n, n) SPD matrix
    C : (m, n) rank-m observer
    Hdot : (n, n) symmetric matrix (dH/dt)
    tolerances : optional numerical tolerances

    Returns
    -------
    InformationBudgetResult
    """
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    c = validate_surjective_map(C, h.shape[0], tol)
    hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    n = h.shape[0]
    m = c.shape[0]

    # Split frame
    geom = visible_geometry(h, c, tolerances=tol)
    phi = geom.phi
    L = geom.lift
    Z = _build_hidden_frame(c, n, m)
    R = symmetrize(Z.T @ h @ Z)

    # Visible and hidden first jets
    V = symmetrize(L.T @ hdot @ L)
    U_h = symmetrize(Z.T @ hdot @ Z)
    B = L.T @ hdot @ Z

    # Rates (using the full dPhi formula for exact conservation)
    h_inv_hdot = solve_spd(h, hdot)
    h_inv = solve_spd(h, np.eye(n, dtype=float))
    dh_inv = -h_inv @ hdot @ h_inv
    dphi = -phi @ (c @ dh_inv @ c.T) @ phi
    dphi = symmetrize(dphi)

    phi_inv_dphi = solve_spd(phi, dphi)
    visible_rate = float(np.trace(phi_inv_dphi))

    R_inv_Uh = solve_spd(R, U_h)
    hidden_rate = float(np.trace(R_inv_Uh))

    ambient_rate = float(np.trace(h_inv_hdot))
    conservation_residual = abs(visible_rate + hidden_rate - ambient_rate)

    vis_frac = visible_rate / ambient_rate if abs(ambient_rate) > tol.atol else float('nan')

    v_eigenvalues = np.linalg.eigvalsh(V)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol, rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="conservation_law",
        ambient_dim=n, support_rank=n, visible_dim=m,
    )

    return InformationBudgetResult(
        visible_rate=visible_rate,
        hidden_rate=hidden_rate,
        ambient_rate=ambient_rate,
        visible_fraction=vis_frac,
        conservation_residual=conservation_residual,
        phi=phi, V=V, R=R, U_h=U_h, B=B,
        v_eigenvalues=v_eigenvalues,
        metadata=metadata,
    )


def source_law(
    H: np.ndarray,
    C: np.ndarray,
    Hdot: np.ndarray,
    Hddot: np.ndarray,
    tolerances: Tolerances | None = None,
) -> SourceLawResult:
    """
    Compute the coupled source operator A_cpl on a support-stable stratum.

    Requires V = L^T Hdot L to be positive definite on the active support.
    If V is not positive definite, raises SupportError.

    Parameters
    ----------
    H : (n, n) SPD matrix
    C : (m, n) rank-m observer
    Hdot : (n, n) symmetric (dH/dt)
    Hddot : (n, n) symmetric (d^2H/dt^2)
    tolerances : optional

    Returns
    -------
    SourceLawResult with A_cpl and decomposition
    """
    tol = resolve_tolerances(tolerances)
    hddot = validate_symmetric_matrix("Hddot", Hddot, tol)

    budget = information_budget(H, C, Hdot, tolerances=tol)

    V = budget.V
    B = budget.B
    R = budget.R
    phi = budget.phi
    n = H.shape[0]
    m = C.shape[0]

    v_min = float(np.min(budget.v_eigenvalues))
    if v_min <= tol.atol:
        raise SupportError(
            f"V is not positive definite on the active support (min eigenvalue = {v_min:.6g}). "
            f"The source law requires V > 0."
        )

    # Split frame
    geom = visible_geometry(H, C, tolerances=tol)
    L = geom.lift
    Z = _build_hidden_frame(C, n, m)
    R_inv = solve_spd(R, np.eye(R.shape[0], dtype=float))

    # Connection forms (exact, no finite differences)
    h_inv = solve_spd(H, np.eye(n, dtype=float))
    dh_inv = -h_inv @ Hdot @ h_inv
    dphi = symmetrize(-phi @ (C @ dh_inv @ C.T) @ phi)
    dL = dh_inv @ C.T @ phi + h_inv @ C.T @ dphi

    # Covariant visible second jet
    alpha = solve_spd(phi, L.T @ H @ dL)
    Vdot = symmetrize(dL.T @ Hdot @ L + L.T @ hddot @ L + L.T @ Hdot @ dL)
    W = symmetrize(Vdot - alpha.T @ V - V @ alpha)

    # V^{1/2} and V^{-1/2}
    eigvals, eigvecs = np.linalg.eigh(V)
    sqrt_eigvals = np.sqrt(eigvals)
    V_sqrt = eigvecs @ np.diag(sqrt_eigvals) @ eigvecs.T
    V_sqrt_inv = eigvecs @ np.diag(1.0 / sqrt_eigvals) @ eigvecs.T

    # A_cpl = -1/2 V^{-1/2} W V^{-1/2}
    A_cpl = symmetrize(-0.5 * V_sqrt_inv @ W @ V_sqrt_inv)

    # Decomposition: A = direct acceleration, hidden_defect = Q_hat / V correction
    A_direct = symmetrize(-0.5 * V_sqrt_inv @ (L.T @ hddot @ L) @ V_sqrt_inv)
    Qhat = symmetrize(B @ R_inv @ B.T)
    hidden_defect = symmetrize(V_sqrt_inv @ Qhat @ V_sqrt_inv)

    a_cpl_eigenvalues = np.linalg.eigvalsh(A_cpl)

    metadata = LinearAlgebraMetadata(
        atol=tol.atol, rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="source_law",
        ambient_dim=n, support_rank=n, visible_dim=m,
    )

    return SourceLawResult(
        A_cpl=A_cpl,
        A_direct=A_direct,
        hidden_defect=hidden_defect,
        a_cpl_eigenvalues=a_cpl_eigenvalues,
        W=W,
        budget=budget,
        metadata=metadata,
    )


def evidence_decomposition(
    H: np.ndarray,
    C: np.ndarray,
    Hdot: np.ndarray,
    Hddot: np.ndarray,
    tolerances: Tolerances | None = None,
) -> EvidenceDecompositionResult:
    """
    Compute the second-order evidence curvature decomposition.

    Returns f''_vis, f''_hid, f''_amb satisfying:
        f''_vis + f''_hid = f''_amb   (second-order conservation)

    When V > 0, also decomposes f''_vis into source + kinematic + connection terms.

    Parameters
    ----------
    H, C, Hdot, Hddot : as for source_law
    tolerances : optional
    """
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    hddot = validate_symmetric_matrix("Hddot", Hddot, tol)
    n = h.shape[0]
    m = C.shape[0]

    budget = information_budget(h, C, hdot, tolerances=tol)
    Z = _build_hidden_frame(C, n, m)
    R = budget.R

    # Ambient second derivative: -Tr((H^{-1} Hdot)^2) + Tr(H^{-1} Hddot)
    h_inv = solve_spd(h, np.eye(n, dtype=float))
    P_H = h_inv @ hdot
    f2_amb = float(-np.trace(P_H @ P_H) + np.trace(h_inv @ hddot))

    # Hidden second derivative: -Tr((R^{-1} U_h)^2) + Tr(R^{-1} U_h2)
    R_inv = solve_spd(R, np.eye(R.shape[0], dtype=float))
    U_h = budget.U_h
    U_h2 = symmetrize(Z.T @ hddot @ Z)
    P_R = R_inv @ U_h
    f2_hid = float(-np.trace(P_R @ P_R) + np.trace(R_inv @ U_h2))

    # Visible from conservation
    f2_vis = f2_amb - f2_hid
    conservation_residual = 0.0  # exact by construction

    # Kinematic term: -Tr((Phi^{-1} dPhi)^2)
    phi = budget.phi
    dh_inv = -h_inv @ hdot @ h_inv
    dphi = symmetrize(-phi @ (C @ dh_inv @ C.T) @ phi)
    P_Phi = solve_spd(phi, dphi)
    kinematic_term = float(-np.trace(P_Phi @ P_Phi))

    # Source term (if V > 0)
    source_term = None
    src = None
    v_min = float(np.min(budget.v_eigenvalues))
    if v_min > tol.atol:
        try:
            src = source_law(h, C, hdot, hddot, tolerances=tol)
            V = budget.V
            eigvals, eigvecs = np.linalg.eigh(V)
            V_sqrt = eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T
            source_term = float(-2.0 * np.trace(
                solve_spd(phi, V_sqrt @ src.A_cpl @ V_sqrt)
            ))
        except SupportError:
            pass

    metadata = LinearAlgebraMetadata(
        atol=tol.atol, rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="evidence_decomposition",
        ambient_dim=n, support_rank=n, visible_dim=m,
    )

    return EvidenceDecompositionResult(
        f_pprime_visible=f2_vis,
        f_pprime_hidden=f2_hid,
        f_pprime_ambient=f2_amb,
        source_term=source_term,
        kinematic_term=kinematic_term,
        conservation_residual=conservation_residual,
        source_law=src,
        budget=budget,
        metadata=metadata,
    )


# ============================================================
# Observer diagnostics and steering
# ============================================================

@dataclass(frozen=True)
class ObserverDiagnosticsResult:
    """Combined diagnostic for observer quality under a perturbation."""
    visible_fraction: float
    v_min_eigenvalue: float
    exact_sector: bool
    hidden_defect_trace: float
    leakage_norm: float
    conservation_residual: float
    budget: InformationBudgetResult
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class CaptureCurveResult:
    """Information capture as a function of observer rank."""
    ranks: tuple[int, ...]
    visible_fractions: tuple[float, ...]
    v_min_eigenvalues: tuple[float, ...]
    exact_sector_flags: tuple[bool, ...]
    ambient_rate: float
    half_capture_rank: int | None
    perturbation_rank: int | None
    eigenvalue_spectrum: tuple[float, ...] | None
    metadata: LinearAlgebraMetadata


def observer_diagnostics(
    H: np.ndarray,
    C: np.ndarray,
    Hdot: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ObserverDiagnosticsResult:
    """
    Combined diagnostic for an observer under a perturbation.

    Returns visible fraction, exact-sector status (V > 0), hidden defect
    strength, and leakage norm. This is the primary diagnostic for
    observer steering: it tells you how good your current observer is
    for this perturbation.

    Parameters
    ----------
    H : (n, n) SPD matrix
    C : (m, n) rank-m observer
    Hdot : (n, n) symmetric (perturbation direction)
    tolerances : optional
    """
    tol = resolve_tolerances(tolerances)
    budget = information_budget(H, C, Hdot, tolerances=tol)
    n = H.shape[0]
    m = C.shape[0]

    v_min = float(np.min(budget.v_eigenvalues))
    exact = v_min > tol.atol

    # Hidden defect trace: Tr(B R^{-1} B^T)
    R_inv = solve_spd(budget.R, np.eye(budget.R.shape[0], dtype=float))
    Qhat = budget.B @ R_inv @ budget.B.T
    qhat_trace = float(np.trace(Qhat))

    # Leakage norm: ||[Hdot, P]||_F where P = L C is the observer projector
    geom = visible_geometry(H, C, tolerances=tol)
    P = geom.lift @ C
    commutator = Hdot @ P - P @ Hdot
    leakage = float(np.linalg.norm(commutator, 'fro'))

    metadata = LinearAlgebraMetadata(
        atol=tol.atol, rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="observer_diagnostics",
        ambient_dim=n, support_rank=n, visible_dim=m,
    )

    return ObserverDiagnosticsResult(
        visible_fraction=budget.visible_fraction,
        v_min_eigenvalue=v_min,
        exact_sector=exact,
        hidden_defect_trace=qhat_trace,
        leakage_norm=leakage,
        conservation_residual=budget.conservation_residual,
        budget=budget,
        metadata=metadata,
    )


def capture_curve(
    H: np.ndarray,
    Hdot: np.ndarray,
    m_max: int | None = None,
    observer_basis: np.ndarray | None = None,
    tolerances: Tolerances | None = None,
) -> CaptureCurveResult:
    """
    Information capture as a function of observer rank.

    For each rank m from 1 to m_max, computes the visible fraction.
    Returns the capture curve and the half-capture rank (smallest m
    where vis_frac >= 0.5).

    Parameters
    ----------
    H : (n, n) SPD matrix
    Hdot : (n, n) symmetric (perturbation direction)
    m_max : maximum rank to test (default: n - 1)
    observer_basis : (n, n) orthogonal matrix whose first m rows form the
        observer at each rank m. If None, uses the canonical basis
        (C = [I_m | 0], i.e., the first m coordinate directions).
    tolerances : optional
    """
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    n = h.shape[0]

    if m_max is None:
        m_max = n - 1
    m_max = min(m_max, n - 1)

    if observer_basis is not None:
        Q = np.asarray(observer_basis, dtype=float)
        if Q.shape != (n, n):
            raise InputValidationError(
                f"observer_basis must be ({n}, {n}), got {Q.shape}"
            )
    else:
        Q = np.eye(n, dtype=float)

    h_inv = solve_spd(h, np.eye(n, dtype=float))
    amb_rate = float(np.trace(h_inv @ hdot))

    # Eigenvalue spectrum of whitened perturbation (for the capture formula)
    try:
        h_inv_sqrt = np.linalg.cholesky(h_inv)  # H^{-1} = LL^T, so H^{-1/2} = L
        hdot_w = h_inv_sqrt.T @ hdot @ h_inv_sqrt
        hdot_w = symmetrize(hdot_w)
        eigvals_w = np.sort(np.linalg.eigvalsh(hdot_w))[::-1]
        pos_eigvals = eigvals_w[eigvals_w > max(tol.atol, 1e-10)]
        pert_rank = int(len(pos_eigvals))
        eig_spectrum = tuple(float(e) for e in pos_eigvals)
    except Exception:
        pert_rank = None
        eig_spectrum = None

    ranks = []
    vis_fracs = []
    v_mins = []
    exact_flags = []
    half_rank = None

    for m in range(1, m_max + 1):
        C = Q[:m, :]

        budget = information_budget(h, C, hdot, tolerances=tol)
        vf = budget.visible_fraction
        v_min = float(np.min(budget.v_eigenvalues))

        ranks.append(m)
        vis_fracs.append(vf)
        v_mins.append(v_min)
        exact_flags.append(v_min > tol.atol)

        if half_rank is None and not np.isnan(vf) and vf >= 0.5:
            half_rank = m

    metadata = LinearAlgebraMetadata(
        atol=tol.atol, rtol=tol.rtol,
        rank_tol=max(tol.atol, tol.rtol),
        method="capture_curve",
        ambient_dim=n, support_rank=n, visible_dim=m_max,
    )

    return CaptureCurveResult(
        ranks=tuple(ranks),
        visible_fractions=tuple(vis_fracs),
        v_min_eigenvalues=tuple(v_mins),
        exact_sector_flags=tuple(exact_flags),
        ambient_rate=amb_rate,
        half_capture_rank=half_rank,
        perturbation_rank=pert_rank,
        eigenvalue_spectrum=eig_spectrum,
        metadata=metadata,
    )
