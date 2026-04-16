"""
Observer steering engine.

Given (H, C, Hdot), score the observer, find a better one, return diagnostics.
This is the control layer that nomocomp and nomoselect both call through
their domain-specific interfaces.

Public API:
    score_observer    -- compare an observer against canonical and PCA baselines
    optimize_observer -- find the best observer for a declared perturbation
    steer             -- top-level: extract + optimize + diagnose
"""
from __future__ import annotations

import numpy as np
from numpy.linalg import svd as _svd

from .core import visible_geometry
from .adapted import closure_adapted_observer
from .source import (
    information_budget,
    observer_diagnostics,
    capture_curve,
    InformationBudgetResult,
    ObserverDiagnosticsResult,
    CaptureCurveResult,
)
from .extract import extract_supervised, SupervisedGeometry
from .validation import (
    Tolerances,
    resolve_tolerances,
    solve_spd,
    symmetrize,
    validate_spd_matrix,
    validate_surjective_map,
    validate_symmetric_matrix,
)
from .exceptions import InputValidationError

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ObserverScoreResult:
    """Comparison of an observer against baselines.

    All fields use plain names. The underlying geometry is hidden.
    """
    # The observer being scored
    observer: np.ndarray
    rank: int

    # Primary diagnostics
    visible_rate: float
    hidden_rate: float
    ambient_rate: float
    visible_fraction: float
    exact_sector: bool
    leakage: float
    hidden_defect: float

    # Baseline comparisons
    canonical_visible_fraction: float
    pca_visible_fraction: float
    advantage_over_canonical: float
    advantage_over_pca: float

    # Conservation check
    conservation_residual: float


@dataclass(frozen=True)
class OptimizeResult:
    """Result of observer optimization."""
    observer: np.ndarray
    rank: int
    diagnostics: ObserverDiagnosticsResult
    score: ObserverScoreResult
    method: str  # "adapted", "adapted_fallback"
    capture: CaptureCurveResult


@dataclass(frozen=True)
class SteerResult:
    """Top-level steering result from raw data."""
    observer: np.ndarray
    rank: int
    visible_fraction: float
    exact_sector: bool
    pca_visible_fraction: float
    advantage_over_pca: float
    half_capture_rank: int | None
    conservation_residual: float
    optimization: OptimizeResult
    geometry: SupervisedGeometry | None


def score_observer(
    H: np.ndarray,
    C: np.ndarray,
    Hdot: np.ndarray,
    tolerances: Tolerances | None = None,
) -> ObserverScoreResult:
    """
    Score an observer against canonical and PCA baselines.

    Parameters
    ----------
    H : (n, n) SPD matrix
    C : (m, n) observer
    Hdot : (n, n) symmetric perturbation
    tolerances : optional

    Returns
    -------
    ObserverScoreResult with diagnostics and baseline comparisons
    """
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    n = h.shape[0]
    m = C.shape[0]

    # Score the given observer
    diag = observer_diagnostics(h, C, hdot, tolerances=tol)

    # Canonical baseline: C = [I_m | 0]
    C_canon = np.zeros((m, n), dtype=float)
    for i in range(m):
        C_canon[i, i] = 1.0
    b_canon = information_budget(h, C_canon, hdot, tolerances=tol)

    # PCA baseline: top m eigenvectors of H
    eigvals_H, eigvecs_H = np.linalg.eigh(h)
    C_pca = eigvecs_H[:, -m:].T
    b_pca = information_budget(h, C_pca, hdot, tolerances=tol)

    return ObserverScoreResult(
        observer=C,
        rank=m,
        visible_rate=diag.budget.visible_rate,
        hidden_rate=diag.budget.hidden_rate,
        ambient_rate=diag.budget.ambient_rate,
        visible_fraction=diag.visible_fraction,
        exact_sector=diag.exact_sector,
        leakage=diag.leakage_norm,
        hidden_defect=diag.hidden_defect_trace,
        canonical_visible_fraction=b_canon.visible_fraction,
        pca_visible_fraction=b_pca.visible_fraction,
        advantage_over_canonical=diag.visible_fraction - b_canon.visible_fraction,
        advantage_over_pca=diag.visible_fraction - b_pca.visible_fraction,
        conservation_residual=diag.conservation_residual,
    )


def optimize_observer(
    H: np.ndarray,
    Hdot: np.ndarray,
    rank: int,
    tolerances: Tolerances | None = None,
) -> OptimizeResult:
    """
    Find the best observer for a declared perturbation.

    Uses the adapted observer (zero hidden coupling, milliseconds).
    Falls back gracefully if the adapted observer fails.

    Parameters
    ----------
    H : (n, n) SPD matrix
    Hdot : (n, n) symmetric perturbation
    rank : int, the observer rank m
    tolerances : optional

    Returns
    -------
    OptimizeResult with observer, diagnostics, score, and capture curve
    """
    tol = resolve_tolerances(tolerances)
    h = validate_spd_matrix("H", H, tol)
    hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    n = h.shape[0]

    if rank < 1 or rank >= n:
        raise InputValidationError(f"rank must be in [1, {n-1}], got {rank}")

    # Tier 1: adapted observer (B = 0, zero hidden coupling)
    method = "adapted"
    try:
        result = closure_adapted_observer(h, [hdot], rank, tolerances=tol)
        C_opt = result.C
    except Exception:
        # Fallback: use PCA (top eigenvectors of H)
        method = "adapted_fallback"
        eigvals_H, eigvecs_H = np.linalg.eigh(h)
        C_opt = eigvecs_H[:, -rank:].T

    diag = observer_diagnostics(h, C_opt, hdot, tolerances=tol)
    sc = score_observer(h, C_opt, hdot, tolerances=tol)
    cc = capture_curve(h, hdot, m_max=min(rank + 5, n - 1), tolerances=tol)

    return OptimizeResult(
        observer=C_opt,
        rank=rank,
        diagnostics=diag,
        score=sc,
        method=method,
        capture=cc,
    )


def steer(
    X: np.ndarray | None = None,
    y: np.ndarray | None = None,
    H: np.ndarray | None = None,
    Hdot: np.ndarray | None = None,
    rank: int = 3,
    task: str = "fisher",
    tolerances: Tolerances | None = None,
) -> SteerResult:
    """
    Top-level steering: extract geometry, optimize observer, return diagnostics.

    Can be called two ways:
    1. With (X, y): extracts supervised geometry automatically
    2. With (H, Hdot): uses raw matrices directly

    Parameters
    ----------
    X : (n_samples, n_features) array, optional
    y : (n_samples,) labels, optional
    H : (n, n) SPD matrix, optional
    Hdot : (n, n) symmetric, optional
    rank : observer rank (default 3)
    task : task family for supervised extraction (default "fisher")
    tolerances : optional

    Returns
    -------
    SteerResult with observer, diagnostics, and comparisons
    """
    tol = resolve_tolerances(tolerances)
    geom = None

    if X is not None and y is not None:
        geom = extract_supervised(X, y, task=task)
        h = geom.H
        hdot = geom.Hdot
    elif H is not None and Hdot is not None:
        h = validate_spd_matrix("H", H, tol)
        hdot = validate_symmetric_matrix("Hdot", Hdot, tol)
    else:
        raise InputValidationError(
            "Provide either (X, y) for supervised data or (H, Hdot) for raw matrices."
        )

    n = h.shape[0]
    rank = min(rank, n - 1)

    opt = optimize_observer(h, hdot, rank, tolerances=tol)

    return SteerResult(
        observer=opt.observer,
        rank=rank,
        visible_fraction=opt.score.visible_fraction,
        exact_sector=opt.diagnostics.exact_sector,
        pca_visible_fraction=opt.score.pca_visible_fraction,
        advantage_over_pca=opt.score.advantage_over_pca,
        half_capture_rank=opt.capture.half_capture_rank,
        conservation_residual=opt.score.conservation_residual,
        optimization=opt,
        geometry=geom,
    )
