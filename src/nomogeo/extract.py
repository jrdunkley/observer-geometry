"""
Geometry extraction: turning data into (H, C, Hdot).

This module bridges the gap between raw data and the observer geometry kernel.
Every downstream function (information_budget, source_law, observer_diagnostics,
capture_curve) takes (H, C, Hdot). This module provides the extraction step
for common data formats.

Public API:
    extract_supervised     -- (X, y, task) → (H, Hdot, task_family)
    extract_covariance     -- (cov, perturbation) → (H, Hdot)
"""
from __future__ import annotations

import numpy as np

from .validation import (
    Tolerances,
    resolve_tolerances,
    symmetrize,
)
from .exceptions import InputValidationError

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SupervisedGeometry:
    """Extracted geometry from supervised (X, y) data.

    H is the regularised within-class precision.
    Hdot is the task-family aggregate (between-class scatter).
    class_scatters are the per-class scatter matrices.
    """
    H: np.ndarray
    Hdot: np.ndarray
    class_scatters: list[np.ndarray]
    class_labels: list
    n_features: int
    n_classes: int
    reg_floor: float


def extract_supervised(
    X: np.ndarray,
    y: np.ndarray,
    task: str = "fisher",
    reg_floor: float = 1e-4,
) -> SupervisedGeometry:
    """
    Extract observer geometry from supervised classification data.

    Parameters
    ----------
    X : (n_samples, n_features) array
        Feature matrix.
    y : (n_samples,) array
        Class labels.
    task : str
        Task family type. One of:
        - "fisher": sample-weighted between-class scatter (default)
        - "equal_weight": equal importance per class
        - "minority": inverse-frequency weighting
    reg_floor : float
        Regularisation floor for the within-class scatter. Default 1e-4.

    Returns
    -------
    SupervisedGeometry
        H (regularised within-class precision), Hdot (task aggregate),
        per-class scatters, labels, dimensions.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y)
    if X.ndim != 2:
        raise InputValidationError("X must be 2D")
    n_samples, n_features = X.shape
    if y.shape != (n_samples,):
        raise InputValidationError("y must have shape (n_samples,)")

    classes = np.unique(y)
    n_classes = len(classes)
    if n_classes < 2:
        raise InputValidationError("Need at least 2 classes")

    # Standardise
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std < 1e-12] = 1.0  # avoid division by zero
    X_s = (X - X_mean) / X_std

    # Within-class scatter
    S_w = np.zeros((n_features, n_features))
    for c in classes:
        Xc = X_s[y == c]
        Xc_centered = Xc - Xc.mean(axis=0)
        S_w += Xc_centered.T @ Xc_centered
    S_w /= n_samples
    S_w = symmetrize(S_w)

    # Regularise and invert
    H = np.linalg.inv(S_w + reg_floor * np.eye(n_features))
    H = symmetrize(H)

    # Between-class scatters
    grand_mean = X_s.mean(axis=0)
    class_scatters = []
    class_counts = []

    for c in classes:
        Xc = X_s[y == c]
        diff = Xc.mean(axis=0) - grand_mean
        class_counts.append(len(Xc))
        class_scatters.append(np.outer(diff, diff))

    # Task weighting
    if task == "fisher":
        weights = np.array([n / n_samples for n in class_counts])
    elif task == "equal_weight":
        weights = np.ones(n_classes) / n_classes
    elif task == "minority":
        inv_freq = np.array([n_samples / max(n, 1) for n in class_counts])
        weights = inv_freq / inv_freq.sum()
    else:
        raise InputValidationError(f"Unknown task: {task!r}. Use 'fisher', 'equal_weight', or 'minority'.")

    weighted_scatters = [w * s for w, s in zip(weights, class_scatters)]
    Hdot = symmetrize(sum(weighted_scatters))

    return SupervisedGeometry(
        H=H,
        Hdot=Hdot,
        class_scatters=weighted_scatters,
        class_labels=classes.tolist(),
        n_features=n_features,
        n_classes=n_classes,
        reg_floor=reg_floor,
    )


def extract_covariance(
    covariance: np.ndarray,
    perturbation: np.ndarray,
    reg_floor: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract (H, Hdot) from a covariance matrix and its perturbation.

    Parameters
    ----------
    covariance : (n, n) symmetric PSD matrix
    perturbation : (n, n) symmetric matrix (the change in covariance)
    reg_floor : float
        Regularisation for inversion.

    Returns
    -------
    H : (n, n) SPD (regularised precision)
    Hdot : (n, n) symmetric (perturbation in precision space)
    """
    cov = np.asarray(covariance, dtype=float)
    pert = np.asarray(perturbation, dtype=float)
    n = cov.shape[0]

    cov = symmetrize(cov)
    pert = symmetrize(pert)

    H = np.linalg.inv(cov + reg_floor * np.eye(n))
    H = symmetrize(H)

    # The perturbation in precision space:
    # If Sigma -> Sigma + delta_Sigma, then H -> H - H delta_Sigma H (to first order)
    Hdot = symmetrize(-H @ pert @ H)

    return H, Hdot
