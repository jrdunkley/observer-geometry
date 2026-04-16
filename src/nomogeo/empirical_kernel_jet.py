"""
Sparse empirical kn-jet extraction — 0.3.3 diagnostic layer.

This module extracts the polynomial jet of a black-box action function
directly in kernel/normal coordinates, using only the sparse subset of
monomials needed by the kernel reducer (Layer 4).

For kernel_dim r and normal_dim p, the full ambient jet to order 4 has
O((r+p)^4) terms.  But the kernel reducer only needs:

  - Pure kernel terms to degree 4:  O(r^4) terms
  - Cubic coupling V_j(u) = coefficient of y_j at degree 3: O(r^2 * p) terms

For the common case r=1, p=O(d), this is O(d) evaluations, not O(d^4).
For r=2 it is O(d), etc.

IMPORTANT: This is a diagnostic/numerical tool, not a certified
mathematical input.  The jets it produces are suitable for empirical
exploration and benchmarking, but NOT for theorem-backed claims.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Mapping

from .kernel_reduction import PolynomialJet

Array = NDArray[np.float64]


def extract_sparse_kn_jet(
    action: callable,
    theta_star: Array,
    kernel_basis: Array,
    normal_basis: Array,
    *,
    order: int = 4,
    h: float = 1e-4,
    metadata: Mapping[str, Any] | None = None,
) -> PolynomialJet:
    """Extract a sparse jet in kernel/normal coordinates.

    Computes only the terms needed by reduce_kernel_action_kn:
      - All pure-kernel monomials up to `order` (degree 2, 3, 4)
      - All cubic coupling terms u^alpha * y_j with |alpha|=2, |e_j|=1
      - The normal Hessian diagonal (y_j^2 terms, degree 2)
      - Cross-normal Hessian terms (y_j * y_k, degree 2)

    Parameters
    ----------
    action : callable
        f(theta) → scalar.  The local action (e.g., -ℓ(θ)/n).
    theta_star : array, shape (D,)
        The stationary point in ambient coordinates.
    kernel_basis : array, shape (D, r)
        Orthonormal kernel eigenvectors.
    normal_basis : array, shape (D, p)
        Orthonormal positive-normal eigenvectors.
    order : int
        Maximum total degree for pure-kernel terms (default 4).
    h : float
        Step size for finite differences (default 1e-4).

    Returns
    -------
    PolynomialJet
        Jet in kn-coordinates (dim = r + p), with only the sparse
        terms needed for kernel reduction.

        IMPORTANT: certified_order is set to 2 (not `order`), because
        the normal sector beyond degree 2 is unprobed.  The pure kernel
        terms and cubic coupling terms are present and may be accurate,
        but the jet cannot claim completeness at higher degrees.  The
        metadata carries explicit provenance marking this as a
        finite-difference diagnostic extraction, not a certified input.
    """
    theta_star = np.asarray(theta_star, dtype=float)
    kernel_basis = np.asarray(kernel_basis, dtype=float)
    normal_basis = np.asarray(normal_basis, dtype=float)

    r = kernel_basis.shape[1]
    p = normal_basis.shape[1]
    d = r + p

    f0 = action(theta_star)

    # Helper: evaluate f at theta_star + sum of displacements in kn-coords.
    def f_at(displacements: dict[int, float]) -> float:
        """Evaluate action at theta_star + sum_i displacements[i] * basis_i.

        Keys 0..r-1 index kernel directions, r..r+p-1 index normal directions.
        """
        delta = np.zeros_like(theta_star)
        for idx, scale in displacements.items():
            if idx < r:
                delta += scale * kernel_basis[:, idx]
            else:
                delta += scale * normal_basis[:, idx - r]
        return action(theta_star + delta)

    terms: dict[tuple[int, ...], float] = {}

    # ══════════════════════════════════════════════════════════════════
    # 1. Pure kernel terms: all monomials in u_1,...,u_r up to `order`
    # ══════════════════════════════════════════════════════════════════
    pure_kernel_indices = _multi_indices(r, max_deg=order, min_deg=2)
    for alpha in pure_kernel_indices:
        coeff = _fd_derivative(f_at, alpha, r, p, h, f0) / _factorial_multi(alpha)
        if abs(coeff) > 1e-15:
            # Pad to full kn dimension
            full_idx = alpha + tuple(0 for _ in range(p))
            terms[full_idx] = coeff

    # ══════════════════════════════════════════════════════════════════
    # 2. Normal Hessian terms: y_j^2 and y_j*y_k (degree 2)
    # ══════════════════════════════════════════════════════════════════
    for j in range(p):
        # y_j^2 term: second derivative in normal direction j
        kn_j = r + j
        fp = f_at({kn_j: h})
        fm = f_at({kn_j: -h})
        d2 = (fp - 2 * f0 + fm) / (h * h)
        coeff = d2 / 2.0  # coefficient of y_j^2 = (1/2!) * d^2f/dy_j^2
        if abs(coeff) > 1e-15:
            idx = tuple(0 for _ in range(r)) + tuple(2 if i == j else 0 for i in range(p))
            terms[idx] = coeff

        # y_j * y_k cross terms (j < k)
        for k in range(j + 1, p):
            kn_k = r + k
            fpp = f_at({kn_j: h, kn_k: h})
            fpm = f_at({kn_j: h, kn_k: -h})
            fmp = f_at({kn_j: -h, kn_k: h})
            fmm = f_at({kn_j: -h, kn_k: -h})
            d2_jk = (fpp - fpm - fmp + fmm) / (4 * h * h)
            # coefficient of y_j * y_k = d^2f / (dy_j dy_k)
            if abs(d2_jk) > 1e-15:
                idx = tuple(0 for _ in range(r)) + tuple(
                    1 if i == j or i == k else 0 for i in range(p)
                )
                terms[idx] = d2_jk

    # ══════════════════════════════════════════════════════════════════
    # 3. Cubic coupling: u^alpha * y_j with |alpha|=2
    # ══════════════════════════════════════════════════════════════════
    # These are ∂^3 Ψ / (∂u^alpha ∂y_j) / alpha! at the origin.
    # For each u-multi-index alpha with |alpha|=2 and each normal j,
    # compute the mixed derivative.
    coupling_u_indices = _multi_indices(r, max_deg=2, min_deg=2)
    for alpha in coupling_u_indices:
        for j in range(p):
            # The mixed derivative ∂^{|alpha|+1} f / (∂u^alpha ∂y_j)
            # We use the fact that for a single y_j derivative:
            #   ∂/∂y_j [g(u)] = [g(u, y_j=h) - g(u, y_j=-h)] / (2h)
            # and g(u) is already differentiated in u.
            coeff = _fd_mixed_coupling(f_at, alpha, j, r, p, h, f0) / _factorial_multi(alpha)
            if abs(coeff) > 1e-15:
                full_idx = alpha + tuple(1 if i == j else 0 for i in range(p))
                terms[full_idx] = coeff

    # ══════════════════════════════════════════════════════════════════
    # Certification bookkeeping
    # ══════════════════════════════════════════════════════════════════
    # The sparse extraction probes:
    #   - Pure kernel terms up to `order` (complete in the kernel sector)
    #   - Normal Hessian (degree 2, complete)
    #   - Cubic coupling u^alpha * y_j with |alpha|=2 (complete)
    #
    # It does NOT probe:
    #   - Higher normal nonlinearities (u*y^2, y^3, u^2*y^2, etc.)
    #   - These can feed into the quartic correction through η₂
    #
    # Therefore the jet is NOT certified complete at degree 3 or 4 in
    # the mixed (kernel × normal) sector.  The reducer's quartic
    # correction may be incomplete if unprobed terms are materially
    # nonzero.
    #
    # We set certified_order = 2 to honestly reflect that only the
    # quadratic normal Hessian is certified complete.  The pure kernel
    # terms and cubic coupling terms ARE present (and may be accurate),
    # but the jet cannot claim completeness at higher degrees because
    # the normal sector is only partially probed.
    #
    # The reducer will see certified_order=2, detect that cubic coupling
    # IS present, and set output certified_order=2 (the conservative
    # branch for "has_cubic_coupling and not certified >= 4").
    # This correctly signals that the reduced germ is diagnostic, not
    # theorem-backed.

    meta = {
        "extraction": "sparse_empirical_kn_jet",
        "provenance": "finite_difference",
        "certified": False,
        "h": h,
        "kernel_dim": r,
        "normal_dim": p,
        "probed_families": [
            "pure_kernel_to_degree_" + str(order),
            "normal_hessian_degree_2",
            "cubic_coupling_u2_yj",
        ],
        "unprobed_families": [
            "higher_normal_nonlinearities (u*y^2, y^3, u^2*y^2, ...)",
        ],
        "note": (
            "Sparse empirical jet: only the term families listed in "
            "'probed_families' were extracted.  The normal sector beyond "
            "degree 2 is unprobed.  The reducer's quartic correction may "
            "be incomplete if unprobed higher-normal terms are materially "
            "nonzero.  Do not use for theorem-backed claims."
        ),
    }
    if metadata:
        meta.update(metadata)

    return PolynomialJet(
        dim=d,
        terms=terms,
        certified_order=2,
        metadata=meta,
    )


# ═══════════════════════════════════════════════════════════════════════
# Finite difference helpers
# ═══════════════════════════════════════════════════════════════════════

def _multi_indices(dim: int, max_deg: int, min_deg: int = 0) -> list[tuple[int, ...]]:
    """Generate all multi-indices in `dim` variables with min_deg ≤ |α| ≤ max_deg."""
    result = []
    for total in range(min_deg, max_deg + 1):
        result.extend(_multi_indices_exact(dim, total))
    return result


def _multi_indices_exact(dim: int, total_deg: int) -> list[tuple[int, ...]]:
    """Generate all multi-indices of exactly total_deg in dim variables."""
    if dim == 0:
        return [()] if total_deg == 0 else []
    if dim == 1:
        return [(total_deg,)]
    result = []
    for first in range(total_deg + 1):
        for rest in _multi_indices_exact(dim - 1, total_deg - first):
            result.append((first,) + rest)
    return result


def _factorial_multi(alpha: tuple[int, ...]) -> float:
    """Compute α! = α₁! · α₂! · ..."""
    result = 1.0
    for a in alpha:
        for i in range(2, a + 1):
            result *= i
    return result


def _fd_derivative(
    f_at: callable,
    alpha: tuple[int, ...],
    r: int,
    p: int,
    h: float,
    f0: float,
) -> float:
    """Compute ∂^|α| f / ∂u^α at the origin via iterated central differences.

    alpha is a multi-index over the first r (kernel) coordinates.
    """
    # Build stencil via iterated central differences
    stencil: dict[tuple, float] = {(): 1.0}

    for i, a_i in enumerate(alpha):
        if a_i == 0:
            continue
        for _ in range(a_i):
            new_stencil: dict[tuple, float] = {}
            for offsets, weight in stencil.items():
                off_dict_plus = dict(offsets) if isinstance(offsets, tuple) and len(offsets) > 0 and isinstance(offsets[0], tuple) else {}
                # Represent offsets as dict of {kn_index: displacement}
                # Actually, let's use a simpler stencil representation
                pass
            # Simpler: use the offset-vector approach
            break
        break

    # Use a cleaner stencil approach: represent as {frozenset of (idx, displacement): weight}
    # Actually, let's just use the recursive approach from the full pipeline but restricted
    # to kernel directions only.
    return _central_diff_pure_kernel(f_at, alpha, r, p, h)


def _central_diff_pure_kernel(
    f_at: callable,
    alpha: tuple[int, ...],
    r: int,
    p: int,
    h: float,
) -> float:
    """Central differences for pure-kernel derivatives.

    Stencil stored as dict mapping {kn_index: displacement} → weight.
    """
    # Start with unit stencil at origin
    stencil: dict[tuple[tuple[int, float], ...], float] = {(): 1.0}

    for i, a_i in enumerate(alpha):
        if a_i == 0:
            continue
        for _ in range(a_i):
            new_stencil: dict[tuple[tuple[int, float], ...], float] = {}
            for offsets, weight in stencil.items():
                off_dict = dict(offsets)
                # +h in direction i
                plus_dict = dict(off_dict)
                plus_dict[i] = plus_dict.get(i, 0.0) + h
                plus_key = tuple(sorted(plus_dict.items()))
                new_stencil[plus_key] = new_stencil.get(plus_key, 0.0) + weight / (2 * h)
                # -h in direction i
                minus_dict = dict(off_dict)
                minus_dict[i] = minus_dict.get(i, 0.0) - h
                minus_key = tuple(sorted(minus_dict.items()))
                new_stencil[minus_key] = new_stencil.get(minus_key, 0.0) - weight / (2 * h)
            stencil = new_stencil

    # Evaluate
    result = 0.0
    for offsets, weight in stencil.items():
        if abs(weight) < 1e-30:
            continue
        off_dict = dict(offsets)
        # Clean up near-zero offsets
        off_dict = {k: v for k, v in off_dict.items() if abs(v) > 1e-15}
        result += weight * f_at(off_dict)

    return result


def _fd_mixed_coupling(
    f_at: callable,
    alpha: tuple[int, ...],
    j: int,
    r: int,
    p: int,
    h: float,
    f0: float,
) -> float:
    """Compute ∂^{|α|+1} f / (∂u^α ∂y_j) at the origin.

    Uses the identity:
        ∂/∂y_j [g(u)] = [g(u, y_j=+h) - g(u, y_j=-h)] / (2h)

    where g(u) = ∂^|α| f / ∂u^α evaluated at y=0 with y_j shifted.
    """
    kn_j = r + j  # index of y_j in the kn coordinate system

    # Build the kernel stencil for ∂^|α| / ∂u^α
    stencil: dict[tuple[tuple[int, float], ...], float] = {(): 1.0}

    for i, a_i in enumerate(alpha):
        if a_i == 0:
            continue
        for _ in range(a_i):
            new_stencil: dict[tuple[tuple[int, float], ...], float] = {}
            for offsets, weight in stencil.items():
                off_dict = dict(offsets)
                plus_dict = dict(off_dict)
                plus_dict[i] = plus_dict.get(i, 0.0) + h
                plus_key = tuple(sorted(plus_dict.items()))
                new_stencil[plus_key] = new_stencil.get(plus_key, 0.0) + weight / (2 * h)
                minus_dict = dict(off_dict)
                minus_dict[i] = minus_dict.get(i, 0.0) - h
                minus_key = tuple(sorted(minus_dict.items()))
                new_stencil[minus_key] = new_stencil.get(minus_key, 0.0) - weight / (2 * h)
            stencil = new_stencil

    # Now apply one central difference in the y_j direction to each stencil point
    result = 0.0
    for offsets, weight in stencil.items():
        if abs(weight) < 1e-30:
            continue
        off_dict = dict(offsets)
        # Evaluate at +h in y_j
        plus_dict = dict(off_dict)
        plus_dict[kn_j] = plus_dict.get(kn_j, 0.0) + h
        f_plus = f_at(plus_dict)
        # Evaluate at -h in y_j
        minus_dict = dict(off_dict)
        minus_dict[kn_j] = minus_dict.get(kn_j, 0.0) - h
        f_minus = f_at(minus_dict)
        result += weight * (f_plus - f_minus) / (2 * h)

    return result
