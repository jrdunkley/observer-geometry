"""
Kernel reduction — 0.3.3 Technical Note §9 (Proposition 9.1).

This module implements Layer 4 of the 5-layer stack: given a
QuadraticKernelSeed from the regime detector and a certified local
action representation, reduce the singular local evidence problem to a
lower-dimensional kernel integral by integrating out the positive normal
directions exactly at the Laplace level.

The key mathematical content is Proposition 9.1:

    ∫_{U∩K} e^{-nΨ(x)} a(x) dx
    = n^{-dim N/2} ∫_{U_R∩K_R} e^{-nΦ(u)} (b(u) + o(1)) du

where:
    - Ψ is the local visible action
    - E_K = R ⊕ N is the kernel ⊕ positive-normal decomposition
    - η: R → N is the critical graph (∂_y Ψ(u, η(u)) = 0)
    - Φ(u) = Ψ(u, η(u)) is the reduced kernel action
    - b(0) = a(0) (2π)^{dim N/2} det(H_N)^{-1/2}

The module supports theorem-backed inputs only:
    - Explicit polynomial jet representations (PolynomialJet)
    - The reduced kernel action as a polynomial on the kernel space

It does NOT support:
    - Finite-difference jet extraction (that is diagnostic, not certified)
    - Arbitrary black-box callables as certified inputs

References
----------
0.3.3 Technical Note:
  - §9  Proposition 9.1 (kernel reduction)
  - §10 Remark (what the reduction achieves)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from .exceptions import InputValidationError
from .regime_types import (
    ConeKind,
    ConeSpec,
    QuadraticKernelSeed,
)

Array = NDArray[np.float64]


# ═══════════════════════════════════════════════════════════════════════
# Polynomial jet representation
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class PolynomialJet:
    """A polynomial jet representation of a local action near a critical point.

    The action Ψ is represented as an explicit polynomial in local coordinates
    centred at the optimum, with Ψ(0) = 0 and ∇Ψ(0) = 0 (critical point).

    The polynomial is stored as a dictionary mapping multi-index tuples to
    coefficients.  For example, for Ψ(x₁,x₂) = 3x₁² + 2x₁x₂ + x₂⁴:

        terms = {(2,0): 3.0, (1,1): 2.0, (0,4): 1.0}

    Only terms of total degree ≥ 2 are meaningful (degree 0 is the action
    value at the critical point, already factored out; degree 1 is zero by
    the critical point condition).

    This representation is theorem-backed: the user certifies that the
    polynomial accurately represents the action to the required order.
    """

    dim: int
    """Number of variables."""

    terms: dict[tuple[int, ...], float]
    """Multi-index → coefficient mapping.  Each key is a tuple of length
    dim giving the exponents.  Example: (2, 0, 1) means x₁²x₃."""

    certified_order: int = 4
    """The total degree to which this jet is certified accurate.
    Higher-order terms may be present but are not guaranteed.
    For kernel reduction, order ≥ 3 is typically needed."""

    metadata: Mapping[str, Any] | None = None

    def evaluate(self, x: Array) -> float:
        """Evaluate the polynomial at a point x."""
        x = np.asarray(x, dtype=float)
        if x.shape != (self.dim,):
            raise InputValidationError(
                f"Expected x of shape ({self.dim},), got {x.shape}"
            )
        result = 0.0
        for multi_idx, coeff in self.terms.items():
            monomial = 1.0
            for i, exp in enumerate(multi_idx):
                if exp > 0:
                    monomial *= x[i] ** exp
            result += coeff * monomial
        return result

    def gradient(self, x: Array) -> Array:
        """Evaluate the gradient at a point x."""
        x = np.asarray(x, dtype=float)
        grad = np.zeros(self.dim, dtype=float)
        for multi_idx, coeff in self.terms.items():
            for j in range(self.dim):
                if multi_idx[j] == 0:
                    continue
                # Derivative with respect to x_j.
                new_exp = list(multi_idx)
                new_exp[j] -= 1
                d_coeff = coeff * multi_idx[j]
                monomial = 1.0
                for i, exp in enumerate(new_exp):
                    if exp > 0:
                        monomial *= x[i] ** exp
                grad[j] += d_coeff * monomial
        return grad

    def hessian_at_origin(self) -> Array:
        """Extract the Hessian matrix ∂²Ψ/∂xᵢ∂xⱼ at the origin.

        For a polynomial with terms stored as multi-index→coefficient,
        the Hessian at the origin picks out exactly the degree-2 terms.
        """
        H = np.zeros((self.dim, self.dim), dtype=float)
        for multi_idx, coeff in self.terms.items():
            total_deg = sum(multi_idx)
            if total_deg != 2:
                continue
            # Find which indices are involved.
            nonzero = [(i, e) for i, e in enumerate(multi_idx) if e > 0]
            if len(nonzero) == 1:
                # Pure quadratic: c * x_i^2 → H[i,i] = 2c
                i, e = nonzero[0]
                H[i, i] += 2.0 * coeff
            elif len(nonzero) == 2:
                # Cross term: c * x_i * x_j → H[i,j] = H[j,i] = c
                i, _ = nonzero[0]
                j, _ = nonzero[1]
                H[i, j] += coeff
                H[j, i] += coeff
        return H

    def restrict_to_subspace(self, basis: Array) -> "PolynomialJet":
        """Restrict the polynomial to a subspace given by columns of basis.

        If basis has shape (dim, k), the restricted polynomial acts on
        k-dimensional coordinates u, with x = basis @ u.

        This is exact for polynomial jets.
        """
        basis = np.asarray(basis, dtype=float)
        if basis.shape[0] != self.dim:
            raise InputValidationError(
                f"basis has {basis.shape[0]} rows, expected {self.dim}"
            )
        k = basis.shape[1]

        # For each term c * x^α, substitute x = B u and collect monomials.
        # This is exact but can be expensive for high-order terms.
        new_terms: dict[tuple[int, ...], float] = {}

        for multi_idx, coeff in self.terms.items():
            # Expand (B₁·u)^α₁ · (B₂·u)^α₂ · ... as a polynomial in u.
            # Use recursive multinomial expansion.
            _expand_term(multi_idx, coeff, basis, k, new_terms)

        # Clean up near-zero terms.
        cleaned = {k: v for k, v in new_terms.items() if abs(v) > 1e-15}
        return PolynomialJet(dim=k, terms=cleaned, certified_order=self.certified_order)


# ═══════════════════════════════════════════════════════════════════════
# True reduced kernel datum (paper's post-reduction object)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ReducedKernelActionDatum:
    """The true reduced kernel datum from Proposition 9.1.

    This is the paper's "reduced kernel datum" — the post-reduction singular
    object that carries:
    - the reduced kernel action germ Φ (as a polynomial jet on R)
    - the residual active cone K_R on the kernel space
    - the normal determinant prefactor from integrating N
    - kernel chart metadata

    This is the input to the certified singular classifier (Layer 5).
    It is NOT the QuadraticKernelSeed (which is only the quadratic-order
    decomposition from the detector).
    """

    kernel_dim: int
    """Dimension of the kernel space R."""

    kernel_action: PolynomialJet
    """The reduced kernel action germ Φ(u) = Ψ(u, η(u)) on R.
    This is a polynomial in kernel coordinates, with Φ(0) = 0 and
    ∇Φ(0) = 0 (since the original critical point restricts)."""

    kernel_cone: ConeSpec
    """The residual admissible cone K_R on the kernel space R."""

    positive_normal_dim: int
    """Dimension of the positive normal complement N."""

    log_normal_prefactor: float
    """log of the normal-direction Laplace prefactor:
    (dim_N / 2) * log(2π) - (1/2) * sum(log(eigenvalues)).
    This is the exact factor from integrating out N."""

    total_dim: int = 0
    """Original active dimension (= kernel_dim + positive_normal_dim)."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Primary API: reduce_kernel_action
# ═══════════════════════════════════════════════════════════════════════

def reduce_kernel_action(
    seed: QuadraticKernelSeed,
    action_jet: PolynomialJet,
    *,
    metadata: Mapping[str, Any] | None = None,
) -> ReducedKernelActionDatum:
    """Perform exact kernel reduction per Proposition 9.1.

    Given the QuadraticKernelSeed (kernel ⊕ normal decomposition from the
    detector) and a certified polynomial jet for the local action Ψ,
    compute the reduced kernel action Φ(u) = Ψ(u, η(u)) where η is the
    critical graph obtained from the implicit function theorem.

    IMPORTANT: this entry point expects the action jet in the ACTIVE
    coordinate system (the ambient eigenbasis).  It internally rotates
    the jet into kernel/normal coordinates via restrict_to_subspace
    before performing the reduction.

    If you already have a jet in kernel/normal coordinates (first r
    coordinates = kernel, next p = normal), use reduce_kernel_action_kn
    instead.  Passing a kn-coordinate jet to this function would
    double-rotate and produce wrong results.

    Parameters
    ----------
    seed : QuadraticKernelSeed
        Output of the quadratic regime detector for an unresolved-kernel case.
    action_jet : PolynomialJet
        Certified polynomial jet for the local action Ψ in the active
        coordinate system.  Must have Ψ(0) = 0, ∇Ψ(0) = 0, and its
        Hessian must be consistent with the seed's eigendecomposition.

    Returns
    -------
    ReducedKernelActionDatum
        The paper's reduced kernel datum: germ Φ on R, residual cone K_R,
        and normal prefactor.
    """
    # Validate dimensions.
    d = seed.kernel_dim + seed.positive_normal_dim
    if action_jet.dim != d:
        raise InputValidationError(
            f"action_jet.dim={action_jet.dim} does not match "
            f"kernel_dim + normal_dim = {d}"
        )

    kernel_basis = seed.kernel_basis       # shape (d, r)
    normal_basis = seed.positive_normal_basis  # shape (d, p)

    # Change to kernel/normal coordinates: x = [kernel_basis | normal_basis] @ z
    # where z = (u, y) with u ∈ R^r (kernel) and y ∈ R^p (normal).
    full_basis = np.hstack([kernel_basis, normal_basis])  # shape (d, d)

    # Express Ψ in the (u, y) coordinate system.
    action_in_kn = action_jet.restrict_to_subspace(full_basis)

    # Delegate to shared core.
    return _reduce_kernel_action_from_kn(seed, action_in_kn, metadata=metadata)


def reduce_kernel_action_kn(
    seed: QuadraticKernelSeed,
    action_kn_jet: PolynomialJet,
    *,
    metadata: Mapping[str, Any] | None = None,
) -> ReducedKernelActionDatum:
    """Perform kernel reduction on a jet already in kernel/normal coordinates.

    This is the sister entry point to reduce_kernel_action.  Use this
    when the action jet was extracted directly in the kernel/normal
    coordinate system (first r coordinates = kernel directions, next p
    coordinates = positive-normal directions).

    This avoids the O(p^order) restrict_to_subspace computation and
    prevents the double-rotation bug that occurs if a kn-coordinate jet
    is passed to reduce_kernel_action.

    Parameters
    ----------
    seed : QuadraticKernelSeed
        Output of the quadratic regime detector for an unresolved-kernel case.
    action_kn_jet : PolynomialJet
        Certified polynomial jet for the local action Ψ already expressed
        in kernel/normal coordinates: z = (u₁,...,uᵣ, y₁,...,yₚ).
        Must have Ψ(0) = 0, ∇Ψ(0) = 0.

    Returns
    -------
    ReducedKernelActionDatum
        The paper's reduced kernel datum: germ Φ on R, residual cone K_R,
        and normal prefactor.
    """
    d = seed.kernel_dim + seed.positive_normal_dim
    if action_kn_jet.dim != d:
        raise InputValidationError(
            f"action_kn_jet.dim={action_kn_jet.dim} does not match "
            f"kernel_dim + normal_dim = {d}"
        )

    return _reduce_kernel_action_from_kn(seed, action_kn_jet, metadata=metadata)


def _reduce_kernel_action_from_kn(
    seed: QuadraticKernelSeed,
    action_in_kn: PolynomialJet,
    *,
    metadata: Mapping[str, Any] | None = None,
) -> ReducedKernelActionDatum:
    """Shared core: perform the reduction given a jet in kn-coordinates.

    Both reduce_kernel_action (active → kn rotation done externally) and
    reduce_kernel_action_kn (jet already in kn-coords) delegate here.
    """
    r = seed.kernel_dim
    p = seed.positive_normal_dim
    d = r + p

    # ── Implicit function theorem: solve ∂_y Ψ(u, η(u)) = 0 ────────
    # In the (u, y) system, the Hessian ∂²Ψ/∂y² at origin is H_N (SPD).
    # By IFT, for small u, the normal equation ∂_y Ψ(u, y) = 0 has a
    # unique solution y = η(u) with η(0) = 0, Dη(0) = 0.
    #
    # For a polynomial jet, we compute η perturbatively:
    # η(u) = η₂(u) + η₃(u) + ... where η_k is homogeneous of degree k.
    # The leading term η₂ comes from matching the cubic terms in Ψ.

    # Extract H_N from the action in (u, y) coordinates.
    H_kn = action_in_kn.hessian_at_origin()
    H_NN = H_kn[r:, r:]  # shape (p, p) — should be SPD
    H_NN_inv = np.linalg.inv(H_NN)

    # Compute the reduced kernel action Φ by eliminating y perturbatively.
    kernel_action = _compute_reduced_kernel_action(
        action_in_kn, r, p, H_NN, H_NN_inv
    )

    return ReducedKernelActionDatum(
        kernel_dim=r,
        kernel_action=kernel_action,
        kernel_cone=seed.kernel_cone,
        positive_normal_dim=p,
        log_normal_prefactor=seed.log_normal_prefactor,
        total_dim=d,
        metadata=metadata,
    )


# ═══════════════════════════════════════════════════════════════════════
# Internal: perturbative kernel action computation
# ═══════════════════════════════════════════════════════════════════════

def _compute_reduced_kernel_action(
    action: PolynomialJet,
    r: int,
    p: int,
    H_NN: Array,
    H_NN_inv: Array,
) -> PolynomialJet:
    """Compute Φ(u) = Ψ(u, η(u)) by perturbative elimination.

    We work in coordinates z = (u₁,...,uᵣ, y₁,...,yₚ) where the first r
    are kernel and the last p are normal.

    The critical graph η(u) satisfies ∂_y Ψ(u, η(u)) = 0.  We solve
    this order by order:

    At quadratic order: η₂ = 0 (since ∂²Ψ/∂u∂y at origin is zero
    by the orthogonality of kernel and normal eigenvectors).

    At cubic order: η₂(u) corrections from cubic mixed terms.

    The reduced action to the certified order collects all pure-u terms
    and the corrections from substituting η.
    """
    d = r + p
    certified = action.certified_order

    # Step 1: Collect terms grouped by normal degree.
    # For each term c * u^α * y^β, the normal degree is |β|.
    pure_kernel_terms: dict[tuple[int, ...], float] = {}
    mixed_terms_by_normal_deg: dict[int, list[tuple[tuple[int, ...], tuple[int, ...], float]]] = {}

    for multi_idx, coeff in action.terms.items():
        u_idx = multi_idx[:r]
        y_idx = multi_idx[r:]
        normal_deg = sum(y_idx)

        if normal_deg == 0:
            # Pure kernel term — goes directly into Φ.
            pure_kernel_terms[u_idx] = pure_kernel_terms.get(u_idx, 0.0) + coeff
        else:
            mixed_terms_by_normal_deg.setdefault(normal_deg, []).append(
                (u_idx, y_idx, coeff)
            )

    # Step 2: Compute η corrections order by order.
    # η₁ = 0 (∇Ψ(0) = 0)
    # η₂: from ∂²Ψ/∂u∂y = 0 at origin (eigenvector orthogonality) → η₂ = 0
    #
    # η comes from solving: for each normal component j,
    #   ∂Ψ/∂yⱼ(u, η(u)) = 0
    #
    # The leading non-trivial contribution to Φ from η comes at order 4:
    # it involves cubic mixed terms (u²y) feeding into η₂, then
    # contributing to Φ at quartic order.
    #
    # Specifically, if Ψ has terms c_{αβ} u^α y^β with |β|=1 and |α|=2,
    # these define the mixed cubic coupling.  Write V_j(u) for the
    # coefficient of y_j in the cubic part:
    #   V_j(u) = Σ_{|α|=2} c_{α, e_j} u^α
    # Then η₂,j(u) = -Σ_k (H_NN^{-1})_{jk} V_k(u).
    # And the quartic correction to Φ is:
    #   ΔΦ₄(u) = -(1/2) Σ_{j,k} (H_NN^{-1})_{jk} V_j(u) V_k(u)

    # Extract V_j(u): coefficients of y_j in cubic terms (degree 3, normal_deg=1).
    V_coeffs = {}  # V_coeffs[j][u_idx] = coefficient
    for u_idx, y_idx, coeff in mixed_terms_by_normal_deg.get(1, []):
        u_deg = sum(u_idx)
        if u_deg < 2:
            continue
        # Find which y_j this term involves.
        for j in range(p):
            if y_idx[j] == 1 and sum(y_idx) == 1:
                V_coeffs.setdefault(j, {})[u_idx] = (
                    V_coeffs.get(j, {}).get(u_idx, 0.0) + coeff
                )

    # Compute the quartic correction ΔΦ₄(u) = -(1/2) V^T H_NN^{-1} V.
    has_cubic_coupling = bool(V_coeffs)
    # Check for higher-order normal nonlinearities: terms with normal_deg >= 2
    # AND total degree > 2.  The pure y² term at the origin (normal_deg=2,
    # total_deg=2) is just the normal Hessian — fully accounted for by the
    # reduction.  Only terms like y³, u·y², u²·y², etc. (total degree > 2)
    # are genuine nonlinearities that may feed uncontrolled corrections.
    has_higher_normal_nonlinearity = False
    for nd, term_list in mixed_terms_by_normal_deg.items():
        if nd < 2:
            continue
        for u_idx, y_idx, coeff in term_list:
            total_deg = sum(u_idx) + sum(y_idx)
            if total_deg > 2 and abs(coeff) > 1e-15:
                has_higher_normal_nonlinearity = True
                break
        if has_higher_normal_nonlinearity:
            break

    if has_cubic_coupling and certified >= 4:
        quartic_correction = _compute_quartic_correction(
            V_coeffs, H_NN_inv, r, p
        )
        for u_idx, coeff in quartic_correction.items():
            pure_kernel_terms[u_idx] = pure_kernel_terms.get(u_idx, 0.0) + coeff

    # ── Certified order bookkeeping ────────────────────────────────────
    # The reducer computes:
    #   - All pure kernel terms (exact, any order)
    #   - The quartic correction from cubic coupling V_j(u) (order 4)
    #
    # It does NOT compute:
    #   - Sextic corrections from quartic coupling or η₃ terms
    #   - Any higher-order corrections from normal nonlinearities
    #
    # Therefore the output certified order must be downgraded when
    # uncontrolled normal nonlinearities are present.  The safe bound:
    #   - If no mixed terms at all: pure kernel restriction is exact,
    #     keep original certified_order.
    #   - If cubic coupling only (normal_deg=1 terms, no higher):
    #     quartic correction is exact, but sextic from η₂·(quartic
    #     coupling) is uncontrolled.  Safe output = min(certified, 4).
    #     Exception: if all cubic coupling terms have u-degree exactly 2
    #     and there are no normal-degree >= 2 terms, then the quartic
    #     correction IS the only correction, so output = min(certified, 4).
    #   - If any normal-degree >= 2 terms exist (quadratic-in-y or higher):
    #     these feed into η₂·(quadratic-in-y) at order 4+, and η₂²·(...)
    #     at order 6+.  The quartic correction from cubic coupling
    #     is still correct, but quartic contributions from normal_deg=2
    #     terms are uncontrolled.  Safe output = min(certified, 4) at best,
    #     but if the quadratic-in-y terms couple with u at total order ≤ 4,
    #     we must downgrade further.  Conservative: output = 4 if we
    #     computed the quartic correction, 2 otherwise.
    if not has_cubic_coupling and not has_higher_normal_nonlinearity:
        # Pure kernel restriction only — exact to input order.
        output_certified_order = certified
    elif has_cubic_coupling and not has_higher_normal_nonlinearity:
        # Quartic correction computed, no higher nonlinearities.
        # The quartic correction is complete.  But sextic terms from
        # η₂ interacting with higher-u-degree cubic terms may be missed.
        # Safe: certify to order 4.
        output_certified_order = min(certified, 4)
    else:
        # Higher normal nonlinearities present.  The quartic correction
        # from cubic coupling is correct, but we may be missing quartic
        # contributions from normal_deg=2 terms feeding through η.
        # Conservative: certify to order 4 only if quartic correction
        # was computed, otherwise order 2.
        if has_cubic_coupling and certified >= 4:
            output_certified_order = 4
        else:
            output_certified_order = 2

    # Clean up near-zero terms but preserve ALL terms regardless of degree.
    # Pure kernel terms are exact (they come from the original jet, not from
    # perturbative elimination).  Correction terms (quartic from cubic coupling)
    # are also exact to the order they were computed.  The certified_order
    # metadata documents the degree to which the jet is certified COMPLETE —
    # terms above that order may be present but are not guaranteed to include
    # all contributions at that degree.  Downstream consumers (the singular
    # classifier) use certified_order to decide what to trust.
    cleaned = {k: v for k, v in pure_kernel_terms.items() if abs(v) > 1e-15}

    return PolynomialJet(
        dim=r,
        terms=cleaned,
        certified_order=output_certified_order,
        metadata={"reduction": "kernel_reduction_prop_9_1",
                   "input_certified_order": certified,
                   "has_cubic_coupling": has_cubic_coupling,
                   "has_higher_normal_nonlinearity": has_higher_normal_nonlinearity},
    )


def _compute_quartic_correction(
    V_coeffs: dict[int, dict[tuple[int, ...], float]],
    H_NN_inv: Array,
    r: int,
    p: int,
) -> dict[tuple[int, ...], float]:
    """Compute the quartic correction ΔΦ₄(u) = -(1/2) V^T H_NN^{-1} V.

    V_j(u) = Σ_α v_{j,α} u^α, so
    ΔΦ₄(u) = -(1/2) Σ_{j,k} (H^{-1})_{jk} Σ_{α,β} v_{j,α} v_{k,β} u^{α+β}
    """
    correction: dict[tuple[int, ...], float] = {}

    for j in range(p):
        for k in range(p):
            h_inv_jk = H_NN_inv[j, k]
            if abs(h_inv_jk) < 1e-15:
                continue
            v_j = V_coeffs.get(j, {})
            v_k = V_coeffs.get(k, {})
            for alpha, c_alpha in v_j.items():
                for beta, c_beta in v_k.items():
                    # Combined index is α + β.
                    combined = tuple(a + b for a, b in zip(alpha, beta))
                    correction[combined] = (
                        correction.get(combined, 0.0)
                        - 0.5 * h_inv_jk * c_alpha * c_beta
                    )

    return correction


# ═══════════════════════════════════════════════════════════════════════
# Polynomial expansion helper
# ═══════════════════════════════════════════════════════════════════════

def _expand_term(
    multi_idx: tuple[int, ...],
    coeff: float,
    basis: Array,
    k: int,
    output: dict[tuple[int, ...], float],
) -> None:
    """Expand a single monomial c * x^α under substitution x = B @ u.

    Collects all resulting monomials in u into the output dict.
    """
    d = basis.shape[0]
    total_deg = sum(multi_idx)

    if total_deg == 0:
        idx = tuple(0 for _ in range(k))
        output[idx] = output.get(idx, 0.0) + coeff
        return

    # For each variable x_i with exponent e_i > 0,
    # x_i = Σ_j B[i,j] u_j, so x_i^{e_i} = (Σ_j B[i,j] u_j)^{e_i}.
    # The full monomial is the product of these.
    #
    # We compute this via iterated polynomial multiplication.
    # Start with the constant polynomial {(): 1.0} and multiply
    # by (Σ_j B[i,j] u_j)^{e_i} for each i with e_i > 0.

    current: dict[tuple[int, ...], float] = {tuple(0 for _ in range(k)): coeff}

    for i, e in enumerate(multi_idx):
        if e == 0:
            continue
        # Build the linear polynomial L_i(u) = Σ_j B[i,j] u_j.
        linear = {}
        for j in range(k):
            if abs(basis[i, j]) > 1e-15:
                idx = tuple(1 if m == j else 0 for m in range(k))
                linear[idx] = basis[i, j]

        # Raise L_i to the e-th power by repeated multiplication.
        power = {tuple(0 for _ in range(k)): 1.0}
        for _ in range(e):
            power = _poly_multiply(power, linear, k)

        # Multiply current by power.
        current = _poly_multiply(current, power, k)

    # Merge into output.
    for idx, val in current.items():
        if abs(val) > 1e-15:
            output[idx] = output.get(idx, 0.0) + val


def _poly_multiply(
    p1: dict[tuple[int, ...], float],
    p2: dict[tuple[int, ...], float],
    k: int,
) -> dict[tuple[int, ...], float]:
    """Multiply two polynomials represented as multi-index dicts."""
    result: dict[tuple[int, ...], float] = {}
    for idx1, c1 in p1.items():
        for idx2, c2 in p2.items():
            combined = tuple(a + b for a, b in zip(idx1, idx2))
            result[combined] = result.get(combined, 0.0) + c1 * c2
    return result
