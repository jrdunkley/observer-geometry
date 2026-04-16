"""
Certified singular classifier — 0.3.3 Technical Note §10–12.

This module implements Layer 5 of the 5-layer stack: given a
ReducedKernelActionDatum from kernel reduction (Layer 4), classify the
kernel germ into one of two certified singular channels or refuse.

Channel 1: Weighted-homogeneous principal germs (Proposition 10.1)
    The kernel action Φ admits a positive weighted-homogeneous principal
    part P with anisotropic weights (ω₁, ..., ωᵣ).  The evidence law is
    n^{-|ω|} · ∫_{K_R} e^{-P(s)} ds.

Channel 2: Branch-restiffening channels (Theorems 11.1–11.3)
    For 2D kernel problems, a moving fibre minimum with gradually
    returning transverse stiffness.  The critical quadratic branch
    channel (c=2, quadratic vanishing) yields n^{-1/2} κ log n.

Refusal: when neither channel can be certified, return an honest
UnresolvedKernelJetRefusal.

References
----------
0.3.3 Technical Note:
  - §10  Definition 10.1 + Proposition 10.1
  - §11  Definition 11.1 + Theorems 11.1–11.3
  - §12  what remains unresolved
"""
from __future__ import annotations

from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from .exceptions import InputValidationError
from .kernel_reduction import PolynomialJet, ReducedKernelActionDatum
from .regime_types import ConeKind
from .singular_types import (
    BranchChannelDatum,
    BranchChannelTemplate,
    CriticalBranchTemplate,
    SingularClassification,
    SingularRegimeKind,
    UnresolvedKernelJetRefusal,
    WeightedHomogeneousTemplate,
)

Array = NDArray[np.float64]


# ═══════════════════════════════════════════════════════════════════════
# Primary API: classify_singular_kernel
# ═══════════════════════════════════════════════════════════════════════

def classify_singular_kernel(
    datum: ReducedKernelActionDatum,
    *,
    branch_straightened: bool = False,
    metadata: Mapping[str, Any] | None = None,
) -> SingularClassification:
    """Classify a reduced kernel action into a certified singular channel.

    Attempts Channel 1 (weighted-homogeneous) first, then Channel 2
    (branch-restiffening, for 2D kernels).  If neither succeeds,
    returns an honest refusal.

    Parameters
    ----------
    datum : ReducedKernelActionDatum
        Output of reduce_kernel_action (Layer 4).
    branch_straightened : bool
        Whether the kernel jet has been expressed in straightened branch
        coordinates (i.e. the branch curve gamma(v) has been analytically
        straightened to gamma = 0 via a coordinate change).  The branch
        channel classifier (Channel 2) assumes this.  If False (default),
        the branch channel is not attempted and 2D kernels that fail WH
        will fall through to refusal.  Set True only when the caller has
        verified the branch chart is straightened.

    Returns
    -------
    SingularClassification
        One of: WeightedHomogeneousTemplate, BranchChannelTemplate,
        CriticalBranchTemplate, or UnresolvedKernelJetRefusal.
    """
    r = datum.kernel_dim
    jet = datum.kernel_action

    # ── Channel 1: Weighted-homogeneous principal germ ──────────────
    wh_result = _try_weighted_homogeneous(jet, datum)
    if wh_result is not None:
        return wh_result

    # ── Channel 2: Branch-restiffening (2D kernels only) ────────────
    # Requires straightened branch coordinates as a precondition.
    if r == 2 and branch_straightened:
        branch_result = _try_branch_channel(jet, datum)
        if branch_result is not None:
            return branch_result

    # ── Refusal ─────────────────────────────────────────────────────
    refusal_text = _refusal_reason(jet, datum)
    if r == 2 and not branch_straightened:
        refusal_text += (
            "; branch channel not attempted (branch_straightened=False)"
        )
    return UnresolvedKernelJetRefusal(
        kernel_dim=r,
        reason=refusal_text,
        log_normal_prefactor=datum.log_normal_prefactor,
        positive_normal_dim=datum.positive_normal_dim,
        metadata=metadata,
    )


# ═══════════════════════════════════════════════════════════════════════
# Channel 1: Weighted-homogeneous classification
# ═══════════════════════════════════════════════════════════════════════

def _try_weighted_homogeneous(
    jet: PolynomialJet,
    datum: ReducedKernelActionDatum,
) -> WeightedHomogeneousTemplate | None:
    """Try to certify the kernel germ as weighted-homogeneous.

    For each variable uᵢ, find the lowest pure-power degree dᵢ
    (the lowest exponent in a term c·uᵢ^{dᵢ} with c > 0 and no
    other variables involved).  Set ωᵢ = 1/dᵢ.

    Then verify that P(u) = Σ cᵢ uᵢ^{dᵢ} (the diagonal principal
    part) is positive on K_R \\{0}, and that all other terms have
    weighted degree ≥ 1 (are subordinate to the principal part).

    This is the separable (product-type) case of Proposition 10.1.
    For non-separable weighted-homogeneous germs, a more general
    Newton polyhedron analysis would be needed.
    """
    r = jet.dim
    if r == 0:
        return None

    # Step 1: For each variable, find the lowest pure-power degree.
    pure_degrees: dict[int, tuple[int, float]] = {}  # var_idx → (degree, coeff)

    for multi_idx, coeff in jet.terms.items():
        # Check if this is a pure power in exactly one variable.
        nonzero = [(i, e) for i, e in enumerate(multi_idx) if e > 0]
        if len(nonzero) != 1:
            continue
        var_idx, deg = nonzero[0]
        if deg < 2:
            continue  # degree-1 terms are zero at a critical point
        if coeff <= 0:
            continue  # need positive coefficient for the principal part
        if var_idx not in pure_degrees or deg < pure_degrees[var_idx][0]:
            pure_degrees[var_idx] = (deg, coeff)

    # All variables must have a pure power.
    if len(pure_degrees) != r:
        return None

    degrees = tuple(pure_degrees[i][0] for i in range(r))
    coeffs = tuple(pure_degrees[i][1] for i in range(r))

    # ── Certified-order enforcement ──────────────────────────────────
    # The principal part lives at degree max(degrees).  If the jet's
    # certified_order is below that, we cannot trust that the principal
    # coefficient is complete (perturbative corrections at that degree
    # may be missing).  Refuse rather than over-certify.
    max_degree = max(degrees)
    if max_degree > jet.certified_order:
        return None

    # On FULL_SPACE cones, odd-degree pure powers are NOT positive
    # (e.g. u³ < 0 for u < 0).  Require even degrees for positivity.
    if datum.kernel_cone.kind == ConeKind.FULL_SPACE:
        if any(d % 2 != 0 for d in degrees):
            return None
    weights = tuple(1.0 / d for d in degrees)
    total_weight = sum(weights)

    # Step 2: Verify strict subordination — all non-principal terms must
    # have weighted degree strictly greater than 1.
    #
    # CRITICAL: mixed terms at weighted degree exactly 1 are part of the
    # principal germ (they live on the Newton boundary).  The separable
    # product integral formula Π_i ∫ e^{-c_i |s_i|^{d_i}} ds_i is only
    # licensed when the principal part is a sum of pure powers.  If any
    # mixed monomial sits at weighted degree 1, the principal germ is
    # non-separable and we cannot certify with the product formula.
    #
    # So we require: non-principal terms have weighted degree > 1 + ε,
    # and any term at weighted degree ≈ 1 that is NOT a principal
    # diagonal term causes refusal.
    for multi_idx, coeff in jet.terms.items():
        if abs(coeff) < 1e-15:
            continue
        weighted_deg = sum(e * w for e, w in zip(multi_idx, weights))
        # Identify the principal diagonal terms: pure power c_i u_i^{d_i}.
        is_principal_diagonal = False
        nonzero = [(i, e) for i, e in enumerate(multi_idx) if e > 0]
        if len(nonzero) == 1:
            var_idx, deg = nonzero[0]
            if var_idx in pure_degrees and deg == pure_degrees[var_idx][0]:
                is_principal_diagonal = True
        if is_principal_diagonal:
            continue  # this is part of the separable principal part

        # Non-principal term: must have weighted degree strictly > 1.
        # A mixed term at weighted degree = 1 is on the Newton boundary
        # and makes the principal germ non-separable → refuse.
        if weighted_deg < 1.0 + 1e-10:
            return None  # subordination fails or principal germ is non-separable

    # Step 3: Verify positivity on the active kernel cone.
    # For full-space cones, positivity of all principal coefficients suffices
    # (since the principal part is a sum of positive pure powers).
    # For non-trivial cones, additional verification is needed.
    if datum.kernel_cone.kind != ConeKind.FULL_SPACE:
        # For non-full-space cones with separable principal parts,
        # positivity still holds if all coefficients are positive.
        # But for restricted cones, some variables may not contribute.
        # For now, accept only FULL_SPACE or ABSTRACT.
        if datum.kernel_cone.kind == ConeKind.ABSTRACT:
            pass  # Trust the abstract structure
        else:
            return None  # Certification on restricted cones not yet implemented

    # Step 4: Compute the principal kernel integral (for FULL_SPACE).
    log_principal_integral = None
    if datum.kernel_cone.kind == ConeKind.FULL_SPACE:
        # For P(u) = Σ cᵢ |uᵢ|^{dᵢ} on all of R^r, the integral factors:
        # ∫ e^{-P(s)} ds = Π_i ∫_{-∞}^{∞} e^{-cᵢ |sᵢ|^{dᵢ}} dsᵢ
        # Each factor = 2 · Γ(1/dᵢ) / (dᵢ · cᵢ^{1/dᵢ})
        from math import lgamma
        log_integral = 0.0
        for i in range(r):
            d_i, c_i = degrees[i], coeffs[i]
            # ∫_{-∞}^{∞} e^{-c·|s|^d} ds = 2·Γ(1/d) / (d·c^{1/d})
            log_factor = (
                np.log(2.0) + lgamma(1.0 / d_i)
                - np.log(float(d_i)) - (1.0 / d_i) * np.log(c_i)
            )
            log_integral += log_factor
        log_principal_integral = log_integral

    return WeightedHomogeneousTemplate(
        kernel_dim=r,
        weights=weights,
        total_weight=total_weight,
        log_principal_integral=log_principal_integral,
        log_normal_prefactor=datum.log_normal_prefactor,
        positive_normal_dim=datum.positive_normal_dim,
        leading_monomial_degrees=degrees,
        metadata={"principal_coeffs": coeffs},
    )


# ═══════════════════════════════════════════════════════════════════════
# Channel 2: Branch-restiffening classification (2D only)
# ═══════════════════════════════════════════════════════════════════════

def _try_branch_channel(
    jet: PolynomialJet,
    datum: ReducedKernelActionDatum,
) -> BranchChannelTemplate | CriticalBranchTemplate | None:
    """Try to certify a 2D kernel germ as a branch-restiffening channel.

    Assumes coordinates (u, v) = (x, y) on the 2D kernel space.

    Step 1: Find the branch.  Look for a one-dimensional fibre minimum:
            ∂_u Φ(γ(v), v) = 0.  For polynomial Φ with ∂_u Φ(0,0) = 0,
            the branch starts at the origin.

    Step 2: Determine the branch order c (first nonzero ∂ᵤᶜ at origin).

    Step 3: Compute the stiffness λ(v) = (1/c!) ∂ᵤᶜ Φ̃(0,v).

    Step 4: Determine the stiffness vanishing order.

    Step 5: For c=2, check the critical quadratic branch channel conditions.
    """
    if jet.dim != 2:
        return None
    if datum.kernel_cone.kind not in (ConeKind.FULL_SPACE, ConeKind.ABSTRACT):
        return None  # Branch channels on restricted cones not yet supported

    # Coordinate convention: u = variable 0, v = variable 1.

    # Step 1: The branch is at x = 0 (since ∂_u Φ(0,0) = 0 by critical point).
    # For a polynomial, the straightened branch is trivially at u = 0 if
    # the linear and quadratic u-terms vanish.  More generally, we need
    # γ(v) from solving ∂_u Φ(γ(v), v) = 0.
    # For simplicity, we assume the branch has been straightened (γ = 0).
    # This is valid when the jet is already expressed in straightened coords.

    # ── Certified-order enforcement ──────────────────────────────────
    # The branch channel requires reading coefficients at degrees
    # up to the branch action degree d and higher transverse degree a.
    # If certified_order < 4, the jet is too incomplete for branch
    # classification (we need at least quartic structure).
    if jet.certified_order < 4:
        return None

    # Step 2: Find branch order c.
    # c is the smallest integer ≥ 2 such that the coefficient of u^c,
    # viewed as a polynomial in v, is not identically zero.  That is,
    # there exists at least one term (c, j) for some j >= 0 with
    # nonzero coefficient.  This is weaker than requiring a pure u^c
    # term at v=0: the stiffness may vanish at v=0 but be non-trivial
    # as a function of v, which is exactly the branch-restiffening
    # scenario (Definition 11.1).
    c = None
    for deg in range(2, jet.certified_order + 1):
        for multi_idx, coeff in jet.terms.items():
            if multi_idx[0] == deg and abs(coeff) > 1e-15:
                c = deg
                break
        if c is not None:
            break

    if c is None:
        return None  # No finite branch order found within certified order

    # Step 3: Compute the stiffness λ(v) = (1/c!) ∂ᵤᶜ Φ̃(0,v).
    # In our polynomial, this is the coefficient of u^c as a function of v:
    # λ(v) = Σ_j coeff_{(c, j)} · v^j
    stiffness_coeffs = {}  # j → coefficient of v^j in λ(v)
    for multi_idx, coeff in jet.terms.items():
        u_exp, v_exp = multi_idx
        if u_exp == c:
            stiffness_coeffs[v_exp] = stiffness_coeffs.get(v_exp, 0.0) + coeff

    # Step 4: Stiffness vanishing order.
    # λ(0) = c_coeff (which may or may not be nonzero).
    # If λ(0) ≠ 0, the stiffness does not vanish and we don't have a
    # branch-restiffening channel in the singular sense.
    # Wait — actually, for the CRITICAL quadratic branch, λ(v) is the
    # coefficient of u² as function of v, and it must vanish as v → 0
    # (i.e. the Hessian ∂²ᵤΦ̃(0,v) vanishes at v=0).
    #
    # In our setup, the pure u^c term at v=0 is C₀ (with c potentially > 2),
    # and the stiffness is ∂²ᵤΦ̃(0,v) = Σ coeff_{(2,j)} · (2·1) · v^j
    # for c=2.  But we're looking at the general case.
    #
    # For branch-restiffening, we need the u² coefficient (when c=2)
    # to vanish at v=0 but grow with v.

    if c == 2:
        # The u² coefficient as function of v.
        H_v_coeffs = {}  # j → coefficient of v^j in (1/2)∂²ᵤΦ̃(0,v)
        for multi_idx, coeff in jet.terms.items():
            u_exp, v_exp = multi_idx
            if u_exp == 2:
                H_v_coeffs[v_exp] = H_v_coeffs.get(v_exp, 0.0) + coeff

        # Check if the u² coefficient vanishes at v=0.
        h_at_zero = H_v_coeffs.get(0, 0.0)
        if abs(h_at_zero) > 1e-15:
            # The Hessian ∂²ᵤΦ is nonzero at origin → this is handled by
            # the normal-direction reduction, not a branch-restiffening channel.
            # (The u-direction is not truly degenerate.)
            return None

        # Find the stiffness vanishing order: smallest j > 0 with H_v_coeffs[j] ≠ 0.
        stiffness_ord = None
        H_0 = None
        for j in sorted(H_v_coeffs.keys()):
            if j == 0:
                continue
            if abs(H_v_coeffs[j]) > 1e-15:
                stiffness_ord = j
                H_0 = H_v_coeffs[j]
                break

        if stiffness_ord is None or H_0 is None or H_0 <= 0:
            return None  # Stiffness doesn't return positively

        # Step 5: Check critical quadratic branch conditions.
        # Need:
        #   (1/2)∂²ᵤΦ̃(0,v) = H₀ v^{stiffness_ord} + O(v^{stiffness_ord+1})
        #   Φ̃(0,v) = B₀ v^d + O(v^{d+1})  with B₀ > 0
        #   (1/a!)∂ᵤᵃΦ̃(0,0) = C₀ > 0  for some a > 2

        # The requirement for critical quadratic branch is stiffness_ord = 2.
        if stiffness_ord != 2:
            # General branch-restiffening but not critical quadratic.
            branch_datum = BranchChannelDatum(
                branch_order=2,
                stiffness_vanishing_order=stiffness_ord,
                branch_channel_exponent=stiffness_ord / 2.0,
            )
            return BranchChannelTemplate(
                kernel_dim=2,
                branch_datum=branch_datum,
                log_normal_prefactor=datum.log_normal_prefactor,
                positive_normal_dim=datum.positive_normal_dim,
            )

        # stiffness_ord = 2.  Now check the remaining conditions.

        # Find B₀, d: the branch action Φ̃(0,v) = B₀ v^d + ...
        B_0 = None
        d_branch = None
        for multi_idx, coeff in jet.terms.items():
            u_exp, v_exp = multi_idx
            if u_exp == 0 and v_exp >= 2 and abs(coeff) > 1e-15:
                if d_branch is None or v_exp < d_branch:
                    d_branch = v_exp
                    B_0 = coeff

        if d_branch is None or B_0 is None or B_0 <= 0:
            return None  # Branch action doesn't have positive leading term

        # Enforce: branch action degree must be within certified range.
        if d_branch > jet.certified_order:
            return None

        # Find C₀, a: higher-order transverse term at v=0.
        # (1/a!)∂ᵤᵃΦ̃(0,0) = C₀ > 0 for some a > 2.
        C_0 = None
        a_trans = None
        for deg in range(3, jet.certified_order + 1):
            idx = (deg, 0)
            coeff = jet.terms.get(idx, 0.0)
            if abs(coeff) > 1e-15 and coeff > 0:
                a_trans = deg
                C_0 = coeff
                break

        if a_trans is None or C_0 is None:
            return None  # No higher transverse term found

        # Compute κ = (d(a-2) - 2a) / (2a).
        kappa = (d_branch * (a_trans - 2) - 2 * a_trans) / (2.0 * a_trans)
        if kappa <= 0:
            return None  # Critical balance condition fails

        # ── Hard scope condition (TN4 Definition 11.3) ──────────────
        # The principal branch-restiffening face is:
        #     C₀ u^a + H₀ u² v² + B₀ v^d
        #
        # Weight system from the outer Newton edge (a,0)-(0,d):
        #     w(i,j) = i/a + j/d
        #
        # Hard scope: every non-principal monomial must have weighted
        # order strictly > 1 + κ.  This ensures no stray term can
        # interfere with the critical branch asymptotic.
        principal_exponents = {(a_trans, 0), (2, 2), (0, d_branch)}
        for multi_idx, coeff in jet.terms.items():
            if abs(coeff) < 1e-15:
                continue
            u_exp, v_exp = multi_idx
            if (u_exp, v_exp) in principal_exponents:
                continue
            weighted_order = u_exp / float(a_trans) + v_exp / float(d_branch)
            if weighted_order <= 1.0 + kappa + 1e-10:
                return None  # Hard scope violated by stray monomial

        # Build the critical quadratic branch template.
        branch_datum = BranchChannelDatum(
            branch_order=2,
            stiffness_vanishing_order=2,
            branch_channel_exponent=1.0,  # ord₀(λ)/c = 2/2 = 1
            branch_action_degree=d_branch,
            higher_transverse_degree=a_trans,
            kappa=kappa,
            stiffness_leading_coeff=H_0,
            branch_leading_coeff=B_0,
            higher_transverse_coeff=C_0,
        )

        # log coefficient = log(√π / d · H₀^{-1/2} · κ)
        log_coeff = (
            0.5 * np.log(np.pi)
            - np.log(float(d_branch))
            - 0.5 * np.log(H_0)
            + np.log(kappa)
        )

        return CriticalBranchTemplate(
            kernel_dim=2,
            branch_datum=branch_datum,
            log_coefficient=log_coeff,
            kappa=kappa,
            log_normal_prefactor=datum.log_normal_prefactor,
            positive_normal_dim=datum.positive_normal_dim,
        )

    else:
        # c > 2: general branch order, not quadratic.
        # Record the branch datum but do not attempt closed-form law.
        stiffness_ord = None
        for j in sorted(stiffness_coeffs.keys()):
            if j == 0:
                continue
            if abs(stiffness_coeffs.get(j, 0.0)) > 1e-15:
                stiffness_ord = j
                break

        if stiffness_ord is None:
            return None

        branch_datum = BranchChannelDatum(
            branch_order=c,
            stiffness_vanishing_order=stiffness_ord,
            branch_channel_exponent=stiffness_ord / float(c),
        )
        return BranchChannelTemplate(
            kernel_dim=2,
            branch_datum=branch_datum,
            log_normal_prefactor=datum.log_normal_prefactor,
            positive_normal_dim=datum.positive_normal_dim,
        )


# ═══════════════════════════════════════════════════════════════════════
# Refusal reason builder
# ═══════════════════════════════════════════════════════════════════════

def _refusal_reason(jet: PolynomialJet, datum: ReducedKernelActionDatum) -> str:
    """Build a human-readable refusal reason."""
    r = jet.dim
    reasons = []

    if r == 0:
        return "kernel dimension is zero"

    # Check if any pure powers exist.
    has_pure = {}
    for multi_idx, coeff in jet.terms.items():
        nonzero = [(i, e) for i, e in enumerate(multi_idx) if e > 0]
        if len(nonzero) == 1 and coeff > 0:
            i, e = nonzero[0]
            if i not in has_pure or e < has_pure[i]:
                has_pure[i] = e

    if len(has_pure) < r:
        missing = [i for i in range(r) if i not in has_pure]
        reasons.append(
            f"no positive pure-power term found for kernel variable(s) {missing}"
        )

    if r > 2:
        reasons.append(
            f"branch-restiffening channels require kernel_dim=2, got {r}"
        )

    if datum.kernel_cone.kind not in (ConeKind.FULL_SPACE, ConeKind.ABSTRACT):
        reasons.append(
            f"singular classification on {datum.kernel_cone.kind} cones "
            f"not yet supported"
        )

    if not reasons:
        reasons.append(
            "no certified singular channel matches the kernel jet structure"
        )

    return "; ".join(reasons)
