"""
Typed local evidence dispatcher — 0.3.3 Technical Note §11.

The evidence dispatcher sits on top of the regime detector.  It takes
a RegimeClassification and returns a typed evidence result:

    RegimeClassification  →  EvidenceResult
        RegularInterior       →  RegularEvidenceTemplate
        RegularCone           →  RegularEvidenceTemplate (cone mass flagged)
        RegularQuotient       →  RegularEvidenceTemplate (orbit volume)
        RegularQuotientCone   →  RegularEvidenceTemplate (both)
        UnresolvedKernel      →  UnresolvedEvidenceRefusal
        IndefiniteStationaryPoint → UnsupportedEvidenceRefusal
        UnsupportedConeReduction → UnsupportedEvidenceRefusal

The dispatcher does NOT produce a single scalar score.  It produces
a typed template whose structural pieces the caller can inspect,
combine with the log-likelihood, or compare across models.

The dispatcher also provides a convenience entry point that composes
hidden elimination → regime detection → evidence dispatch in a single
call, for callers who want the full pipeline.

References
----------
0.3.3 Technical Note §7 (Corollary 7.1: typed local evidence template),
§11 (practical consequence: the 5-layer stack).
"""
from __future__ import annotations

import numpy as np

from .regime import classify_regime, classify_from_hessian
from .regime_types import (
    ConeKind,
    ConeSpec,
    EvidenceResult,
    IndefiniteStationaryPoint,
    OrbitSpec,
    ReducedLocalDatum,
    RegimeClassification,
    RegimeKind,
    RegularCone,
    RegularEvidenceTemplate,
    RegularInterior,
    RegularQuotient,
    RegularQuotientCone,
    SingularEvidenceTemplate,
    UnresolvedEvidenceRefusal,
    UnresolvedKernel,
    UnsupportedConeReduction,
    UnsupportedEvidenceRefusal,
)
from .singular_types import (
    SingularClassification,
    WeightedHomogeneousTemplate,
    BranchChannelTemplate,
    CriticalBranchTemplate,
    UnresolvedKernelJetRefusal,
)
from .validation import Tolerances, resolve_tolerances


# ═══════════════════════════════════════════════════════════════════════
# Primary API: dispatch_evidence
# ═══════════════════════════════════════════════════════════════════════

def dispatch_evidence(
    classification: RegimeClassification,
) -> EvidenceResult:
    """Dispatch a regime classification to a typed evidence result.

    This is the core evidence layer.  It interprets the structural
    invariants from the regime detector and produces a typed evidence
    template (for regular regimes) or a typed refusal (for unresolved
    or unsupported cases).

    Parameters
    ----------
    classification : RegimeClassification
        Output of classify_regime or classify_from_hessian.

    Returns
    -------
    EvidenceResult
        RegularEvidenceTemplate, UnresolvedEvidenceRefusal,
        or UnsupportedEvidenceRefusal.
    """
    if isinstance(classification, RegularInterior):
        return _regular_interior_evidence(classification)

    if isinstance(classification, RegularCone):
        return _regular_cone_evidence(classification)

    if isinstance(classification, RegularQuotient):
        return _regular_quotient_evidence(classification)

    if isinstance(classification, RegularQuotientCone):
        return _regular_quotient_cone_evidence(classification)

    if isinstance(classification, UnresolvedKernel):
        return UnresolvedEvidenceRefusal(
            kernel=classification.kernel,
            kernel_dim=classification.kernel.kernel_dim,
        )

    if isinstance(classification, IndefiniteStationaryPoint):
        return UnsupportedEvidenceRefusal(
            cone_kind=ConeKind.FULL_SPACE,
            reason=(
                f"indefinite Hessian ({classification.negative_dim} negative "
                f"eigenvalue(s), min λ = {classification.min_eigenvalue:.2e}): "
                f"stationary point is not a local minimum; "
                f"suggested action: {classification.suggested_action}"
            ),
            metadata={"indefinite_classification": classification},
        )

    if isinstance(classification, UnsupportedConeReduction):
        return UnsupportedEvidenceRefusal(
            cone_kind=classification.cone_kind,
            reason=classification.reason,
        )

    raise TypeError(
        f"Unknown RegimeClassification variant: {type(classification)}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Singular evidence dispatch (Layer 5 → evidence template)
# ═══════════════════════════════════════════════════════════════════════

def dispatch_singular_evidence(
    singular: SingularClassification,
) -> EvidenceResult:
    """Dispatch a singular classification to a typed evidence result.

    This is the Layer 5 → evidence bridge.  It takes the output of
    classify_singular_kernel and wraps it into a SingularEvidenceTemplate
    (for certified channels) or an UnresolvedEvidenceRefusal (for
    refusals).

    Parameters
    ----------
    singular : SingularClassification
        Output of classify_singular_kernel (Layer 5).

    Returns
    -------
    EvidenceResult
        SingularEvidenceTemplate for certified channels,
        UnresolvedEvidenceRefusal for unresolved kernel jets.
    """
    if isinstance(singular, UnresolvedKernelJetRefusal):
        # The kernel jet was not resolved.  Return an evidence-level
        # refusal that carries the structural data forward.
        return UnresolvedEvidenceRefusal(
            # We don't have the original QuadraticKernelSeed here,
            # so we build a minimal stub.  The caller retains the
            # full singular classification for inspection.
            kernel=None,  # type: ignore[arg-type]
            kernel_dim=singular.kernel_dim,
            reason=(
                f"singular kernel jet unresolved: {singular.reason}"
            ),
            metadata={"singular_refusal": singular},
        )

    # Certified singular channel — wrap in SingularEvidenceTemplate.
    if isinstance(singular, WeightedHomogeneousTemplate):
        return SingularEvidenceTemplate(
            regime=RegimeKind.SINGULAR_WH,
            active_dim=singular.kernel_dim + singular.positive_normal_dim,
            kernel_dim=singular.kernel_dim,
            positive_normal_dim=singular.positive_normal_dim,
            log_normal_prefactor=singular.log_normal_prefactor,
            singular_classification=singular,
            metadata={"channel": "weighted_homogeneous"},
        )

    if isinstance(singular, CriticalBranchTemplate):
        return SingularEvidenceTemplate(
            regime=RegimeKind.SINGULAR_CRITICAL_BRANCH,
            active_dim=singular.kernel_dim + singular.positive_normal_dim,
            kernel_dim=singular.kernel_dim,
            positive_normal_dim=singular.positive_normal_dim,
            log_normal_prefactor=singular.log_normal_prefactor,
            singular_classification=singular,
            metadata={"channel": singular.kind.value},
        )

    if isinstance(singular, BranchChannelTemplate):
        return SingularEvidenceTemplate(
            regime=RegimeKind.SINGULAR_BRANCH,
            active_dim=singular.kernel_dim + singular.positive_normal_dim,
            kernel_dim=singular.kernel_dim,
            positive_normal_dim=singular.positive_normal_dim,
            log_normal_prefactor=singular.log_normal_prefactor,
            singular_classification=singular,
            metadata={"channel": singular.kind.value},
        )

    raise TypeError(
        f"Unknown SingularClassification variant: {type(singular)}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Convenience: classify + dispatch in one call
# ═══════════════════════════════════════════════════════════════════════

def local_evidence(
    datum: ReducedLocalDatum,
    tolerances: Tolerances | None = None,
) -> EvidenceResult:
    """Classify regime and dispatch to evidence in one call.

    Composes: regime detection → evidence dispatch.

    Parameters
    ----------
    datum : ReducedLocalDatum
        Pre-classification reduced local object.
    tolerances : Tolerances, optional

    Returns
    -------
    EvidenceResult
    """
    classification = classify_regime(datum, tolerances=tolerances)
    return dispatch_evidence(classification)


def local_evidence_from_hessian(
    H: np.ndarray,
    cone: ConeSpec | None = None,
    orbit: OrbitSpec | None = None,
    tolerances: Tolerances | None = None,
) -> EvidenceResult:
    """Classify regime from a raw Hessian and dispatch to evidence.

    Composes: datum construction → regime detection → evidence dispatch.

    Parameters
    ----------
    H : array, shape (d, d)
        Symmetric PSD Hessian of the visible action at the optimum.
        For quotient cases (orbit is not None), this must be the
        transverse Hessian H_⊥ on the symmetry slice.
    cone : ConeSpec, optional
        Admissible tangent cone.  Defaults to FULL_SPACE.
    orbit : OrbitSpec, optional
        Orbit/slice data.  Defaults to None.
    tolerances : Tolerances, optional

    Returns
    -------
    EvidenceResult
    """
    classification = classify_from_hessian(
        H, cone=cone, orbit=orbit, tolerances=tolerances
    )
    return dispatch_evidence(classification)


# ═══════════════════════════════════════════════════════════════════════
# Internal dispatch helpers
# ═══════════════════════════════════════════════════════════════════════

def _log_quadratic_det(active_dim: int, log_det_h: float) -> float:
    """Compute (m/2) log(2π) - (1/2) log det(H_active).

    This is the ordinary Laplace determinant contribution to the
    log local evidence.
    """
    return 0.5 * active_dim * np.log(2.0 * np.pi) - 0.5 * log_det_h


def _regular_interior_evidence(
    c: RegularInterior,
) -> RegularEvidenceTemplate:
    """Evidence for the regular interior regime.

    Z_loc ≈ exp(-S_vis) · (2π)^{m/2} · det(H)^{-1/2}
    """
    m = c.datum.active_dim
    return RegularEvidenceTemplate(
        regime=RegimeKind.REGULAR_INTERIOR,
        active_dim=m,
        log_quadratic_det_term=_log_quadratic_det(m, c.log_det_h_active),
    )


def _regular_cone_evidence(
    c: RegularCone,
) -> RegularEvidenceTemplate:
    """Evidence for the regular cone regime.

    Z_loc ≈ exp(-S_vis) · (2π)^{m/2} · det(H)^{-1/2} · P(Z ∈ K)

    The cone mass P(Z ∈ K) depends on the cone kind:
    - ORTHANT: computable via multivariate normal CDF (not yet implemented
      in 0.3.3 — flagged as None with cone_mass_known_exact=True since the
      formula is known).
    - ABSTRACT: not computable (None, cone_mass_known_exact=True since the
      limitation is extraction, not approximation).
    - POLYHEDRAL (full-span): same situation as ORTHANT.

    In all cases the determinant term is exact.  The cone mass is flagged
    as needing separate computation.
    """
    m = c.datum.active_dim
    cone_kind = c.datum.cone.kind

    # The cone mass is a separate numerical problem.
    # For 0.3.3, we return None (not yet computed) and flag
    # cone_mass_known_exact=True (the formula is known to be exact;
    # computation is deferred to a future auxiliary layer).
    log_cone_mass = None
    cone_mass_known_exact = True

    if cone_kind == ConeKind.ORTHANT:
        # For orthant cones, the Gaussian mass is computable via
        # multivariate normal CDF.  Deferred to a future auxiliary
        # layer — not part of the core exact theory.
        log_cone_mass = None

    return RegularEvidenceTemplate(
        regime=RegimeKind.REGULAR_CONE,
        active_dim=m,
        log_quadratic_det_term=_log_quadratic_det(m, c.log_det_h_active),
        log_cone_mass=log_cone_mass,
        cone_mass_known_exact=cone_mass_known_exact,
    )


def _regular_quotient_evidence(
    c: RegularQuotient,
) -> RegularEvidenceTemplate:
    """Evidence for the regular quotient regime.

    Z_loc ≈ exp(-S_vis) · Vol(G·v̂) · (2π)^{m/2} · det(H_⊥)^{-1/2}

    where m = active_dim = transverse dimension after slice reduction.

    The orbit volume and slice Jacobian are carried from the OrbitSpec
    as external multiplicative factors.  The Gaussian integral
    (determinant and (2π) prefactor) operates entirely on the
    transverse slice, so gaussian_dim = active_dim.

    CONVENTION: The ReducedLocalDatum must already carry the transverse
    Hessian H_⊥ and transverse dimension m as active_dim.  The orbit
    factor multiplies from outside — it is not mixed into the
    determinant computation.
    """
    m = c.datum.active_dim
    orbit = c.datum.orbit

    log_orbit_vol = None
    log_slice_jac = None
    if orbit is not None:
        log_orbit_vol = orbit.log_orbit_volume
        log_slice_jac = orbit.log_slice_jacobian

    return RegularEvidenceTemplate(
        regime=RegimeKind.REGULAR_QUOTIENT,
        active_dim=m,
        log_quadratic_det_term=_log_quadratic_det(m, c.log_det_h_active),
        log_orbit_volume=log_orbit_vol,
        log_slice_jacobian=log_slice_jac,
    )


def _regular_quotient_cone_evidence(
    c: RegularQuotientCone,
) -> RegularEvidenceTemplate:
    """Evidence for the regular quotient-cone regime.

    Z_loc ≈ exp(-S_vis) · Vol(G·v̂) · (2π)^{m/2}
          · det(H_{⊥,K})^{-1/2} · P(Z ∈ K_⊥)

    where m = active_dim = transverse dimension after slice reduction.
    Same convention as _regular_quotient_evidence: the datum carries the
    already-reduced transverse Hessian, and the orbit factor is external.
    """
    m = c.datum.active_dim
    orbit = c.datum.orbit

    log_orbit_vol = None
    log_slice_jac = None
    if orbit is not None:
        log_orbit_vol = orbit.log_orbit_volume
        log_slice_jac = orbit.log_slice_jacobian

    return RegularEvidenceTemplate(
        regime=RegimeKind.REGULAR_QUOTIENT_CONE,
        active_dim=m,
        log_quadratic_det_term=_log_quadratic_det(m, c.log_det_h_active),
        log_cone_mass=None,  # deferred, same as regular cone
        cone_mass_known_exact=True,
        log_orbit_volume=log_orbit_vol,
        log_slice_jacobian=log_slice_jac,
    )
