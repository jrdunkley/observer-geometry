"""
Typed local evidence geometry — core data objects.

These types implement the 0.3.3 Technical Note's architecture:

    local reduction  →  ReducedLocalDatum  →  detector  →  RegimeClassification
                                                             ├── RegularInterior
                                                             ├── RegularCone
                                                             ├── RegularQuotient
                                                             ├── RegularQuotientCone
                                                             └── UnresolvedKernel
                                                                   └── QuadraticKernelSeed

Design invariants
-----------------
1.  ReducedLocalDatum is the *pre-classification* object.  It carries the
    reduced active space, active Hessian (PSD, not forced SPD), cone spec,
    and optional orbit/slice data.  It does NOT carry a kernel datum — that
    is derived by the detector only when the reduced quadratic form is
    singular on the active object.

2.  The detector produces a RegimeClassification (tagged union).  Each
    variant carries the datum it was derived from, plus the structural
    invariants needed by a downstream evidence dispatcher.

3.  QuadraticKernelSeed is the output of the quadratic regime detector
    when the reduced Hessian has a nontrivial kernel.  It carries the
    kernel basis, residual cone, positive-normal eigenvalues, and normal
    prefactor — i.e. the quadratic-order decomposition E_K = R ⊕ N.
    This is the *input* to kernel reduction (Layer 4), not the paper's
    "reduced kernel datum" (which is the *output* of kernel reduction:
    the germ Φ, the residual cone K_R, and the normal prefactor).
    The name "seed" makes the distinction explicit.

4.  ConeSpec and OrbitSpec are value objects, not validation layers.  The
    regime detector and evidence dispatcher interpret them; the types
    themselves are inert carriers.

References
----------
0.3.3 Technical Note §7 (quadratic regime detector),
§8 (kernel reduction), §9 (anisotropic kernel law),
§10 (kernel-jet problem), §11 (practical consequence).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Union

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


# ═══════════════════════════════════════════════════════════════════════
# Cone specification
# ═══════════════════════════════════════════════════════════════════════

class ConeKind(str, Enum):
    """Classification of the admissible tangent cone at a visible optimum."""

    FULL_SPACE = "full_space"
    """No boundary constraint — the full tangent space is admissible."""

    ORTHANT = "orthant"
    """Non-negativity constraints on some coordinates.
    The admissible set is an orthant of the active span."""

    POLYHEDRAL = "polyhedral"
    """Linear inequality constraints: {x : A x >= 0} or {x : A x >= b}.
    Includes orthant as a special case but allows general half-space
    intersections."""

    ABSTRACT = "abstract"
    """Cone structure is known to exist but no concrete numerical
    representation is available.  The exact form of the law
    (determinant × cone mass) still holds, but the cone mass
    cannot yet be evaluated numerically.  This lets the theory
    proceed even when extraction is incomplete."""


@dataclass(frozen=True)
class ConeSpec:
    """Specification of the admissible tangent cone at a visible optimum.

    A tangent cone is homogeneous by definition: it is a cone in the
    tangent space at the optimum, not an affine feasible set.  Affine
    offsets belong to the pre-localisation feasible-set description,
    not to this object.

    This is a value object — it carries data but does not validate or
    compute anything.  The regime detector and evidence dispatcher
    interpret the cone; this type is an inert carrier.
    """

    kind: ConeKind
    ambient_dim: int

    # ── Optional concrete data for numerical layers ──
    #
    # For ORTHANT: signs[i] = +1 means x_i >= 0, -1 means x_i <= 0.
    # Length must equal ambient_dim.  None means the positive orthant.
    orthant_signs: Array | None = None

    # For POLYHEDRAL: the homogeneous cone {x : A x >= 0}.
    # A has shape (num_constraints, ambient_dim).
    # No affine offsets — tangent cones are homogeneous.
    halfspace_normals: Array | None = None

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Orbit / slice specification
# ═══════════════════════════════════════════════════════════════════════

class JacobianConvention(str, Enum):
    """Convention used for the slice-chart Jacobian in a quotient reduction."""

    NONE = "none"
    """No quotient reduction has been applied."""

    SLICE_LEBESGUE = "slice_lebesgue"
    """Slice chart with Lebesgue measure on the transverse directions.
    The Jacobian absorbs the volume distortion from the slice embedding."""

    SLICE_RIEMANNIAN = "slice_riemannian"
    """Slice chart with Riemannian measure inherited from the ambient
    Fisher-Rao metric.  The Jacobian includes the metric determinant."""

    USER_DECLARED = "user_declared"
    """User-supplied Jacobian value — the system trusts but does not
    derive the convention."""


@dataclass(frozen=True)
class OrbitSpec:
    """Specification of a compact orbit to be quotiented out.

    Carries the orbit dimension, Jacobian convention, and optional
    pre-computed orbit volume.  The slice Jacobian value is carried
    separately because it may depend on the specific chart, not just
    the orbit structure.
    """

    orbit_dim: int
    jacobian_convention: JacobianConvention
    log_orbit_volume: float | None = None
    log_slice_jacobian: float | None = None
    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Reduced local datum (PRE-classification)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ReducedLocalDatum:
    """The reduced local object after hidden elimination and any chart reduction.

    This is the pre-classification input to the regime detector.  It carries
    the reduced active space, the active Hessian (which may be PSD, not
    necessarily SPD — that is the whole point of regime detection), cone
    data if a boundary is present, and orbit data if a symmetry slice has
    been applied.

    CRITICAL CONVENTION (0.3.3 Technical Note §7):

    This object represents the **already-reduced** local datum.  For
    quotient cases (orbit is not None), the caller must have already
    performed the exact symmetry-slice reduction before constructing
    this object.  That means:

      - active_dim is the **transverse** dimension (d - orbit_dim),
        not the ambient visible dimension.
      - h_active is the **transverse Hessian** H_⊥ on the slice,
        shape (active_dim, active_dim), not the full visible Hessian.
      - OrbitSpec carries the orbit volume and slice Jacobian as
        separated external factors — they multiply the evidence law
        from outside the Gaussian integral.

    The evidence formula for the quotient regime is then:

        Z_loc ≈ exp(-S_vis) · Vol(G·v̂) · (2π)^{m/2} · det(H_⊥)^{-1/2}

    where m = active_dim = transverse dimension, and H_⊥ = h_active.
    If you pass a full (d × d) Hessian with orbit_dim > 0, the
    evidence formula will mix incompatible dimensions.

    The detector examines h_active, projects onto the cone span if needed,
    computes the reduced kernel, and returns a RegimeClassification.

    h_active is allowed to be PSD (not forced SPD).  If you force SPD at
    this layer, you destroy the purpose of the regime detector.
    """

    active_dim: int
    """Dimension of the active space after all reductions.  For quotient
    cases this is the transverse dimension (d - orbit_dim), not the
    ambient visible dimension."""

    h_active: Array
    """Symmetric PSD matrix on the active span, shape (active_dim, active_dim).
    For quotient cases this is the transverse Hessian H_⊥ on the slice.
    It may have a nontrivial kernel — that is what the detector classifies."""

    cone: ConeSpec
    """Admissible tangent cone at the visible optimum."""

    orbit: OrbitSpec | None = None
    """Orbit/slice data, present only when a compact symmetry has been
    quotiented out.  None means no quotient reduction was applied.
    When present, active_dim and h_active must already be the
    transverse (post-slice) quantities."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Reduced kernel datum (derived by detector, NOT pre-classification)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class QuadraticKernelSeed:
    """The quadratic-order kernel decomposition from the regime detector.

    Produced by the detector when the active Hessian has a nontrivial
    kernel on the active cone span.  This is the *input* to kernel
    reduction (Layer 4 of the 5-layer stack), not the paper's "reduced
    kernel datum" (which is the *output* of kernel reduction and carries
    the germ Φ, residual cone K_R, and normal prefactor).

    The seed carries:
    - the kernel subspace R = ker(H_K),
    - the positive normal complement N,
    - the eigenvalues and determinant prefactor for N,
    - the residual cone restricted to R.

    This is the E_K = R ⊕ N decomposition from Proposition 9.1,
    but without the actual kernel action germ Φ.  To obtain Φ you
    need the full local action representation, not just the Hessian.
    """

    kernel_dim: int
    """Dimension of the reduced kernel R_K = ker(H_K)."""

    kernel_basis: Array
    """Orthonormal basis for R_K, shape (active_dim, kernel_dim)."""

    kernel_cone: ConeSpec
    """The residual admissible cone restricted to the kernel directions.
    If the original cone was FULL_SPACE, this is also FULL_SPACE on the
    kernel.  If the original cone was ORTHANT or POLYHEDRAL, this is the
    intersection of that cone with the kernel subspace."""

    positive_normal_dim: int
    """Dimension of the positive normal complement N (= active_dim - kernel_dim)."""

    positive_normal_basis: Array
    """Orthonormal basis for the positive normal complement N,
    shape (active_dim, positive_normal_dim).  Carrying this avoids
    a recomputation if Layer 5 needs the normal projector for
    reduced germ extraction."""

    positive_normal_eigenvalues: Array
    """Eigenvalues of H_K restricted to the positive normal complement N,
    shape (positive_normal_dim,).  All positive."""

    log_normal_prefactor: float
    """The log of the normal-direction Laplace prefactor:
    (dim_N / 2) * log(2π) - (1/2) * sum(log(eigenvalues)).
    This is the exact contribution from integrating out the positive
    normal directions (Proposition 9.1)."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Regime classification (tagged union — detector output)
# ═══════════════════════════════════════════════════════════════════════

class RegimeKind(str, Enum):
    """Typed regime classification from the quadratic detector."""

    REGULAR_INTERIOR = "regular_interior"
    """Full space, SPD Hessian.  Ordinary Laplace determinant law."""

    REGULAR_CONE = "regular_cone"
    """Boundary cone, SPD Hessian on cone span.  Determinant × cone mass."""

    REGULAR_QUOTIENT = "regular_quotient"
    """Quotient slice, SPD transverse Hessian, no boundary.
    Orbit volume × transverse determinant."""

    REGULAR_QUOTIENT_CONE = "regular_quotient_cone"
    """Quotient slice with boundary cone on transverse directions.
    Orbit volume × transverse determinant × cone mass."""

    UNRESOLVED_KERNEL = "unresolved_kernel"
    """Nontrivial kernel in the active Hessian.  No purely quadratic
    rule can determine the leading local evidence law.  Higher-order
    data (jet, principal germ) is needed."""

    INDEFINITE_STATIONARY = "indefinite_stationary"
    """The active Hessian has materially negative eigenvalues.  The
    stationary point is not a local minimum — it is a saddle point or
    other non-minimum critical point.  The theorem-backed regime
    detector requires a PSD Hessian (local minimum), so this case is
    refused before regime classification.  Typical cause: EM or other
    optimiser converged to a saddle point."""

    UNSUPPORTED_REDUCTION = "unsupported_reduction"
    """The cone or slice structure requires a reduction step that is not
    yet implemented.  This is an implementation boundary, not a
    mathematical obstruction.  The datum is returned unchanged so the
    caller can inspect or fall back to ABSTRACT cone handling."""

    # ── Singular regime kinds (from Layer 5 classifier) ──────────────
    SINGULAR_WH = "singular_weighted_homogeneous"
    """Weighted-homogeneous kernel principal germ.  The evidence law is
    governed by the anisotropic kernel integral: Z ~ n^{-(p/2+|ω|)}."""

    SINGULAR_BRANCH = "singular_branch_restiffening"
    """General branch-restiffening channel.  Structural data certified
    but no closed-form asymptotic law yet available."""

    SINGULAR_CRITICAL_BRANCH = "singular_critical_quadratic_branch"
    """Critical quadratic branch channel (c=2, quadratic stiffness
    vanishing).  The evidence law carries a logarithmic correction:
    Z ~ n^{-(p/2+1/2)} · κ log n."""


@dataclass(frozen=True)
class RegularInterior:
    """Regular interior regime — ordinary Laplace determinant law.

    The active Hessian is SPD on the full space.  No boundary, no
    quotient.  The local evidence is:

        Z_loc ≈ exp(-S_vis) · (2π)^{m/2} · det(H)^{-1/2}
    """

    kind: RegimeKind = field(default=RegimeKind.REGULAR_INTERIOR, init=False)
    datum: ReducedLocalDatum = field(repr=False)
    log_det_h_active: float = 0.0
    """log det(H_active).  The ordinary Laplace determinant."""


@dataclass(frozen=True)
class RegularCone:
    """Regular cone regime — determinant × Gaussian cone mass.

    The active Hessian is SPD on the cone span.  The local evidence is:

        Z_loc ≈ exp(-S_vis) · (2π)^{m/2} · det(H)^{-1/2} · P(Z ∈ K)
    """

    kind: RegimeKind = field(default=RegimeKind.REGULAR_CONE, init=False)
    datum: ReducedLocalDatum = field(repr=False)
    log_det_h_active: float = 0.0
    """log det(H_active restricted to cone span)."""


@dataclass(frozen=True)
class RegularQuotient:
    """Regular quotient regime — orbit volume × transverse determinant.

    The transverse Hessian (after removing orbit directions) is SPD.
    No boundary cone.  The local evidence is:

        Z_loc ≈ exp(-S_vis) · Vol(G·v̂) · (2π)^{(d-r)/2} · det(H_⊥)^{-1/2}
    """

    kind: RegimeKind = field(default=RegimeKind.REGULAR_QUOTIENT, init=False)
    datum: ReducedLocalDatum = field(repr=False)
    log_det_h_active: float = 0.0
    """log det(H_⊥) — Hessian restricted to transverse slice."""
    orbit_dim: int = 0
    """Dimension of the compact orbit that has been quotiented out."""


@dataclass(frozen=True)
class RegularQuotientCone:
    """Regular quotient-cone regime — orbit × transverse determinant × cone mass.

    Both quotient and boundary are present.  The local evidence is:

        Z_loc ≈ exp(-S_vis) · Vol(G·v̂) · (2π)^{m/2}
              · det(H_{⊥,K})^{-1/2} · P(Z ∈ K_⊥)
    """

    kind: RegimeKind = field(default=RegimeKind.REGULAR_QUOTIENT_CONE, init=False)
    datum: ReducedLocalDatum = field(repr=False)
    log_det_h_active: float = 0.0
    """log det(H_{⊥,K}) on the transverse cone span."""
    orbit_dim: int = 0
    """Dimension of the compact orbit that has been quotiented out."""


@dataclass(frozen=True)
class UnresolvedKernel:
    """Unresolved kernel regime — quadratic data alone cannot determine the law.

    The active Hessian has a nontrivial kernel on the active cone span.
    The quadratic kernel seed carries the kernel basis, residual cone,
    and normal prefactor.  A downstream kernel reduction (Layer 4) and
    singular classifier (Layer 5) must inspect higher-order data.
    """

    kind: RegimeKind = field(default=RegimeKind.UNRESOLVED_KERNEL, init=False)
    datum: ReducedLocalDatum = field(repr=False)
    kernel: QuadraticKernelSeed = field(repr=False)


@dataclass(frozen=True)
class UnsupportedConeReduction:
    """Implementation-boundary refusal — the cone requires a reduction
    step that is not yet implemented in this version.

    This is distinct from both regular quadratic control (where the law
    is known) and unresolved kernel (where the mathematics says more data
    is needed).  Here the mathematics may well be tractable, but the
    software has not yet implemented the required reduction.

    The datum is returned unchanged so the caller can inspect the
    situation, fall back to ABSTRACT cone handling, or defer.
    """

    kind: RegimeKind = field(default=RegimeKind.UNSUPPORTED_REDUCTION, init=False)
    datum: ReducedLocalDatum = field(repr=False)

    cone_kind: ConeKind = ConeKind.POLYHEDRAL
    """Which cone kind triggered the refusal."""

    detected_span_dim: int = 0
    """The detected span dimension of the cone (less than active_dim)."""

    active_dim: int = 0
    """The active dimension the cone was expected to span."""

    reason: str = ""
    """Human-readable explanation of the implementation limitation."""

    workaround: str = (
        "Use ConeSpec(kind=ConeKind.ABSTRACT, ...) to proceed with "
        "regime classification only (cone mass will not be computable)."
    )
    """Suggested workaround for the caller."""


@dataclass(frozen=True)
class IndefiniteStationaryPoint:
    """The active Hessian has materially negative eigenvalues.

    This is not a local minimum of the visible action — the stationary
    point is a saddle or other non-minimum critical point.  The
    theorem-backed regime detector (which requires a PSD Hessian) does
    not apply.  This is a pre-detector refusal: the spectral analysis
    runs before regime classification.

    Typical causes:
      - EM converged to a saddle point (common for overfit mixtures)
      - Numerical optimisation did not reach a local minimum
      - The action landscape has no local minimum in this region

    The datum is returned so the caller can inspect the eigenstructure,
    refit with different initialisation, or discard the stationary point.
    """

    kind: RegimeKind = field(default=RegimeKind.INDEFINITE_STATIONARY, init=False)
    datum: ReducedLocalDatum = field(repr=False)

    negative_dim: int = 0
    """Number of eigenvalues that are materially negative."""

    zero_dim: int = 0
    """Number of eigenvalues within the spectral tolerance of zero."""

    positive_dim: int = 0
    """Number of eigenvalues that are materially positive."""

    min_eigenvalue: float = 0.0
    """Most negative eigenvalue of the active Hessian."""

    spectral_gap: float = 0.0
    """Gap between the most negative eigenvalue and the spectral tolerance.
    Large gap = clearly indefinite; small gap = borderline."""

    eigenvalues: Array | None = None
    """Full eigenvalue array for inspection, if carried."""

    suggested_action: str = "refit"
    """Suggested next step: 'refit', 'polish', or 'discard'."""

    metadata: Mapping[str, Any] | None = None


# Tagged union type for pattern matching / isinstance dispatch.
RegimeClassification = Union[
    RegularInterior,
    RegularCone,
    RegularQuotient,
    RegularQuotientCone,
    UnresolvedKernel,
    IndefiniteStationaryPoint,
    UnsupportedConeReduction,
]


# ═══════════════════════════════════════════════════════════════════════
# Tower elimination record
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class TowerStageRecord:
    """One stage of a sequential affine-hidden elimination tower.

    Records the block that was eliminated, its dimension, the action
    contribution from the Schur complement, and optional metadata for
    inspection and audit.
    """

    block_name: str
    """Human-readable name for the eliminated hidden block."""

    hidden_dim: int
    """Dimension of the eliminated hidden block at this stage."""

    half_log_det_increment: float
    """(1/2) log det(D_ee) — the fibre-volume action contribution from
    this stage's Schur complement.  This is a half-logdet, not a raw
    logdet, because the Gaussian integral contributes (1/2) log det to
    the effective action.  The full action shift for the stage also
    includes the coupling term -(1/2) j_e^T D_ee^{-1} j_e, which is
    folded into final_action but not stored separately here."""

    schur_condition_number: float | None = None
    """Condition number of the eliminated block's precision, if computed."""

    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class TowerEliminationRecord:
    """Complete record of a sequential affine-hidden elimination tower.

    Carries per-stage records, the accumulated half-logdet, and the
    final reduced datum ready for the regime detector.
    """

    stages: tuple[TowerStageRecord, ...]
    """Per-stage elimination records, in elimination order."""

    total_half_log_det: float
    """Accumulated (1/2) log det across all stages:
    sum of half_log_det_increment for each stage.

    This is the total fibre-volume contribution to the effective action,
    not a raw log-determinant.  To recover the raw log-determinant,
    multiply by 2.  To bridge to the evidence formula, note that the
    evidence dispatcher's log_quadratic_det_term also uses the (1/2)
    convention internally."""

    final_action: float
    """The visible action after all hidden blocks have been eliminated.
    Includes both fibre-volume (half-logdet) and coupling contributions."""

    final_coupling: Array
    """Residual coupling vector (zero-length if all hidden eliminated)."""

    final_precision: Array
    """Residual precision matrix (0×0 if all hidden eliminated)."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Typed evidence template (evidence.py output)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class RegularEvidenceTemplate:
    """Typed evidence template for regimes where the quadratic law applies.

    This is the output of the evidence dispatcher when the regime is
    regular (interior, cone, quotient, or quotient-cone).  It carries
    the log-evidence contribution decomposed into its structural pieces.

    The dispatcher does NOT produce a single scalar score.  It produces
    a typed template whose pieces the caller can inspect, combine with
    the log-likelihood, or compare across models.
    """

    regime: RegimeKind
    """Which regular regime this evidence was computed under."""

    active_dim: int
    """Dimension of the active space used in the evidence computation."""

    log_quadratic_det_term: float
    """(m/2) log(2π) - (1/2) log det(H_active).
    The ordinary Laplace determinant contribution."""

    log_cone_mass: float | None = None
    """log P(Z ∈ K) in the inverse-Hessian geometry.
    None means the cone mass has not been computed — either because
    the regime has no boundary (interior, quotient-without-boundary)
    or because computation is deferred (0.3.3 does not yet evaluate
    multivariate normal cone probabilities).
    For ABSTRACT cones, this will always be None."""

    cone_mass_known_exact: bool = True
    """Whether the cone mass formula is known to be exact (True) or
    would require approximation e.g. by Monte Carlo (False).

    This flag describes the *mathematical status* of the cone mass,
    not whether it has been computed.  When log_cone_mass is None,
    this flag indicates whether a future computation would be exact
    (True) or approximate (False).  When log_cone_mass is populated,
    True means it was computed exactly, False means it was estimated.

    For orthant and full-dimensional polyhedral cones the formula is
    known exactly (multivariate normal CDF).  For ABSTRACT cones the
    formula is also exact in principle but not extractable."""

    log_orbit_volume: float | None = None
    """log Vol(G · v̂).  None if no quotient was applied."""

    log_slice_jacobian: float | None = None
    """log of the slice-chart Jacobian.  None if no quotient was applied."""

    metadata: Mapping[str, Any] | None = None

    # ── Certified local asymptotic comparator API ────────────────────
    #
    # These methods assemble the full local evidence formula from the
    # stored n-independent structural pieces plus (n, S_vis) at query
    # time.  The general form is:
    #
    #   log Z_loc(n) = -n·S_vis - λ·log(n) + (m-1)·log(log(n)) + C
    #
    # For all regular templates: λ = active_dim/2, m = 1, and C is
    # the sum of the stored log-domain structural terms.
    #
    # NAMING DISCIPLINE: these are *local* quantities — the learning
    # exponent and evidence contribution of a single certified local
    # chart.  They become global only when localisation/dominance
    # conditions are verified externally.

    @property
    def local_learning_exponent(self) -> float:
        """Local learning exponent λ = active_dim / 2.

        This is the n-exponent in Z_loc ~ n^{-λ}: the rate at which
        the local Laplace integral concentrates with sample size.
        For regular regimes, λ = (number of active directions) / 2.
        """
        return self.active_dim / 2.0

    @property
    def local_multiplicity(self) -> int:
        """Local multiplicity m in Z_loc ~ n^{-λ} (log n)^{m-1}.

        For all regular regimes, m = 1 (no logarithmic correction).
        """
        return 1

    def log_local_evidence(self, n: int, S_vis: float) -> float:
        """Certified local asymptotic log-evidence at sample size n.

        Assembles the full formula:
            log Z_loc(n) = -n·S_vis - λ·log(n) + C

        where C = log_quadratic_det_term + log_cone_mass + log_orbit_volume
                + log_slice_jacobian (each omitted if None).

        Parameters
        ----------
        n : int
            Sample size (number of observations).
        S_vis : float
            Visible action evaluated at the local optimum (per sample).

        Returns
        -------
        float
            The log of the certified local evidence contribution.

        Raises
        ------
        ValueError
            If n < 1.
        """
        if n < 1:
            raise ValueError(f"Sample size must be ≥ 1, got {n}")
        import math
        lam = self.local_learning_exponent
        log_n = math.log(n)
        return -n * S_vis - lam * log_n + self._log_constant()

    def _log_constant(self) -> float:
        """The n-independent constant C in log Z_loc = -nS - λ log n + C."""
        c = self.log_quadratic_det_term
        if self.log_cone_mass is not None:
            c += self.log_cone_mass
        if self.log_orbit_volume is not None:
            c += self.log_orbit_volume
        if self.log_slice_jacobian is not None:
            c += self.log_slice_jacobian
        return c


@dataclass(frozen=True)
class UnresolvedEvidenceRefusal:
    """Honest refusal when the regime is outside quadratic control.

    This is not an error — it is a typed statement that the evidence law
    cannot be determined from quadratic data alone.  The kernel seed is
    carried forward so that a downstream kernel reduction (Layer 4) and
    singular classifier (Layer 5) can attempt further resolution.
    """

    kernel: QuadraticKernelSeed | None
    """The quadratic kernel seed from the detector.  None when the
    refusal comes from the singular classifier (Layer 5) rather than
    from the quadratic regime detector (Layer 3)."""

    kernel_dim: int
    """Dimension of the unresolved kernel.  Repeated here for
    quick inspection without unpacking the full datum."""

    reason: str = "nontrivial reduced kernel: quadratic data insufficient"
    """Human-readable explanation of why evidence was refused."""

    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class UnsupportedEvidenceRefusal:
    """Refusal when the detector returned an implementation-boundary case.

    The cone or slice requires a reduction step not yet implemented.
    This is distinct from UnresolvedEvidenceRefusal (mathematical
    obstruction).  The datum is carried for inspection.
    """

    cone_kind: ConeKind
    """Which cone kind triggered the refusal."""

    reason: str = ""
    """Human-readable explanation from the detector."""

    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class SingularEvidenceTemplate:
    """Evidence template for singular (non-quadratic) regimes.

    Produced when kernel reduction (Layer 4) and singular classification
    (Layer 5) succeed in certifying a singular channel.  Carries the
    structural invariants needed to compute the full evidence law.

    The full certified local asymptotic log-evidence is:

        log Z_loc(n) = -n·S_vis - λ·log(n) + (m-1)·log(log(n)) + C

    where:
      - WH: λ = p/2 + |ω|, m = 1, C = log_normal_prefactor + log_principal_integral
      - Critical branch: λ = p/2 + 1/2, m = 2, C = log_normal_prefactor + log_coefficient
      - General branch: λ and m not yet in closed form (structural data only)

    NAMING DISCIPLINE: log_normal_prefactor = (p/2)log(2π) - (1/2)log det(H_N)
    is the n-independent normal-direction constant.  The -p/2·log(n) factor
    from the normal Laplace integral is subsumed into λ (which includes p/2).
    """

    regime: RegimeKind
    """Which regime produced this template.  Should be one of the
    SINGULAR_* variants of RegimeKind."""

    active_dim: int
    """Total active dimension (kernel + normal)."""

    kernel_dim: int
    """Dimension of the unresolved kernel R."""

    positive_normal_dim: int
    """Dimension of the positive normal complement N."""

    log_normal_prefactor: float
    """Normal-direction Laplace prefactor: (dim_N/2)log(2π) - (1/2)log det(H_N).
    This is the n-independent constant only; the -p/2·log(n) factor is
    absorbed into local_learning_exponent."""

    singular_classification: Any = None
    """The full SingularClassification object from the Layer 5 classifier.
    Typed as Any to avoid circular import with singular_types."""

    log_quadratic_det_term: float | None = None
    """If the regular directions also contribute a Laplace determinant
    (e.g. from non-kernel positive-normal dims), this carries it."""

    metadata: Mapping[str, Any] | None = None

    # ── Certified local asymptotic comparator API ────────────────────

    @property
    def local_learning_exponent(self) -> float:
        """Local learning exponent λ.

        WH:              λ = positive_normal_dim/2 + total_weight
        Critical branch: λ = positive_normal_dim/2 + 1/2
        General branch:  raises ValueError (no closed form)
        """
        from .singular_types import (
            WeightedHomogeneousTemplate,
            CriticalBranchTemplate,
            BranchChannelTemplate,
        )
        sc = self.singular_classification
        p_half = self.positive_normal_dim / 2.0
        if isinstance(sc, WeightedHomogeneousTemplate):
            return p_half + sc.total_weight
        if isinstance(sc, CriticalBranchTemplate):
            return p_half + 0.5
        if isinstance(sc, BranchChannelTemplate):
            raise ValueError(
                "General branch channel: no certified closed-form "
                "learning exponent.  Structural data only."
            )
        raise TypeError(
            f"Unknown singular classification type: {type(sc)}"
        )

    @property
    def local_multiplicity(self) -> int:
        """Local multiplicity m in Z_loc ~ n^{-λ} (log n)^{m-1}.

        WH:              m = 1  (no logarithmic correction)
        Critical branch: m = 2  (one power of log n)
        General branch:  raises ValueError (no closed form)
        """
        from .singular_types import (
            WeightedHomogeneousTemplate,
            CriticalBranchTemplate,
            BranchChannelTemplate,
        )
        sc = self.singular_classification
        if isinstance(sc, WeightedHomogeneousTemplate):
            return 1
        if isinstance(sc, CriticalBranchTemplate):
            return 2
        if isinstance(sc, BranchChannelTemplate):
            raise ValueError(
                "General branch channel: no certified closed-form "
                "multiplicity.  Structural data only."
            )
        raise TypeError(
            f"Unknown singular classification type: {type(sc)}"
        )

    def log_local_evidence(self, n: int, S_vis: float) -> float | None:
        """Certified local asymptotic log-evidence at sample size n.

        Assembles the full formula:
            log Z_loc(n) = -n·S_vis - λ·log(n) + (m-1)·log(log(n)) + C

        Returns None if the n-independent constant C cannot be computed
        (e.g. principal integral not available for non-FULL_SPACE cones,
        or general branch channel with no closed-form law).

        Parameters
        ----------
        n : int
            Sample size (number of observations).
        S_vis : float
            Visible action evaluated at the local optimum (per sample).

        Returns
        -------
        float or None
            The log of the certified local evidence contribution,
            or None if the constant term is not fully determined.
        """
        if n < 1:
            raise ValueError(f"Sample size must be ≥ 1, got {n}")
        c = self._log_constant()
        if c is None:
            return None
        import math
        lam = self.local_learning_exponent
        m = self.local_multiplicity
        log_n = math.log(n)
        result = -n * S_vis - lam * log_n + c
        if m > 1:
            result += (m - 1) * math.log(log_n)
        return result

    def _log_constant(self) -> float | None:
        """The n-independent constant C in log Z_loc = -nS - λ log n + (m-1)log(log n) + C."""
        from .singular_types import (
            WeightedHomogeneousTemplate,
            CriticalBranchTemplate,
            BranchChannelTemplate,
        )
        sc = self.singular_classification
        if isinstance(sc, WeightedHomogeneousTemplate):
            if sc.log_principal_integral is None:
                return None
            return self.log_normal_prefactor + sc.log_principal_integral
        if isinstance(sc, CriticalBranchTemplate):
            return self.log_normal_prefactor + sc.log_coefficient
        if isinstance(sc, BranchChannelTemplate):
            return None  # No closed-form constant for general branch
        return None


# Tagged union for the evidence dispatcher output.
EvidenceResult = Union[
    RegularEvidenceTemplate,
    SingularEvidenceTemplate,
    UnresolvedEvidenceRefusal,
    UnsupportedEvidenceRefusal,
]


# ═══════════════════════════════════════════════════════════════════════
# Deprecation alias — remove in 0.3.4
# ═══════════════════════════════════════════════════════════════════════

ReducedKernelDatum = QuadraticKernelSeed
"""Deprecated alias.  Use QuadraticKernelSeed instead.

The old name was misleading: the detector's output is the quadratic-order
decomposition E_K = R ⊕ N (a seed for kernel reduction), not the paper's
"reduced kernel datum" (which is the post-reduction germ Φ on R).  The
rename aligns the code with the 0.3.3 Technical Note §9–10."""
