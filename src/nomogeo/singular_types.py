"""
Typed singular evidence templates — 0.3.3 Technical Note §10–12.

These types carry the output of the certified singular classifier (Layer 5).
Two exact channels are supported:

    1. Weighted-homogeneous principal germs (Proposition 10.1)
       — the kernel action Φ admits a positive principal germ P under
         anisotropic dilations, yielding learning exponent dim_N/2 + |ω|.

    2. Branch-restiffening channels (Theorems 11.1–11.3)
       — a moving one-dimensional fibre minimum with gradually returning
         transverse stiffness.  The critical quadratic branch channel
         yields an explicit logarithmic correction.

When neither channel is certified, the classifier returns an honest
refusal (UnresolvedKernelJetRefusal).

References
----------
0.3.3 Technical Note:
  - §10  Proposition 10.1 (anisotropic kernel law)
  - §11  Theorems 11.1–11.3 (branch-restiffening channels)
  - §12  Definition 12.1 (reduced kernel datum)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Union

import numpy as np
from numpy.typing import NDArray

from .kernel_reduction import PolynomialJet, ReducedKernelActionDatum

Array = NDArray[np.float64]


# ═══════════════════════════════════════════════════════════════════════
# Singular regime kinds
# ═══════════════════════════════════════════════════════════════════════

class SingularRegimeKind(str, Enum):
    """Typed singular regime classification from the kernel classifier."""

    WEIGHTED_HOMOGENEOUS = "weighted_homogeneous"
    """The kernel germ admits a positive weighted-homogeneous principal
    part.  The evidence law is governed by the anisotropic kernel
    integral with exponent |ω| = ω₁ + ... + ωᵣ."""

    BRANCH_RESTIFFENING = "branch_restiffening"
    """The kernel germ admits an analytic branch channel datum with
    a certified branch-restiffening structure.  General (non-critical)
    branch channels."""

    CRITICAL_QUADRATIC_BRANCH = "critical_quadratic_branch"
    """The kernel germ is in the critical quadratic branch channel
    (c=2, quadratic stiffness vanishing).  The evidence law has an
    explicit logarithmic correction: n^{-1/2} κ log n."""

    UNRESOLVED_KERNEL_JET = "unresolved_kernel_jet"
    """The kernel jet has not been resolved into any certified
    singular channel.  Honest refusal."""


# ═══════════════════════════════════════════════════════════════════════
# Weighted-homogeneous singular template
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class WeightedHomogeneousTemplate:
    """Evidence template for the weighted-homogeneous kernel sector.

    Proposition 10.1: if the reduced kernel action Φ admits a positive
    weighted-homogeneous principal germ P of degree 1 with weights
    (ω₁, ..., ωᵣ), then the kernel integral contributes n^{-|ω|} times
    the principal kernel integral ∫_{K_R} e^{-P(s)} ds.

    The full local evidence exponent is dim_N/2 + |ω|.
    """

    kind: SingularRegimeKind = field(
        default=SingularRegimeKind.WEIGHTED_HOMOGENEOUS, init=False
    )

    kernel_dim: int
    """Dimension of the kernel space R."""

    weights: tuple[float, ...]
    """Anisotropic weights (ω₁, ..., ωᵣ).  Each ωᵢ > 0.
    The learning exponent from the kernel is |ω| = Σ ωᵢ."""

    total_weight: float
    """|ω| = ω₁ + ... + ωᵣ.  The kernel learning exponent."""

    log_principal_integral: float | None = None
    """log ∫_{K_R} e^{-P(s)} ds.  None if not yet computed
    (e.g. because cone integration is deferred)."""

    log_normal_prefactor: float = 0.0
    """Normal-direction prefactor from kernel reduction."""

    positive_normal_dim: int = 0
    """Dimension of the positive normal complement N."""

    leading_monomial_degrees: tuple[int, ...] | None = None
    """The degrees (d₁, ..., dᵣ) such that ωᵢ = 1/dᵢ.  Present
    when the principal germ is a sum of pure powers Σ cᵢ uᵢ^{dᵢ}."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Branch-restiffening singular template
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class BranchChannelDatum:
    """Data structure for an analytic branch channel.

    Records the branch straightening parameters from Definition 11.1:
    the branch order c, stiffness vanishing order, and intrinsic
    branch-channel exponent β_br = ord₀(λ) / c.
    """

    branch_order: int
    """c ≥ 2: the order of the first non-vanishing transverse derivative
    along the straightened branch.  c=2 is the quadratic branch case."""

    stiffness_vanishing_order: int
    """ord₀(λ): the vanishing order of the branch restoring stiffness
    λ(v) = (1/c!) ∂ᵤᶜ Φ̃(0,v) as v → 0.  This is intrinsic by
    Theorem 11.1."""

    branch_channel_exponent: float
    """β_br = ord₀(λ) / c.  The intrinsic branch-channel exponent."""

    branch_action_degree: int | None = None
    """d: degree of the leading branch action term B₀ vᵈ along
    the minimum curve Φ̃(0,v) = B₀ vᵈ + O(vᵈ⁺¹).  Present for
    the critical quadratic branch channel."""

    higher_transverse_degree: int | None = None
    """a: degree of the higher-order transverse term C₀ uᵃ at the
    origin.  Present for the critical quadratic branch channel."""

    kappa: float | None = None
    """κ = (d(a-2) - 2a) / (2a).  The critical-balance parameter.
    Present (and > 0) only in the critical quadratic branch channel."""

    stiffness_leading_coeff: float | None = None
    """H₀ > 0: leading coefficient of the quadratic stiffness
    (1/2)∂²ᵤ Φ̃(0,v) = H₀ v² + O(v³).  For critical quadratic branch."""

    branch_leading_coeff: float | None = None
    """B₀ > 0: leading coefficient of the branch action."""

    higher_transverse_coeff: float | None = None
    """C₀ > 0: leading coefficient of the higher transverse term."""

    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class BranchChannelTemplate:
    """Evidence template for a general branch-restiffening channel.

    For non-critical branch channels, the exact asymptotic law depends
    on the specific balance of terms.  This template records the
    structural data but may not carry an explicit closed-form law.
    """

    kind: SingularRegimeKind = field(
        default=SingularRegimeKind.BRANCH_RESTIFFENING, init=False
    )

    kernel_dim: int
    """Dimension of the kernel space R (= 2 for branch channels)."""

    branch_datum: BranchChannelDatum
    """The certified branch channel structural data."""

    log_normal_prefactor: float = 0.0
    """Normal-direction prefactor from kernel reduction."""

    positive_normal_dim: int = 0
    """Dimension of the positive normal complement N."""

    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class CriticalBranchTemplate:
    """Evidence template for the critical quadratic branch channel.

    Theorem 11.3: the reduced kernel integral obeys

        ∫ e^{-nΦ(z)} b(z) dz = n^{-1/2} [b(0) (√π/d) H₀^{-1/2} κ log n + O(1)]

    so the full local evidence contribution is

        n^{-dim_N/2 - 1/2} · det(H_N)^{-1/2} · (√π/d) H₀^{-1/2} κ log n

    The logarithmic correction is the hallmark of this sector.
    """

    kind: SingularRegimeKind = field(
        default=SingularRegimeKind.CRITICAL_QUADRATIC_BRANCH, init=False
    )

    kernel_dim: int
    """Dimension of the kernel space R (= 2 for branch channels)."""

    branch_datum: BranchChannelDatum
    """The certified branch channel structural data."""

    log_coefficient: float
    """log of the coefficient in front of log(n):
    log(b(0) · (√π/d) · H₀^{-1/2} · κ).
    This coefficient multiplies log(n) in the asymptotic law."""

    kappa: float
    """κ = (d(a-2) - 2a) / (2a).  The critical-balance parameter."""

    log_normal_prefactor: float = 0.0
    """Normal-direction prefactor from kernel reduction."""

    positive_normal_dim: int = 0
    """Dimension of the positive normal complement N."""

    metadata: Mapping[str, Any] | None = None


# ═══════════════════════════════════════════════════════════════════════
# Unresolved kernel jet refusal
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class UnresolvedKernelJetRefusal:
    """Honest refusal when the kernel jet cannot be certified.

    The kernel action germ does not fall into any of the supported
    certified singular channels.  Possible reasons:

    - The principal jet is not positive on the active kernel cone.
    - Multiple competing faces have different dominant sectors.
    - The branch chart is not uniquely certifiable.
    - The kernel dimension is > 2 (branch channels require dim = 2).
    - No stable dominant face has been identified.

    The package never silently guesses a singular law.
    """

    kind: SingularRegimeKind = field(
        default=SingularRegimeKind.UNRESOLVED_KERNEL_JET, init=False
    )

    kernel_dim: int
    """Dimension of the unresolved kernel."""

    reason: str = "kernel jet not resolved into any certified singular channel"
    """Human-readable explanation of the refusal."""

    log_normal_prefactor: float = 0.0
    """Normal-direction prefactor (still valid from kernel reduction)."""

    positive_normal_dim: int = 0
    """Dimension of the positive normal complement."""

    metadata: Mapping[str, Any] | None = None


# Tagged union for the singular classifier output.
SingularClassification = Union[
    WeightedHomogeneousTemplate,
    BranchChannelTemplate,
    CriticalBranchTemplate,
    UnresolvedKernelJetRefusal,
]
