"""
Quadratic regime detector — 0.3.3 Technical Note §7–8.

The detector is the centrepiece of the typed local evidence geometry.
It takes a ReducedLocalDatum (the pre-classification object produced by
hidden elimination and any chart/slice reduction) and returns a
RegimeClassification (tagged union).

The classification is exact:
    R_K = {0}  →  one of the four regular regimes
    R_K ≠ {0}  →  UnresolvedKernel with a full QuadraticKernelSeed

The detector does NOT score evidence.  It classifies the regime and
provides the structural invariants needed by the evidence dispatcher
(evidence.py).  Keeping these layers separate is a design invariant.

0.3.3 scope
-----------
Supported cone kinds:

  - FULL_SPACE — no boundary, full active space.
  - ORTHANT — non-negativity constraints; cone spans the full space.
  - ABSTRACT — cone structure acknowledged but no numerical cone-mass
    computation; regime classification proceeds on the full active space.
  - POLYHEDRAL with full-dimensional cone — the cone {x : Ax >= 0}
    spans the full active space (detected via LP feasibility of strict
    inequalities, not rank(A)).

Refused (returns UnsupportedConeReduction):

  - POLYHEDRAL cones requiring genuine span reduction, i.e., cones
    whose linear hull is a proper subspace of the active space.

In all supported cases the cone span equals the full active space,
so the active datum IS the reduced datum.  The RegimeClassification
variants therefore carry the original ReducedLocalDatum unchanged.
When genuine polyhedral span reduction is implemented, the
classification will need to carry a genuinely reduced datum instead.

References
----------
0.3.3 Technical Note:
  - §7  Definition 8.1 (reduced quadratic datum)
  - §7  Theorem 8.1  (quadratic regime detector)
  - §8  Proposition 9.1 (kernel reduction)
"""
from __future__ import annotations

import numpy as np

from .exceptions import InputValidationError
from .regime_types import (
    ConeKind,
    ConeSpec,
    IndefiniteStationaryPoint,
    JacobianConvention,
    OrbitSpec,
    QuadraticKernelSeed,
    ReducedLocalDatum,
    RegimeClassification,
    RegimeKind,
    RegularCone,
    RegularInterior,
    RegularQuotient,
    RegularQuotientCone,
    UnresolvedKernel,
    UnsupportedConeReduction,
)
from .validation import Tolerances, resolve_tolerances


# ═══════════════════════════════════════════════════════════════════════
# Primary API: classify_regime
# ═══════════════════════════════════════════════════════════════════════

def classify_regime(
    datum: ReducedLocalDatum,
    tolerances: Tolerances | None = None,
) -> RegimeClassification:
    """Classify the local evidence regime of a reduced local datum.

    This is the main entry point for the quadratic regime detector.
    It implements Theorem 8.1 of the 0.3.3 Technical Note.

    Parameters
    ----------
    datum : ReducedLocalDatum
        The pre-classification reduced local object.  Its h_active may
        be PSD (not necessarily SPD) — that is the whole point.
    tolerances : Tolerances, optional
        Numerical tolerances for eigenvalue classification.

    Returns
    -------
    RegimeClassification
        One of six typed outcomes:

        - RegularInterior, RegularCone, RegularQuotient,
          RegularQuotientCone — quadratic control established.
        - UnresolvedKernel — mathematical obstruction; the reduced
          kernel is nontrivial and higher-order data is needed.
        - UnsupportedConeReduction — implementation boundary; the
          cone requires a span reduction not yet built in 0.3.3.

    Notes
    -----
    The detector does NOT compute evidence.  It returns the structural
    invariants needed by the evidence dispatcher.

    The classification boundary is sharp (Remark after Theorem 8.1):
    when the reduced kernel is nontrivial, germs with the same quadratic
    datum can have different leading exponents.
    """
    tol = resolve_tolerances(tolerances)
    _validate_datum(datum, tol)

    # Step 0: Three-way spectral pre-check.
    # Before the theorem-backed regime detector (which requires a PSD
    # Hessian), classify the eigenvalue spectrum into positive / zero /
    # negative bands.  A materially negative eigenvalue means the
    # stationary point is not a local minimum — refuse before the
    # detector sees it.
    indefinite = _check_indefinite(datum, tol)
    if indefinite is not None:
        return indefinite

    # Step 1: Compute the reduced quadratic datum.
    # Project H onto the cone span if cone ≠ full space.
    projection_result = _project_to_cone_span(datum, tol)
    if isinstance(projection_result, UnsupportedConeReduction):
        return projection_result
    h_k, active_dim_k = projection_result

    # Step 2: Eigendecompose H_K on the active span.
    eigenvalues, eigenvectors = np.linalg.eigh(h_k)
    rank_tol = _rank_tolerance(eigenvalues, tol)

    # Step 3: Classify the kernel R_K = ker(H_K).
    kernel_mask = eigenvalues <= rank_tol
    kernel_dim = int(np.sum(kernel_mask))

    if kernel_dim == 0:
        # R_K = {0} → regular regime.
        log_det = float(np.sum(np.log(eigenvalues)))
        return _classify_regular(datum, log_det)
    else:
        # R_K ≠ {0} → unresolved kernel.
        return _build_unresolved(
            datum, eigenvalues, eigenvectors, kernel_mask,
            kernel_dim, active_dim_k, tol,
        )


# ═══════════════════════════════════════════════════════════════════════
# Convenience constructor: from raw Hessian + optional cone/orbit
# ═══════════════════════════════════════════════════════════════════════

def classify_from_hessian(
    H: np.ndarray,
    cone: ConeSpec | None = None,
    orbit: OrbitSpec | None = None,
    tolerances: Tolerances | None = None,
) -> RegimeClassification:
    """Convenience wrapper: build a ReducedLocalDatum and classify.

    Parameters
    ----------
    H : array, shape (d, d)
        Symmetric PSD Hessian of the visible action at the optimum.
        For quotient cases (orbit is not None), this must be the
        **transverse Hessian** H_⊥ on the symmetry slice, not the
        full visible Hessian.  The dimension d should equal the
        transverse dimension (ambient_dim - orbit_dim).
    cone : ConeSpec, optional
        Admissible tangent cone.  Defaults to FULL_SPACE.
    orbit : OrbitSpec, optional
        Orbit/slice data.  Defaults to None (no quotient).
    tolerances : Tolerances, optional
        Numerical tolerances.

    Returns
    -------
    RegimeClassification
    """
    tol = resolve_tolerances(tolerances)
    h = np.asarray(H, dtype=float)
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise InputValidationError("H must be a square matrix")
    d = h.shape[0]

    if cone is None:
        cone = ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=d)

    datum = ReducedLocalDatum(
        active_dim=d,
        h_active=h,
        cone=cone,
        orbit=orbit,
    )
    return classify_regime(datum, tolerances=tol)


# ═══════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════

def _validate_datum(datum: ReducedLocalDatum, tol: Tolerances) -> None:
    """Basic structural validation of a ReducedLocalDatum.

    Checks shape, dimension consistency, and symmetry.  Does NOT check
    positive semidefiniteness — that is handled by the three-way spectral
    pre-check (_check_indefinite), which returns a typed result instead
    of raising an exception.
    """
    h = datum.h_active
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise InputValidationError("h_active must be a square matrix")
    if h.shape[0] != datum.active_dim:
        raise InputValidationError(
            f"h_active shape {h.shape} does not match active_dim={datum.active_dim}"
        )
    # Symmetry check (not SPD — that's the detector's job).
    if not np.allclose(h, h.T, atol=tol.atol, rtol=tol.rtol):
        raise InputValidationError("h_active must be symmetric")


def _check_indefinite(
    datum: ReducedLocalDatum,
    tol: Tolerances,
) -> IndefiniteStationaryPoint | None:
    """Three-way spectral pre-check: positive / zero / negative.

    Eigenvalues of h_active are classified against a single spectral
    tolerance ε = max(atol, rtol · max(1, ||H||_op)):

        λ > ε   → positive
        |λ| ≤ ε → zero (kernel candidate)
        λ < −ε  → negative (indefinite)

    If any eigenvalue is materially negative, returns an
    IndefiniteStationaryPoint refusal.  Otherwise returns None and
    the caller proceeds to the PSD-based regime detector.
    """
    eigenvalues = np.linalg.eigvalsh(datum.h_active)
    max_abs = float(np.max(np.abs(eigenvalues)))
    eps = max(tol.atol, tol.rtol * max(1.0, max_abs))

    negative_mask = eigenvalues < -eps
    negative_dim = int(np.sum(negative_mask))

    if negative_dim == 0:
        return None

    zero_mask = np.abs(eigenvalues) <= eps
    positive_mask = eigenvalues > eps

    min_eig = float(eigenvalues.min())
    spectral_gap = abs(min_eig) - eps

    # Suggest action based on severity
    if spectral_gap < 10 * eps:
        suggested = "polish"   # borderline — gradient polishing may fix it
    elif negative_dim <= 2:
        suggested = "refit"    # a few bad directions — try different init
    else:
        suggested = "discard"  # severely indefinite — not salvageable

    return IndefiniteStationaryPoint(
        datum=datum,
        negative_dim=negative_dim,
        zero_dim=int(np.sum(zero_mask)),
        positive_dim=int(np.sum(positive_mask)),
        min_eigenvalue=min_eig,
        spectral_gap=spectral_gap,
        eigenvalues=eigenvalues.copy(),
        suggested_action=suggested,
    )


def _project_to_cone_span(
    datum: ReducedLocalDatum,
    tol: Tolerances,
) -> tuple[np.ndarray, int] | UnsupportedConeReduction:
    """Project H_active onto the span of the cone if cone ≠ FULL_SPACE.

    For all supported cone kinds (FULL_SPACE, ORTHANT, ABSTRACT, and
    full-dimensional POLYHEDRAL), the cone spans the full active space,
    so H_active is returned unchanged.

    For POLYHEDRAL cones whose linear hull is a proper subspace of the
    active space, returns an UnsupportedConeReduction typed refusal
    (not an exception).  Full-dimensionality is tested via LP: the cone
    is full-dimensional iff it has nonempty interior, i.e. there exists
    x with Ax > 0 strictly.

    Returns
    -------
    (h_k, active_dim_k) : tuple
        The Hessian and active dimension (unchanged for supported cones), OR
    UnsupportedConeReduction
        Typed refusal for non-full-span polyhedral cones.
    """
    if datum.cone.kind == ConeKind.FULL_SPACE:
        return datum.h_active.copy(), datum.active_dim

    if datum.cone.kind == ConeKind.ABSTRACT:
        # Abstract cone: we cannot project numerically, so we operate
        # on the full active space.  The detector can still classify
        # kernel vs non-kernel; the cone mass just cannot be computed.
        return datum.h_active.copy(), datum.active_dim

    if datum.cone.kind == ConeKind.ORTHANT:
        # For orthant cones, the span is the full active space
        # (the cone {x : x_i >= 0} spans R^d).  The restriction
        # is to the same Hessian — the cone changes the integral
        # domain but not the span.
        return datum.h_active.copy(), datum.active_dim

    if datum.cone.kind == ConeKind.POLYHEDRAL:
        # For polyhedral cones K = {x : A x >= 0}, the span may be a
        # proper subspace of the active space.  The genuine
        # span-reduction construction (project H onto the cone span)
        # is deferred beyond 0.3.3.
        #
        # For 0.3.3, we refuse polyhedral cones whose span is not the
        # full active space.  This keeps the detector theorem-exact on
        # what it claims to handle.
        if datum.cone.halfspace_normals is None:
            raise InputValidationError(
                "POLYHEDRAL cone requires halfspace_normals"
            )
        A = datum.cone.halfspace_normals
        # Check whether K is full-dimensional (has nonempty interior).
        #
        # IMPORTANT: rank(A) does NOT determine this.  Example:
        # A = [[1, 0]] in R^2 gives K = {x : x_1 >= 0}, which is
        # full-dimensional despite rank(A) = 1 < 2.
        #
        # The correct test: K has full span iff it has nonempty
        # interior, i.e. there exists x with Ax > 0 (strictly).
        # We detect this by finding implicit equalities — rows i
        # where a_i^T x = 0 for ALL x in K.  If there are none,
        # the cone is full-dimensional.
        cone_span_dim = _polyhedral_cone_span_dim(A, datum.active_dim, tol)
        if cone_span_dim < datum.active_dim:
            return UnsupportedConeReduction(
                datum=datum,
                cone_kind=ConeKind.POLYHEDRAL,
                detected_span_dim=cone_span_dim,
                active_dim=datum.active_dim,
                reason=(
                    f"POLYHEDRAL cone has span of dimension "
                    f"{cone_span_dim} < active_dim={datum.active_dim}.  "
                    f"Non-full-span polyhedral cone reduction is not "
                    f"yet implemented in 0.3.3."
                ),
            )
        return datum.h_active.copy(), datum.active_dim

    raise InputValidationError(f"Unknown cone kind: {datum.cone.kind}")


def _rank_tolerance(eigenvalues: np.ndarray, tol: Tolerances) -> float:
    """Compute the rank tolerance for eigenvalue classification."""
    max_eig = float(np.max(np.abs(eigenvalues)))
    return max(tol.atol, tol.rtol * max(1.0, max_eig))


def _classify_regular(
    datum: ReducedLocalDatum,
    log_det: float,
) -> RegimeClassification:
    """Dispatch to the correct regular regime variant."""
    has_cone = datum.cone.kind != ConeKind.FULL_SPACE
    has_orbit = datum.orbit is not None

    if not has_cone and not has_orbit:
        return RegularInterior(datum=datum, log_det_h_active=log_det)
    elif has_cone and not has_orbit:
        return RegularCone(datum=datum, log_det_h_active=log_det)
    elif not has_cone and has_orbit:
        return RegularQuotient(
            datum=datum,
            log_det_h_active=log_det,
            orbit_dim=datum.orbit.orbit_dim,
        )
    else:
        return RegularQuotientCone(
            datum=datum,
            log_det_h_active=log_det,
            orbit_dim=datum.orbit.orbit_dim,
        )


def _build_unresolved(
    datum: ReducedLocalDatum,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    kernel_mask: np.ndarray,
    kernel_dim: int,
    active_dim_k: int,
    tol: Tolerances,
) -> UnresolvedKernel:
    """Build an UnresolvedKernel with the full QuadraticKernelSeed.

    Implements the decomposition E_K = R ⊕ N from Proposition 9.1:
    - R = ker(H_K) — the kernel directions
    - N = complement — the positive normal directions

    This is the quadratic-order seed, not yet the paper's "reduced
    kernel datum" (which requires the full action germ Φ).
    """
    # Kernel basis: eigenvectors corresponding to near-zero eigenvalues.
    kernel_basis = eigenvectors[:, kernel_mask].copy()

    # Positive normal eigenvalues and basis.
    normal_mask = ~kernel_mask
    normal_basis = eigenvectors[:, normal_mask].copy()
    normal_eigenvalues = eigenvalues[normal_mask].copy()
    positive_normal_dim = int(np.sum(normal_mask))

    # Normal prefactor: (dim_N / 2) log(2π) - (1/2) Σ log(λ_i)
    log_normal_prefactor = (
        0.5 * positive_normal_dim * np.log(2.0 * np.pi)
        - 0.5 * float(np.sum(np.log(normal_eigenvalues)))
    )

    # Residual kernel cone: restrict the original cone to the kernel.
    kernel_cone = _restrict_cone_to_kernel(datum.cone, kernel_basis, kernel_dim)

    kernel_datum = QuadraticKernelSeed(
        kernel_dim=kernel_dim,
        kernel_basis=kernel_basis,
        kernel_cone=kernel_cone,
        positive_normal_dim=positive_normal_dim,
        positive_normal_basis=normal_basis,
        positive_normal_eigenvalues=normal_eigenvalues,
        log_normal_prefactor=log_normal_prefactor,
    )

    return UnresolvedKernel(datum=datum, kernel=kernel_datum)


def _restrict_cone_to_kernel(
    cone: ConeSpec,
    kernel_basis: np.ndarray,
    kernel_dim: int,
) -> ConeSpec:
    """Restrict a cone specification to the kernel subspace."""
    if cone.kind == ConeKind.FULL_SPACE:
        return ConeSpec(kind=ConeKind.FULL_SPACE, ambient_dim=kernel_dim)

    if cone.kind == ConeKind.ABSTRACT:
        return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=kernel_dim)

    if cone.kind == ConeKind.ORTHANT:
        # The intersection of the positive orthant with a subspace is
        # in general a polyhedral cone in the kernel coordinates.
        # For now, mark as ABSTRACT — exact orthant restriction is a
        # future numerical extension.
        return ConeSpec(
            kind=ConeKind.ABSTRACT,
            ambient_dim=kernel_dim,
            metadata={"derived_from": "orthant_kernel_restriction"},
        )

    if cone.kind == ConeKind.POLYHEDRAL:
        if cone.halfspace_normals is not None:
            # Project the half-space normals onto the kernel.
            # A_kernel = A @ kernel_basis, giving constraints in kernel coords.
            A_kernel = cone.halfspace_normals @ kernel_basis
            return ConeSpec(
                kind=ConeKind.POLYHEDRAL,
                ambient_dim=kernel_dim,
                halfspace_normals=A_kernel,
                metadata={"derived_from": "polyhedral_kernel_restriction"},
            )
        return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=kernel_dim)

    return ConeSpec(kind=ConeKind.ABSTRACT, ambient_dim=kernel_dim)


def _polyhedral_cone_span_dim(
    A: np.ndarray,
    ambient_dim: int,
    tol: Tolerances,
) -> int:
    """Compute the dimension of span(K) for K = {x : Ax >= 0}.

    The span of a homogeneous polyhedral cone is the ambient space
    minus the implicit equalities.  Row i of A is an implicit equality
    if a_i^T x = 0 for ALL x in K, i.e., there is no x in K with
    a_i^T x > 0.

    We detect implicit equalities via LP: for each row i, solve

        max  a_i^T x
        s.t. Ax >= 0,  sum_j (Ax)_j <= 1

    The bounding constraint sum_j (Ax)_j <= 1 keeps the LP bounded
    (since K is a cone, without it the objective is either 0 or +inf).

    If the optimal value is 0 (within tolerance), row i is an implicit
    equality.  The span dimension is then:

        ambient_dim - rank(A_eq)

    where A_eq is the matrix of implicit-equality rows.

    IMPORTANT: rank(A) alone does NOT determine this.  The cone
    {x in R^d : a^T x >= 0} for a single nonzero a has rank(A) = 1
    but is always full-dimensional (span = R^d).
    """
    from scipy.optimize import linprog

    m, d = A.shape
    assert d == ambient_dim

    # Constraints: -A x <= 0  (i.e. Ax >= 0)
    #              1^T A x <= 1  (bounding — prevents unbounded cone)
    neg_A = -A
    ones_A = np.sum(A, axis=0, keepdims=True)  # shape (1, d)
    A_ub = np.vstack([neg_A, ones_A])  # shape (m+1, d)
    b_ub = np.zeros(m + 1)
    b_ub[-1] = 1.0  # sum(Ax) <= 1

    implicit_rows = []
    cutoff = max(tol.atol, tol.rtol)

    for i in range(m):
        c = -A[i]  # minimize -a_i^T x  <==>  maximize a_i^T x
        result = linprog(
            c, A_ub=A_ub, b_ub=b_ub,
            bounds=(None, None),  # x is free (unbounded below)
            method="highs",
        )
        if result.success:
            opt_val = -result.fun  # max a_i^T x
            if opt_val <= cutoff:
                implicit_rows.append(i)
        else:
            # If LP is infeasible, the cone is {0} — all rows are
            # implicit equalities.  But this shouldn't happen for a
            # well-formed cone (x=0 is always feasible).  If it does,
            # treat the row as implicit equality conservatively.
            implicit_rows.append(i)

    if not implicit_rows:
        return ambient_dim

    A_eq = A[implicit_rows]
    eq_rank = int(np.linalg.matrix_rank(
        A_eq,
        tol=max(tol.atol, tol.rtol * max(1.0, float(np.max(np.abs(A_eq)))))
    ))
    return ambient_dim - eq_rank
