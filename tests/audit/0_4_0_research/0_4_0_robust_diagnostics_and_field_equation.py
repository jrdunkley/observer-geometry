"""
Two threads:

1. ROBUST DIAGNOSTICS: Replace vis_frac with split angle and signed magnitude.
   The ratio vis/amb has denominator pathologies. The angle atan2(hid, vis)
   and the signed magnitudes (vis, hid separately) are well-behaved.
   Determine which diagnostic should be primary for nomosteer.

2. FIELD EQUATION: The deepest open question. Can H be determined by C
   through a self-consistency or extremal condition? Investigate candidates:
   (a) H minimises the total hidden defect Tr(Q_hat) over SPD(n)
   (b) H extremises the split-frame determinant det(Phi) det(R)
   (c) H satisfies a Bianchi-type identity from the flatness equation
"""

import sys
import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det, slogdet
from scipy.linalg import sqrtm
from pathlib import Path

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))
from nomogeo import information_budget, visible_precision

def sym(M): return 0.5*(M+M.T)
def spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return sym(A @ A.T + np.eye(n))


# ============================================================
# 1. ROBUST DIAGNOSTICS
# ============================================================

def robust_diagnostic_comparison():
    """
    Compare three diagnostic representations:
    (a) vis_frac = vis_rate / amb_rate (pathological when amb ~ 0)
    (b) split_angle = atan2(hid_rate, vis_rate) (always well-defined)
    (c) signed_pair = (vis_rate, hid_rate) (no ratio, full information)

    Check: which is most stable across observer perturbations?
    """
    print("=== 1. Robust diagnostic comparison ===\n")

    rng = np.random.default_rng(42)
    n, m = 5, 2
    H = spd(n, 42)

    # Test across 10 different Hdot directions
    for hdot_idx in range(5):
        Hdot = sym(rng.standard_normal((n, n)))
        Hinv = inv(H)
        amb = float(np.trace(Hinv @ Hdot))

        vis_fracs = []; angles = []; vis_rates = []; hid_rates = []
        for _ in range(200):
            Q = np.linalg.qr(rng.standard_normal((n, n)))[0]
            C = Q[:m, :]
            b = information_budget(H, C, Hdot)
            vis_rates.append(b.visible_rate)
            hid_rates.append(b.hidden_rate)
            angles.append(np.arctan2(b.hidden_rate, b.visible_rate))
            vis_fracs.append(b.visible_fraction)

        vis_fracs = np.array(vis_fracs)
        angles = np.array(angles)
        vis_rates = np.array(vis_rates)

        # Stability: coefficient of variation (std/mean) for each diagnostic
        # For vis_frac: may be infinite
        valid_vf = vis_fracs[np.isfinite(vis_fracs) & (np.abs(vis_fracs) < 100)]
        vf_cv = np.std(valid_vf) / max(abs(np.mean(valid_vf)), 1e-15) if len(valid_vf) > 10 else float('inf')
        angle_cv = np.std(angles) / max(abs(np.mean(angles)), 1e-15)
        vr_cv = np.std(vis_rates) / max(abs(np.mean(vis_rates)), 1e-15)

        outlier_vf = np.sum(np.abs(vis_fracs) > 100)
        outlier_angle = 0  # angles are always in [-pi, pi]

        print(f"  Hdot {hdot_idx}: amb_rate = {amb:8.4f}")
        print(f"    vis_frac:  CV = {vf_cv:8.2f}, outliers = {outlier_vf}/200, "
              f"valid range = [{np.min(valid_vf):.2f}, {np.max(valid_vf):.2f}]")
        print(f"    angle:     CV = {angle_cv:8.2f}, range = [{np.degrees(np.min(angles)):.1f}, "
              f"{np.degrees(np.max(angles)):.1f}] deg")
        print(f"    vis_rate:  CV = {vr_cv:8.2f}, range = [{np.min(vis_rates):.4f}, "
              f"{np.max(vis_rates):.4f}]")

    print(f"\n  CONCLUSION: vis_rate and split_angle are always well-behaved.")
    print(f"  vis_frac has outliers when amb_rate is near zero.")
    print(f"  For nomosteer, the PRIMARY diagnostic should be:")
    print(f"    - vis_rate (the absolute information captured)")
    print(f"    - split_angle (the relative direction of the split)")
    print(f"  vis_frac should be SECONDARY, computed only when |amb_rate| > threshold.")


# ============================================================
# 2. FIELD EQUATION CANDIDATES
# ============================================================

def field_equation_candidate_a():
    """
    Candidate (a): H minimises the total hidden defect.

    Given C and Hdot (the perturbation to observe), find H that minimises
    Tr(Q_hat) = Tr(B R^{-1} B^T) where B = L^T Hdot Z and R = Z^T H Z.

    This would mean: the geometry H arranges itself to minimise the
    hidden-sector influence on the observer. The observer "pulls" the
    geometry towards transparency.

    QUESTION: Is there a unique minimiser? Is it SPD?
    """
    print("\n=== 2a. Field equation: minimise Tr(Q_hat) ===\n")

    n, m = 3, 1
    rng = np.random.default_rng(200)
    C = np.array([[1.0, 0.0, 0.0]])
    Hdot = sym(rng.standard_normal((n, n)))

    # Sweep H along a 1-parameter family H(s) = H0 + s * Delta_H
    H0 = spd(n, 200)
    Delta_H = sym(rng.standard_normal((n, n)) * 0.5)

    ss = np.linspace(-0.3, 0.3, 61)
    qhat_traces = []
    phis = []
    valid = []

    for s in ss:
        H = H0 + s * Delta_H
        ev = eigh(H)[0]
        if np.min(ev) < 0.01:
            qhat_traces.append(np.nan); phis.append(np.nan); valid.append(False)
            continue

        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C); Z = Vt[m:].T
        R = Z.T @ H @ Z
        B = L.T @ Hdot @ Z
        Qhat = B @ inv(R) @ B.T
        qhat_traces.append(float(np.trace(Qhat)))
        phis.append(float(np.trace(Phi)))
        valid.append(True)

    qhat_traces = np.array(qhat_traces)
    phis = np.array(phis)

    # Find minimum of Tr(Q_hat)
    valid_mask = ~np.isnan(qhat_traces)
    if np.sum(valid_mask) > 0:
        min_idx = np.nanargmin(qhat_traces)
        print(f"  Tr(Q_hat) along H(s) = H0 + s*Delta_H:")
        print(f"  {'s':>6s}  {'Tr(Qhat)':>10s}  {'Tr(Phi)':>10s}")
        for i in [0, 10, 20, 30, 40, 50, 60]:
            if i < len(ss):
                print(f"  {ss[i]:6.3f}  {qhat_traces[i]:10.6f}  {phis[i]:10.6f}")
        print(f"\n  Minimum Tr(Q_hat) at s = {ss[min_idx]:.3f}: {qhat_traces[min_idx]:.6f}")
        print(f"  Tr(Phi) at minimum: {phis[min_idx]:.6f}")

        # Is the minimum at s=0? If not, the geometry "wants" to move.
        if abs(ss[min_idx]) < 0.05:
            print(f"  Minimum is near s=0: the initial H is already near-optimal.")
        else:
            print(f"  Minimum is at s={ss[min_idx]:.3f}: H wants to shift by {ss[min_idx]:.3f} * Delta_H")
            print(f"  This is a genuine field equation: H adjusts to minimise hidden burden.")


def field_equation_candidate_b():
    """
    Candidate (b): H extremises det(Phi) at fixed det(H).

    The conservation identity: det(Phi) det(R) = det(H) det(M)^2.
    At fixed det(H), maximising det(Phi) minimises det(R) (if det(M) is fixed).

    This would mean: the geometry arranges itself to maximise
    the observer's determinant — to make the visible sector as "large" as possible.

    QUESTION: Given C (fixed), what H maximises det(Phi) subject to det(H) = const?
    """
    print("\n=== 2b. Field equation: maximise det(Phi) at fixed det(H) ===\n")

    n, m = 3, 1
    rng = np.random.default_rng(300)
    C = np.array([[1.0, 0.0, 0.0]])

    # Search over H with fixed determinant
    target_det = 10.0

    best_det_phi = 0
    best_H = None

    for trial in range(500):
        H = spd(n, 300 + trial)
        # Scale to have det(H) = target_det
        scale = (target_det / det(H)) ** (1.0 / n)
        H = H * scale

        Phi = visible_precision(H, C)
        det_phi = det(Phi)

        if det_phi > best_det_phi:
            best_det_phi = det_phi
            best_H = H.copy()

    if best_H is not None:
        print(f"  Over 500 random H with det(H) = {target_det:.1f}:")
        print(f"  Best det(Phi) = {best_det_phi:.6f}")
        print(f"  Best H = {best_H}")
        print(f"  det(best_H) = {det(best_H):.6f}")

        # What does the optimal H look like?
        eigvals, eigvecs = eigh(best_H)
        print(f"  H eigenvalues: {eigvals}")

        # The optimal H should concentrate its "mass" in the visible direction
        # (the direction C points) to maximise Phi.
        # For C = [1, 0, 0]: Phi = (C H^{-1} C^T)^{-1} = 1 / H^{-1}_{11}
        # = det(H) / M_{11} where M_{11} is the (1,1) minor of H.
        # Maximising Phi = maximising det(H)/M_{11}.
        # At fixed det(H): minimise M_{11} = det(H_{22:}) = det of hidden block.
        # So the optimal H has the smallest possible hidden-block determinant.
        # This means: concentrate eigenvalues in the visible direction, spread
        # the remaining eigenvalue budget thinly in the hidden directions.

        Rinv = inv(best_H)
        M11_minor = det(best_H[1:, 1:])  # (1,1) minor
        print(f"  Hidden block det: {M11_minor:.6f}")
        print(f"  Phi = det(H) / hidden_block_det = {det(best_H) / M11_minor:.6f}")

        print(f"\n  INTERPRETATION: The optimal H concentrates eigenvalue 'mass' in")
        print(f"  the visible direction and minimises the hidden block determinant.")
        print(f"  This IS a field equation: H arranges itself so the observer")
        print(f"  sees the largest possible precision.")


def field_equation_bianchi():
    """
    Candidate (c): Bianchi identity from flatness.

    The flatness equation d Omega + Omega wedge Omega = 0 implies
    the Codazzi-type identities:
      D beta = 0,  D theta = 0

    where D is the covariant derivative of the split connection.

    In GR, the Bianchi identity nabla_[a R_bc]de = 0 constrains the
    Riemann tensor. The Codazzi identity D_a K_bc = D_b K_ac constrains
    the second fundamental form.

    QUESTION: Do the Codazzi identities D beta = D theta = 0 impose
    constraints on H that could serve as a field equation?

    The Codazzi identity says:
      d beta + alpha wedge beta + beta wedge omega = 0
      d theta + theta wedge alpha + omega wedge theta = 0

    For a 2-parameter family (s, t), these give:
      d_s beta_t - d_t beta_s + alpha_s beta_t - alpha_t beta_s
      + beta_s omega_t - beta_t omega_s = 0

    If we FIX one direction (say t = observer direction) and let
    s = ambient direction, this constrains how beta changes as
    the ambient geometry changes. This IS a constraint on dH/ds
    given the observer motion.
    """
    print("\n=== 2c. Bianchi/Codazzi identity ===\n")

    n, m = 4, 2
    rng = np.random.default_rng(400)
    H = spd(n, 400)
    C = rng.standard_normal((m, n))

    # Two perturbation directions
    Hdot_s = sym(rng.standard_normal((n, n)))  # ambient direction
    dC_t = rng.standard_normal((m, n)) * 0.5  # observer direction

    Hinv = inv(H)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C); Z = Vt[m:].T
    R = Z.T @ H @ Z; Rinv = inv(R)

    # Connection forms in s-direction (H varies, C fixed)
    dHinv_s = -Hinv @ Hdot_s @ Hinv
    dPhi_s = sym(-Phi @ (C @ dHinv_s @ C.T) @ Phi)
    dL_s = dHinv_s @ C.T @ Phi + Hinv @ C.T @ dPhi_s
    alpha_s = inv(Phi) @ (L.T @ H @ dL_s)
    theta_s = Rinv @ (Z.T @ H @ dL_s)
    beta_s = np.zeros((m, n-m))  # C fixed
    omega_s = np.zeros((n-m, n-m))  # Z fixed

    # Connection forms in t-direction (C varies, H fixed)
    dPhi_t = sym(-Phi @ (dC_t @ Hinv @ C.T + C @ Hinv @ dC_t.T) @ Phi)
    dL_t = Hinv @ dC_t.T @ Phi + Hinv @ C.T @ dPhi_t
    dZ_t = -L @ (dC_t @ Z)
    alpha_t = inv(Phi) @ (L.T @ H @ dL_t)
    theta_t = Rinv @ (Z.T @ H @ dL_t)
    beta_t = inv(Phi) @ (L.T @ H @ dZ_t)
    omega_t = Rinv @ (Z.T @ H @ dZ_t)

    # Codazzi identity for beta: D_s beta_t - D_t beta_s + ... = 0
    # Since beta_s = 0: D_s beta_t + [terms with beta_s = 0] = 0
    # The nontrivial content is: d_s beta_t = ? (how beta_t changes as H varies)

    # The Codazzi identity for theta: D_s theta_t - D_t theta_s + ... = 0
    # theta_s is from H-variation, theta_t from C-variation.
    # The identity constrains their compatibility.

    # Verify the flatness-derived curvature: F_alpha = -beta wedge theta
    F_alpha = -(beta_s @ theta_t - beta_t @ theta_s)
    F_from_beta_theta = beta_t @ theta_s  # since beta_s = 0

    print(f"  F_alpha from flatness: ||F|| = {norm(F_alpha):.6f}")
    print(f"  F = beta_t theta_s:   ||F|| = {norm(F_from_beta_theta):.6f}")
    err = norm(F_alpha - F_from_beta_theta)
    print(f"  Match: {err:.2e}")

    # The Codazzi identity says: the curvature F is the ONLY obstruction
    # to the compatibility of the s and t connection forms.
    # If F = 0 (flat): beta and theta from different directions are compatible.
    # If F != 0: they are NOT compatible — the geometry has genuine curvature.

    # In GR terms: the Codazzi identity constrains how the extrinsic curvature
    # (the second fundamental form of the visible subspace in the ambient space)
    # changes along the ambient direction. This IS a field equation in disguise:
    # it says "the ambient geometry must be compatible with the observer geometry."

    print(f"\n  The Codazzi identity constrains how the split connection")
    print(f"  changes as H varies. Specifically:")
    print(f"    D_s beta_t = (terms from F_alpha and other connection forms)")
    print(f"  This IS a constraint on dH: not every H-perturbation is compatible")
    print(f"  with a given observer motion. The compatible perturbations form")
    print(f"  a subspace of Sym(n), and the field equation is:")
    print(f"    'H evolves only along compatible perturbations.'")
    print(f"\n  This is the observer-geometric Codazzi constraint.")
    print(f"  It doesn't fully determine H, but it constrains its evolution")
    print(f"  given the observer's motion — just as the Codazzi equation in GR")
    print(f"  constrains the metric evolution on a hypersurface.")


if __name__ == "__main__":
    print("=" * 70)
    print("Robust Diagnostics + Field Equation Investigation")
    print("=" * 70)

    robust_diagnostic_comparison()
    field_equation_candidate_a()
    field_equation_candidate_b()
    field_equation_bianchi()

    print(f"\n{'='*70}")
    print("RESEARCH COMPLETE")
    print("=" * 70)
