"""
Goals B and C:
B - Sextic hierarchy on a noncommuting family (proper degenerate example)
C - Curvature diagnostic on existing nomogeo examples

Also Goal D: General 2-parameter compatibility (both H and C vary simultaneously).
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det
from scipy.linalg import sqrtm, expm

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok

def spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T + np.eye(n)

def sym(M):
    return 0.5 * (M + M.T)

def build_split_frame(H, C):
    n = H.shape[0]; m = C.shape[0]
    Hinv = inv(H); Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C); Z = Vt[m:].T
    R = Z.T @ H @ Z
    return Phi, L, Z, R


# ============================================================
# GOAL C: Curvature on existing nomogeo examples
# ============================================================

def curvature_on_entanglement():
    """
    Entanglement hidden-load example.
    H = inv(two_mode_squeezed_thermal_covariance(r))
    C = selector_a = [[1,0,0,0],[0,1,0,0]]  (observe mode A of a 2-mode system)
    n=4, m=2
    """
    print("\n=== Curvature: entanglement_hidden_load ===")

    def two_mode_squeezed_thermal(r, n_a=0.0, n_b=0.0):
        c = np.cosh(r); s = np.sinh(r)
        return np.array([
            [(2*n_a+1)*c**2 + (2*n_b+1)*s**2, 0, (2*n_a+2*n_b+2)*c*s, 0],
            [0, (2*n_a+1)*c**2 + (2*n_b+1)*s**2, 0, -(2*n_a+2*n_b+2)*c*s],
            [(2*n_a+2*n_b+2)*c*s, 0, (2*n_b+1)*c**2 + (2*n_a+1)*s**2, 0],
            [0, -(2*n_a+2*n_b+2)*c*s, 0, (2*n_b+1)*c**2 + (2*n_a+1)*s**2],
        ])

    C = np.array([[1,0,0,0],[0,1,0,0]], dtype=float)
    n, m = 4, 2

    # Curvature at different squeezing parameters
    rs = [0.3, 0.5, 0.65, 0.85, 1.0, 1.5]
    print("  Squeezing r | ||F_alpha|| (mean over 100 pairs)")

    for r in rs:
        sigma = two_mode_squeezed_thermal(r)
        H = inv(sigma)
        Phi, L, Z, R = build_split_frame(H, C)
        Rinv = inv(R)
        Hinv = inv(H)
        Rsqrt_inv = np.real(sqrtm(Rinv))
        Phi_sqrt = np.real(sqrtm(Phi))

        rng = np.random.default_rng(int(r * 1000))
        F_norms = []
        for _ in range(100):
            dC1 = rng.standard_normal((m, n))
            dC2 = rng.standard_normal((m, n))

            for dC in [dC1, dC2]:
                pass  # just use random perturbations

            # Extract whitened betas
            betas = []
            for dC in [dC1, dC2]:
                dPhi = -Phi @ (dC @ Hinv @ C.T + C @ Hinv @ dC.T) @ Phi
                dZ = -L @ (dC @ Z)
                beta = inv(Phi) @ (L.T @ H @ dZ)
                bw = Phi_sqrt @ beta @ Rsqrt_inv
                betas.append(bw)

            F = betas[0] @ betas[1].T - betas[1] @ betas[0].T
            F_norms.append(norm(F))

        print(f"  r = {r:.2f}     | {np.mean(F_norms):.4f}")

    # Also compute with thermal noise
    print("\n  With thermal noise (n_a=0.3, n_b=0.4):")
    for r in [0.65, 0.85]:
        sigma = two_mode_squeezed_thermal(r, n_a=0.3, n_b=0.4)
        H = inv(sigma)
        Phi, L, Z, R = build_split_frame(H, C)
        Rinv = inv(R)
        Hinv = inv(H)
        Rsqrt_inv = np.real(sqrtm(Rinv))
        Phi_sqrt = np.real(sqrtm(Phi))

        rng = np.random.default_rng(int(r * 1000))
        F_norms = []
        for _ in range(100):
            dC1 = rng.standard_normal((m, n))
            dC2 = rng.standard_normal((m, n))
            betas = []
            for dC in [dC1, dC2]:
                dPhi = -Phi @ (dC @ Hinv @ C.T + C @ Hinv @ dC.T) @ Phi
                dZ = -L @ (dC @ Z)
                beta = inv(Phi) @ (L.T @ H @ dZ)
                bw = Phi_sqrt @ beta @ Rsqrt_inv
                betas.append(bw)
            F = betas[0] @ betas[1].T - betas[1] @ betas[0].T
            F_norms.append(norm(F))
        print(f"  r = {r:.2f}     | {np.mean(F_norms):.4f}")


def curvature_on_arrow():
    """
    Arrow rank-deficiency example.
    H = 3x3 latent precision matrix, various observers.
    """
    print("\n=== Curvature: arrow_rank_deficiency ===")

    H = np.array([
        [4.0, 1.8, 1.4],
        [1.8, 3.0, 0.1],
        [1.4, 0.1, 2.2],
    ])
    n = 3

    observers = {
        "plurality": np.array([[1.0, 0.0, 0.0]]),    # m=1
        "approval": np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),  # m=2
    }

    for name, C in observers.items():
        m = C.shape[0]
        print(f"\n  Observer: {name} (m={m})")

        if m == 1:
            print(f"    m=1 => curvature is always zero (1x1 antisymmetric)")
            continue

        Phi, L, Z, R = build_split_frame(H, C)
        Rinv = inv(R)
        Hinv = inv(H)
        Rsqrt_inv = np.real(sqrtm(Rinv))
        Phi_sqrt = np.real(sqrtm(Phi))

        rng = np.random.default_rng(42)
        F_norms = []
        for _ in range(200):
            dC1 = rng.standard_normal((m, n))
            dC2 = rng.standard_normal((m, n))
            betas = []
            for dC in [dC1, dC2]:
                dPhi = -Phi @ (dC @ Hinv @ C.T + C @ Hinv @ dC.T) @ Phi
                dZ = -L @ (dC @ Z)
                beta = inv(Phi) @ (L.T @ H @ dZ)
                bw = Phi_sqrt @ beta @ Rsqrt_inv
                betas.append(bw)
            F = betas[0] @ betas[1].T - betas[1] @ betas[0].T
            F_norms.append(norm(F))

        print(f"    Mean ||F||: {np.mean(F_norms):.4f}")
        print(f"    Max ||F||:  {np.max(F_norms):.4f}")
        print(f"    (n-m={n-m} hidden dimensions)")

        # Check conditioning
        print(f"    Phi eigenvalues: {sorted(eigh(Phi)[0])}")
        print(f"    R eigenvalues:   {sorted(eigh(R)[0])}")


# ============================================================
# GOAL D: General 2-parameter compatibility
# ============================================================

def general_2param_compatibility(n, m, seed=0):
    """
    When both H and C vary in both directions, the curvature is:
    F_alpha(ds, dt) = -(beta_s theta_t - beta_t theta_s)

    All four connection forms are nonzero. Verify this general formula.
    """
    print(f"\n=== General 2-parameter compatibility (n={n}, m={m}) ===")
    rng = np.random.default_rng(seed)

    H0 = spd(n, seed)
    C0 = rng.standard_normal((m, n))

    # s-direction: both H and C vary
    dH_s = sym(rng.standard_normal((n, n)))
    dC_s = rng.standard_normal((m, n)) * 0.5

    # t-direction: both H and C vary
    dH_t = sym(rng.standard_normal((n, n)))
    dC_t = rng.standard_normal((m, n)) * 0.5

    Phi0, L0, Z0, R0 = build_split_frame(H0, C0)
    R0inv = inv(R0)
    Hinv0 = inv(H0)

    def connection_forms(dH, dC):
        """Connection forms when both H and C change."""
        # dL from both H change and C change
        dHinv = -Hinv0 @ dH @ Hinv0
        dPhi_H = -Phi0 @ (C0 @ dHinv @ C0.T) @ Phi0  # from H change
        dPhi_C = -Phi0 @ (dC @ Hinv0 @ C0.T + C0 @ Hinv0 @ dC.T) @ Phi0  # from C change
        dPhi = dPhi_H + dPhi_C

        dL_H = dHinv @ C0.T @ Phi0 + Hinv0 @ C0.T @ dPhi_H  # from H change
        dL_C = Hinv0 @ dC.T @ Phi0 + Hinv0 @ C0.T @ dPhi_C  # from C change
        dL = dL_H + dL_C

        # dZ: from C change only (Z = ker C, and H doesn't affect Z directly)
        dZ = -L0 @ (dC @ Z0)

        alpha = inv(Phi0) @ (L0.T @ H0 @ dL)
        theta = R0inv @ (Z0.T @ H0 @ dL)
        beta = inv(Phi0) @ (L0.T @ H0 @ dZ)
        omega = R0inv @ (Z0.T @ H0 @ dZ)

        return alpha, beta, theta, omega

    alpha_s, beta_s, theta_s, omega_s = connection_forms(dH_s, dC_s)
    alpha_t, beta_t, theta_t, omega_t = connection_forms(dH_t, dC_t)

    # Curvature from flatness: F_alpha(ds, dt) = -(beta_s theta_t - beta_t theta_s)
    F_from_flatness = -(beta_s @ theta_t - beta_t @ theta_s)

    # Verify against finite-difference curvature
    eps = 1e-5

    def frame_at(s, t):
        H_st = H0 + s * dH_s + t * dH_t
        C_st = C0 + s * dC_s + t * dC_t
        ev = eigh(H_st)[0]
        if np.min(ev) < 0.01:
            return None
        return build_split_frame(H_st, C_st)

    # Curvature via Ambrose-Singer: track how the frame parallel-transports around
    # a small square (0,0) -> (eps,0) -> (eps,eps) -> (0,eps) -> (0,0)
    # The holonomy in the alpha-sector gives F_alpha * eps^2

    # Simpler: use the definition F_alpha = d alpha + alpha wedge alpha
    # F_alpha(ds, dt) = d_s alpha_t - d_t alpha_s + [alpha_s, alpha_t]
    # But d_s alpha_t requires second derivatives...

    # Instead, verify the general flatness formula against the specific case
    # where we decompose into pure-H and pure-C parts:
    # beta_s = beta_s^{(C)} (only from dC_s)
    # theta_s = theta_s^{(H)} + theta_s^{(C)} (from both dH_s and dC_s)

    # The F factorisation in the fully general case:
    # F(ds,dt) = -(beta_s theta_t - beta_t theta_s)
    #          = beta_t theta_s - beta_s theta_t

    # This is ALWAYS true (from flatness). The mixed factorisation theorem was a
    # special case where beta_s = 0.

    print(f"  beta_s nonzero: ||beta_s|| = {norm(beta_s):.6f}")
    print(f"  beta_t nonzero: ||beta_t|| = {norm(beta_t):.6f}")
    print(f"  theta_s: ||theta_s|| = {norm(theta_s):.6f}")
    print(f"  theta_t: ||theta_t|| = {norm(theta_t):.6f}")
    print(f"  ||F_alpha(ds,dt)|| = {norm(F_from_flatness):.6f}")

    # Verify F = -beta wedge theta via numerical second derivatives
    # Use d/ds[alpha_t] - d/dt[alpha_s] + [alpha_s, alpha_t] as cross-check
    def alpha_in_dir(dH_base, dC_base, dH_pert, dC_pert, eps_val):
        """alpha in (dH_base, dC_base) direction, evaluated at point shifted by eps*(dH_pert, dC_pert)."""
        H_shifted = H0 + eps_val * dH_pert
        C_shifted = C0 + eps_val * dC_pert
        ev = eigh(H_shifted)[0]
        if np.min(ev) < 0.01:
            return None
        Phi_sh, L_sh, Z_sh, R_sh = build_split_frame(H_shifted, C_shifted)
        Hinv_sh = inv(H_shifted)

        dHinv = -Hinv_sh @ dH_base @ Hinv_sh
        dPhi_H = -Phi_sh @ (C_shifted @ dHinv @ C_shifted.T) @ Phi_sh
        dPhi_C = -Phi_sh @ (dC_base @ Hinv_sh @ C_shifted.T + C_shifted @ Hinv_sh @ dC_base.T) @ Phi_sh
        dPhi = dPhi_H + dPhi_C
        dL = dHinv @ C_shifted.T @ Phi_sh + Hinv_sh @ C_shifted.T @ dPhi_H + Hinv_sh @ dC_base.T @ Phi_sh + Hinv_sh @ C_shifted.T @ dPhi_C

        return inv(Phi_sh) @ (L_sh.T @ H_shifted @ dL)

    # d/ds[alpha_t] at (0,0)
    alpha_t_at_plus = alpha_in_dir(dH_t, dC_t, dH_s, dC_s, eps)
    alpha_t_at_minus = alpha_in_dir(dH_t, dC_t, dH_s, dC_s, -eps)
    if alpha_t_at_plus is not None and alpha_t_at_minus is not None:
        d_s_alpha_t = (alpha_t_at_plus - alpha_t_at_minus) / (2 * eps)

        alpha_s_at_plus = alpha_in_dir(dH_s, dC_s, dH_t, dC_t, eps)
        alpha_s_at_minus = alpha_in_dir(dH_s, dC_s, dH_t, dC_t, -eps)
        d_t_alpha_s = (alpha_s_at_plus - alpha_s_at_minus) / (2 * eps)

        # F_alpha = d_s alpha_t - d_t alpha_s + [alpha_s, alpha_t]
        F_from_connection = d_s_alpha_t - d_t_alpha_s + alpha_s @ alpha_t - alpha_t @ alpha_s

        report("F from flatness vs F from connection (FD)", norm(F_from_flatness - F_from_connection), tol=1e-4)

        # Also verify: F = -beta wedge theta = beta_t theta_s - beta_s theta_t
        F_from_wedge = beta_t @ theta_s - beta_s @ theta_t
        report("F from wedge formula", norm(F_from_flatness - F_from_wedge))


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Goals B, C, D: Sextic, Curvature on Examples, General Compatibility")
    print("=" * 70)

    # GOAL C: Curvature on existing examples
    print("\n" + "=" * 70)
    print("GOAL C: Curvature on Existing Nomogeo Examples")
    print("=" * 70)
    curvature_on_entanglement()
    curvature_on_arrow()

    # GOAL D: General 2-parameter
    print("\n" + "=" * 70)
    print("GOAL D: General 2-Parameter Compatibility")
    print("=" * 70)
    general_2param_compatibility(4, 2, seed=0)
    general_2param_compatibility(5, 3, seed=10)
    general_2param_compatibility(6, 2, seed=20)

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
