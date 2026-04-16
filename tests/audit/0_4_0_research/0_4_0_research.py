"""
Research code for 0.4.0 follow-up directions.

Direction 1: Compatibility identity — exact relation between curvature and
             hidden-space Gram geometry.
Direction 2: Finite-epsilon corrections to the fast-hidden lift.
Direction 3: Connection to the TN4 typed evidence stack.
Direction 4: Curvature as observer quality metric.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det
from scipy.linalg import sqrtm, expm

TOL = 1e-10
passes = 0
fails = 0

def spd(n, seed=None):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T + np.eye(n)

def sym(M):
    return 0.5 * (M + M.T)

def report(name, err, tol=TOL):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok

def build_split_frame(H, C):
    n = H.shape[0]; m = C.shape[0]
    Hinv = inv(H); Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    U, s, Vt = svd(C); Z = Vt[m:].T
    R = Z.T @ H @ Z
    return Phi, L, Z, R

def compute_connection_exact(H, C, Hdot, Z):
    Hinv = inv(H); Phi = inv(C @ Hinv @ C.T); L = Hinv @ C.T @ Phi; R = Z.T @ H @ Z
    dHinv = -Hinv @ Hdot @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi
    alpha = inv(Phi) @ (L.T @ H @ dL)
    theta = inv(R) @ (Z.T @ H @ dL)
    return alpha, theta, dL, dPhi


# ============================================================
# DIRECTION 1: Source-Curvature Compatibility
# ============================================================

def verify_compatibility_identity(num_trials=1000):
    """
    THEOREM (Curvature-Gram Identity, whitened pure-observer branch):

    Let beta_1, beta_2 in R^{m x (n-m)} be tangent directions at a reference
    observer in Gr(m,n). In whitened gauge (Phi=I, R=I):

      F_alpha = beta_1 beta_2^T - beta_2 beta_1^T   (m x m, antisymmetric)

    Define hidden-space Gram matrices:
      G_i = beta_i^T beta_i   ((n-m) x (n-m), PSD)
      C   = beta_1^T beta_2   ((n-m) x (n-m))

    Then:
      ||F_alpha||_F^2 = 2(Tr(G_1 G_2) - Tr(C^2))

    PROOF:
      ||F||^2 = Tr(F^T F) = Tr((B-A)(A-B)^T) with A = beta_1 beta_2^T, B = beta_2 beta_1^T.
      Expanding and using cyclic trace:
        = 2 Tr(beta_2 beta_1^T beta_1 beta_2^T) - 2 Tr((beta_1^T beta_2)^2)
        = 2 Tr(G_1 G_2) - 2 Tr(C^2).                                       QED.

    COROLLARY (Cauchy-Schwarz bound):
      ||F_alpha||_F^2 <= 2 Tr(G_1) Tr(G_2) = 2 ||beta_1||_F^2 ||beta_2||_F^2 = 2 Tr(O_1) Tr(O_2)

    where O_i = beta_i beta_i^T (m x m) is the observer tensor.
    """
    print("--- Compatibility identity: ||F||^2 = 2(Tr(G1 G2) - Tr(C^2)) ---")

    max_err = 0
    max_ineq_violation = 0
    for trial in range(num_trials):
        rng = np.random.default_rng(trial)
        n = rng.integers(3, 10)
        m = rng.integers(1, n)

        beta1 = rng.standard_normal((m, n - m))
        beta2 = rng.standard_normal((m, n - m))

        F = beta1 @ beta2.T - beta2 @ beta1.T
        F_sq = np.trace(F.T @ F)

        G1 = beta1.T @ beta1
        G2 = beta2.T @ beta2
        C = beta1.T @ beta2

        identity_rhs = 2.0 * (np.trace(G1 @ G2) - np.trace(C @ C))
        max_err = max(max_err, abs(F_sq - identity_rhs))

        # Cauchy-Schwarz bound
        bound = 2.0 * np.trace(G1) * np.trace(G2)
        max_ineq_violation = max(max_ineq_violation, F_sq - bound)

    report(f"Identity over {num_trials} trials", max_err)
    report(f"Cauchy-Schwarz bound (violation)", max(0, max_ineq_violation))
    return max_err < 1e-10


def verify_curvature_vanishing():
    """
    F_alpha = 0 iff Tr(G_1 G_2) = Tr(C^2).

    Special case: when m = 1, F is always zero (1x1 antisymmetric = 0).
    For m >= 2, F = 0 requires a specific alignment between hidden Gram matrices.
    """
    print("\n--- Curvature vanishing analysis ---")

    # m=1: always flat
    rng = np.random.default_rng(42)
    n, m = 5, 1
    max_F = 0
    for _ in range(100):
        b1 = rng.standard_normal((m, n-m))
        b2 = rng.standard_normal((m, n-m))
        F = b1 @ b2.T - b2 @ b1.T
        max_F = max(max_F, norm(F))
    report(f"m=1 always flat (max ||F|| = {max_F:.3e})", max_F)

    # m=2, n=4: construct commuting betas (same right singular vectors)
    n, m = 4, 2
    V_base = np.linalg.qr(rng.standard_normal((n-m, n-m)))[0]
    s1 = np.array([2.0, 1.0])
    s2 = np.array([3.0, 0.5])
    U1 = np.linalg.qr(rng.standard_normal((m, m)))[0]
    U2 = U1  # same left singular vectors
    beta1 = U1 @ np.diag(s1) @ V_base.T
    beta2 = U2 @ np.diag(s2) @ V_base.T

    G1 = beta1.T @ beta1; G2 = beta2.T @ beta2; C = beta1.T @ beta2
    F = beta1 @ beta2.T - beta2 @ beta1.T
    deficit = np.trace(G1 @ G2) - np.trace(C @ C)
    print(f"  Same left-SVD: ||F|| = {norm(F):.6f}, deficit = {deficit:.6f}")

    # Commuting G1, G2 with same eigenbasis: C = G1^{1/2} Q G2^{1/2} for some Q
    # If beta_i share the SAME right singular vectors, then C = V diag(s1*s2) V^T which is symmetric.
    # Then C^2 has eigenvalues (s1_i * s2_i)^2 and G1 G2 has eigenvalues s1_i^2 * s2_i^2 (if same eigenbasis).
    # So Tr(C^2) = sum (s1_i s2_i)^2 = Tr(G1 G2). F = 0!

    # But to get F = 0 for m >= 2, we also need beta1 beta2^T = beta2 beta1^T.
    # beta1 beta2^T = U1 diag(s1) V^T V diag(s2) U2^T = U1 diag(s1*s2) U1^T (since U2=U1)
    # beta2 beta1^T = U1 diag(s2*s1) U1^T. Same! So F = 0.

    report("Commuting betas give F=0", norm(F))

    # Now perturb to break commutativity
    U2_pert = U2 @ expm(0.1 * np.array([[0, 1], [-1, 0]]))
    beta2_pert = U2_pert @ np.diag(s2) @ V_base.T
    F_pert = beta1 @ beta2_pert.T - beta2_pert @ beta1.T
    print(f"  Perturbed (noncommuting): ||F|| = {norm(F_pert):.6f}")
    report("Noncommuting betas give F != 0", 0.0 if norm(F_pert) > 0.01 else 1.0)


def verify_general_gauge(num_trials=500):
    """
    In the general gauge, define whitened betas:
      beta_i^{(w)} = Phi^{1/2} beta_i R^{-1/2}

    Then F_alpha = Phi^{-1/2} F^{(w)} Phi^{1/2} where
    F^{(w)} = beta_1^{(w)} beta_2^{(w)T} - beta_2^{(w)} beta_1^{(w)T}

    and the identity applies in whitened coordinates:
      ||F^{(w)}||_F^2 = 2(Tr(G_1^{(w)} G_2^{(w)}) - Tr(C_w^2))

    where G_i^{(w)} = beta_i^{(w)T} beta_i^{(w)} and C_w = beta_1^{(w)T} beta_2^{(w)}.
    """
    print(f"\n--- General gauge identity ({num_trials} trials) ---")

    max_err = 0
    for trial in range(num_trials):
        rng = np.random.default_rng(trial + 50000)
        n = rng.integers(4, 8); m = rng.integers(2, n-1)
        H = spd(n, trial + 50000)
        C_obs = rng.standard_normal((m, n))
        Phi, L, Z, R = build_split_frame(H, C_obs)
        Rsqrt_inv = np.real(sqrtm(inv(R)))
        Phi_sqrt = np.real(sqrtm(Phi))

        beta1 = rng.standard_normal((m, n - m))
        beta2 = rng.standard_normal((m, n - m))

        # Whitened betas
        bw1 = Phi_sqrt @ beta1 @ Rsqrt_inv
        bw2 = Phi_sqrt @ beta2 @ Rsqrt_inv

        # Whitened curvature
        Fw = bw1 @ bw2.T - bw2 @ bw1.T
        Fw_sq = np.trace(Fw.T @ Fw)

        Gw1 = bw1.T @ bw1; Gw2 = bw2.T @ bw2; Cw = bw1.T @ bw2
        rhs = 2.0 * (np.trace(Gw1 @ Gw2) - np.trace(Cw @ Cw))
        max_err = max(max_err, abs(Fw_sq - rhs))

    report(f"Whitened identity in general gauge", max_err)


# ============================================================
# DIRECTION 2: Finite-Epsilon Fast-Hidden Lift
# ============================================================

def finite_epsilon_analysis(n, m, seed=0):
    """
    For constant-coefficient system (A, R, Theta fixed):
      dX/dt = 1/2 A X - 1/2 R Y
      eps dY/dt = -Y + Theta X

    Slow manifold: Y = Theta X + eps Theta G X + O(eps^2) where G = 1/2(A - R Theta).
    Corrected generator: G_eff = G - eps/2 R Theta G.
    Corrected Gram: Pdot = G_eff + G_eff^T = A_cpl + eps * Delta_1 + O(eps^2).
    Delta_1 = -(R Theta G + G^T Theta^T R^T) / 2.
    """
    print(f"\n--- Finite-eps (n={n}, m={m}) ---")
    rng = np.random.default_rng(seed + 300)

    for attempt in range(100):
        H0 = spd(n, seed + attempt * 500)
        H1 = sym(rng.standard_normal((n, n)) * 5.0)
        H2 = sym(rng.standard_normal((n, n)))
        C = rng.standard_normal((m, n))
        Phi, L, Z, R_mat = build_split_frame(H0, C)
        V = L.T @ H1 @ L
        if np.min(eigh(V)[0]) > 0.5:
            break
    else:
        print("  Skipped: V > 0 not found"); return

    B = L.T @ H1 @ Z
    alpha_t, theta_t, dL, dPhi = compute_connection_exact(H0, C, H1, Z)
    Vsqrt = np.real(sqrtm(V)); Vsqrt_inv = inv(Vsqrt)

    A_vis = -0.5 * Vsqrt_inv @ (L.T @ (2*H2) @ L) @ Vsqrt_inv
    R_term = Vsqrt_inv @ B
    Theta_term = theta_t @ Vsqrt_inv

    G = 0.5 * (A_vis - R_term @ Theta_term)
    A_cpl = G + G.T

    # First-order correction
    Delta_1 = -0.5 * (R_term @ Theta_term @ G + G.T @ Theta_term.T @ R_term.T)

    epsilons = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
    err_0 = []
    err_1 = []

    for eps_val in epsilons:
        X = np.eye(m); Y = Theta_term.copy()
        dt_sim = min(eps_val * 0.005, 5e-6)
        T_final = 5e-5
        n_steps = max(int(T_final / dt_sim), 20)
        dt_sim = T_final / n_steps

        for _ in range(n_steps):
            def f(XX, YY):
                return (0.5 * A_vis @ XX - 0.5 * R_term @ YY,
                        (-YY + Theta_term @ XX) / eps_val)
            k1x, k1y = f(X, Y)
            k2x, k2y = f(X + .5*dt_sim*k1x, Y + .5*dt_sim*k1y)
            k3x, k3y = f(X + .5*dt_sim*k2x, Y + .5*dt_sim*k2y)
            k4x, k4y = f(X + dt_sim*k3x, Y + dt_sim*k3y)
            X += dt_sim/6*(k1x + 2*k2x + 2*k3x + k4x)
            Y += dt_sim/6*(k1y + 2*k2y + 2*k3y + k4y)

        Pdot = sym((X @ X.T - np.eye(m)) / T_final)
        err_0.append(norm(Pdot - A_cpl))
        err_1.append(norm(Pdot - A_cpl - eps_val * Delta_1))

    print("  eps     ||Pdot-A_cpl||  ||Pdot-A_cpl-eps*D1||")
    for i, eps in enumerate(epsilons):
        print(f"  {eps:.3f}  {err_0[i]:.4e}       {err_1[i]:.4e}")

    # Convergence rates
    print("  Rates:")
    for i in range(1, len(epsilons)):
        if err_0[i] > 1e-14 and err_0[i-1] > 1e-14:
            r0 = np.log(err_0[i]/err_0[i-1]) / np.log(epsilons[i]/epsilons[i-1])
            r1 = np.log(max(err_1[i],1e-16)/max(err_1[i-1],1e-16)) / np.log(epsilons[i]/epsilons[i-1])
            print(f"  eps={eps:.3f}: uncorrected ~eps^{r0:.1f}, corrected ~eps^{r1:.1f}")

    return err_0, err_1


# ============================================================
# DIRECTION 3: Evidence curvature
# ============================================================

def evidence_curvature(n, m, seed=0):
    """
    f(t) = log det Phi(t) is the log-evidence-volume along a path H(t).
    f'(0) = Tr(Phi^{-1} dPhi) = Tr(Phi^{-1} V) + 2 Tr(alpha).
    f''(0) involves W and hence A_cpl.
    """
    print(f"\n--- Evidence curvature (n={n}, m={m}) ---")
    rng = np.random.default_rng(seed + 500)

    H0 = spd(n, seed + 500)
    H1 = sym(rng.standard_normal((n, n)) * 2.0)
    H2 = sym(rng.standard_normal((n, n)))
    C = rng.standard_normal((m, n))

    ts = np.linspace(-0.15, 0.15, 301)
    f_vals = []
    for t in ts:
        H_t = H0 + t*H1 + t**2*H2
        if np.min(eigh(H_t)[0]) < 0.01:
            f_vals.append(np.nan); continue
        f_vals.append(np.log(det(inv(C @ inv(H_t) @ C.T))))
    f_vals = np.array(f_vals)

    dt = ts[1]-ts[0]; mid = len(ts)//2
    fp_num = (f_vals[mid+1]-f_vals[mid-1])/(2*dt)
    fpp_num = (f_vals[mid+1]-2*f_vals[mid]+f_vals[mid-1])/dt**2

    Phi, L, Z, R = build_split_frame(H0, C)
    dHinv = -inv(H0) @ H1 @ inv(H0)
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    fp_exact = np.trace(inv(Phi) @ dPhi)

    report(f"f'(0) numerical vs exact", abs(fp_num - fp_exact), tol=1e-4)
    print(f"  f'(0) = {fp_exact:.6f}")
    print(f"  f''(0) = {fpp_num:.6f}")

    V = L.T @ H1 @ L
    alpha_t, theta_t, dL, _ = compute_connection_exact(H0, C, H1, Z)
    fp_split = np.trace(inv(Phi) @ V) + 2*np.trace(alpha_t)
    report(f"f'(0) = Tr(Phi^-1 V) + 2Tr(alpha)", abs(fp_exact - fp_split))

    return fp_exact, fpp_num


# ============================================================
# DIRECTION 4: Curvature spectrum
# ============================================================

def curvature_spectrum(n, m, n_dirs=40, seed=0):
    """
    Compute the distribution of curvature norms across random observer
    perturbation pairs. This characterises the geometry of the observer space.
    """
    print(f"\n--- Curvature spectrum (n={n}, m={m}, {n_dirs} directions) ---")
    rng = np.random.default_rng(seed + 700)

    H = spd(n, seed + 700)
    C0 = rng.standard_normal((m, n))
    Phi0, L0, Z0, R0 = build_split_frame(H, C0)
    Rsqrt_inv = np.real(sqrtm(inv(R0)))
    Phi_sqrt = np.real(sqrtm(Phi0))
    Hinv = inv(H)

    # Extract whitened betas for random perturbation directions
    betas_w = []
    for _ in range(n_dirs):
        dC = rng.standard_normal((m, n))
        dPhi = -Phi0 @ (dC @ Hinv @ C0.T + C0 @ Hinv @ dC.T) @ Phi0
        dZ = -L0 @ (dC @ Z0)
        beta = inv(Phi0) @ (L0.T @ H @ dZ)
        bw = Phi_sqrt @ beta @ Rsqrt_inv
        betas_w.append(bw)

    # Curvature for all pairs
    from itertools import combinations
    F_norms = []
    for i, j in combinations(range(n_dirs), 2):
        F = betas_w[i] @ betas_w[j].T - betas_w[j] @ betas_w[i].T
        F_norms.append(norm(F))

    F_norms = np.array(F_norms)
    n_pairs = len(F_norms)
    n_flat = np.sum(F_norms < 1e-8)

    print(f"  {n_pairs} pairs: mean ||F|| = {np.mean(F_norms):.4f}, "
          f"max = {np.max(F_norms):.4f}, min = {np.min(F_norms):.4f}")
    print(f"  Flat pairs: {n_flat}/{n_pairs} ({100*n_flat/n_pairs:.0f}%)")
    print(f"  Percentiles: p25={np.percentile(F_norms,25):.4f}, "
          f"p50={np.percentile(F_norms,50):.4f}, p75={np.percentile(F_norms,75):.4f}")

    return F_norms


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("0.4.0 Research — Numerical Results")
    print("=" * 70)

    print("\n" + "=" * 70)
    print("DIRECTION 1: Curvature-Gram Compatibility Identity")
    print("=" * 70)
    verify_compatibility_identity(2000)
    verify_curvature_vanishing()
    verify_general_gauge(500)

    print("\n" + "=" * 70)
    print("DIRECTION 2: Finite-Epsilon Correction to Fast-Hidden Lift")
    print("=" * 70)
    finite_epsilon_analysis(4, 2, seed=0)
    finite_epsilon_analysis(5, 3, seed=10)

    print("\n" + "=" * 70)
    print("DIRECTION 3: Source Law and Evidence Geometry")
    print("=" * 70)
    evidence_curvature(4, 2, seed=0)
    evidence_curvature(5, 3, seed=10)

    print("\n" + "=" * 70)
    print("DIRECTION 4: Curvature Spectrum")
    print("=" * 70)
    curvature_spectrum(5, 2, 40, seed=0)
    curvature_spectrum(6, 3, 40, seed=10)

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
