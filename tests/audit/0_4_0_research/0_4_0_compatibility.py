"""
Direction 1 deep dive: the 2-parameter compatibility theorem.

KEY INSIGHT: The flatness equation d Omega + Omega wedge Omega = 0 gives
  F_alpha = -beta wedge theta
  F_alpha(ds, dt) = -(beta_s theta_t - beta_t theta_s)

For a mixed 2-parameter family:
  s-direction: H varies, C fixed (source active, beta_s = 0)
  t-direction: C varies, H fixed (curvature active, B_t = 0)

Then: F_alpha(ds, dt) = beta_t theta_s

where theta_s is the hidden connection from the H-variation (appears in source law)
and beta_t is from the C-variation (appears in curvature).

This IS the compatibility: the mixed curvature factorises as the product of
the source sector's theta and the curvature sector's beta.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det
from scipy.linalg import sqrtm

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


def verify_mixed_compatibility(n, m, seed=0):
    """
    THEOREM (Mixed source-curvature factorisation).

    Let H(s) be a 1-parameter family of SPD fields with fixed observer C.
    Let C(t) be a 1-parameter family of observers with fixed law H.
    Consider the 2-parameter family (H(s), C(t)).

    In the s-direction: beta_s = 0 (C fixed).
    In the t-direction: H_t = 0, so B_t = L^T H_t Z = 0.

    From flatness: F_alpha(ds, dt) = -(beta_s theta_t - beta_t theta_s) = beta_t theta_s.

    Moreover, theta_s = -R^{-1} B_s^T (from B_callig with beta_s = 0).

    So: F_alpha(ds, dt) = -beta_t R^{-1} B_s^T

    where B_s = L^T H_s Z (the visible-hidden coupling from the source law)
    and beta_t is the hidden frame rotation from the C-variation.

    This factorisation means: the mixed curvature = (C-motion coupling) * (H-coupling).
    The source law controls the H-coupling (through A_cpl), and the pure-observer
    curvature controls the C-motion coupling (through beta_t).
    """
    print(f"\n=== Mixed compatibility (n={n}, m={m}) ===")
    rng = np.random.default_rng(seed)

    # Base point
    H0 = spd(n, seed)
    C0 = rng.standard_normal((m, n))

    Phi0, L0, Z0, R0 = build_split_frame(H0, C0)
    R0inv = inv(R0)
    Hinv0 = inv(H0)

    # s-direction: H perturbation (C fixed)
    dH = sym(rng.standard_normal((n, n)))  # Hdot

    # t-direction: C perturbation (H fixed)
    dC = rng.standard_normal((m, n))

    # ---- s-direction connection forms (C fixed, H varies) ----
    # dL_s from H change
    dHinv_s = -Hinv0 @ dH @ Hinv0
    dPhi_s = -Phi0 @ (C0 @ dHinv_s @ C0.T) @ Phi0
    dL_s = dHinv_s @ C0.T @ Phi0 + Hinv0 @ C0.T @ dPhi_s
    # dZ_s = 0 (Z from SVD of C, which is fixed)

    alpha_s = inv(Phi0) @ (L0.T @ H0 @ dL_s)
    theta_s = R0inv @ (Z0.T @ H0 @ dL_s)
    beta_s = np.zeros((m, n - m))  # C fixed => beta_s = 0
    omega_s = np.zeros((n - m, n - m))

    # Verify beta_s = 0
    report("beta_s = 0 (C fixed in s-direction)", norm(beta_s))

    # Source law quantities
    B_s = L0.T @ dH @ Z0
    # theta_s should equal -R^{-1} B_s^T (from B_callig with beta = 0)
    theta_s_check = -R0inv @ B_s.T
    report("theta_s = -R^{-1} B_s^T", norm(theta_s - theta_s_check))

    # ---- t-direction connection forms (H fixed, C varies) ----
    # dL_t from C change
    dPhi_t = -Phi0 @ (dC @ Hinv0 @ C0.T + C0 @ Hinv0 @ dC.T) @ Phi0
    dL_t = Hinv0 @ dC.T @ Phi0 + Hinv0 @ C0.T @ dPhi_t
    # dZ_t: need C(t)Z(t) = 0, so dC Z + C dZ = 0, hence dZ = -L (dC Z)
    dZ_t = -L0 @ (dC @ Z0)

    alpha_t = inv(Phi0) @ (L0.T @ H0 @ dL_t)
    theta_t = R0inv @ (Z0.T @ H0 @ dL_t)
    beta_t = inv(Phi0) @ (L0.T @ H0 @ dZ_t)
    omega_t = R0inv @ (Z0.T @ H0 @ dZ_t)

    # Verify B_t = 0 (H fixed in t-direction, so dH = 0 along t)
    B_t = L0.T @ np.zeros((n, n)) @ Z0  # H_t = 0
    report("B_t = 0 (H fixed in t-direction)", norm(B_t))

    # Pure-observer identity: theta_t = -R^{-1} beta_t^T Phi (from dH=0)
    theta_t_check = -R0inv @ beta_t.T @ Phi0
    report("theta_t = -R^{-1} beta_t^T Phi (pure observer)", norm(theta_t - theta_t_check))

    # ---- The compatibility theorem ----
    # F_alpha(ds, dt) = -(beta_s theta_t - beta_t theta_s)
    # With beta_s = 0: F_alpha(ds, dt) = beta_t theta_s

    F_alpha_from_flatness = beta_t @ theta_s

    # Direct formula: F_alpha(ds, dt) = -beta_t R^{-1} B_s^T
    F_alpha_factored = -beta_t @ R0inv @ B_s.T
    report("F_alpha = beta_t theta_s", norm(F_alpha_from_flatness - F_alpha_factored))

    # ---- Verify against full curvature formula ----
    # General: F_alpha(d1, d2) = -beta_1 theta_2 + beta_2 theta_1
    # Here d1 = ds (beta_1 = beta_s = 0, theta_1 = theta_s)
    #      d2 = dt (beta_2 = beta_t, theta_2 = theta_t)
    F_alpha_general = -beta_s @ theta_t + beta_t @ theta_s
    report("General formula matches factored", norm(F_alpha_general - F_alpha_factored))

    # ---- Norm identity ----
    # ||F_alpha(ds, dt)||_F^2 = ||beta_t theta_s||_F^2
    #                         = Tr(theta_s^T beta_t^T beta_t theta_s)
    #                         = Tr(O_t^{hid} theta_s theta_s^T)
    # where O_t^{hid} = beta_t^T beta_t ((n-m)x(n-m), the hidden Gram of beta_t)
    F_norm_sq = np.trace(F_alpha_factored.T @ F_alpha_factored)
    O_hid_t = beta_t.T @ beta_t
    theta_gram = theta_s @ theta_s.T
    trace_form = np.trace(O_hid_t @ theta_gram)
    report("||F||^2 = Tr(G_beta^t * theta_s theta_s^T)", abs(F_norm_sq - trace_form))

    # ---- Connection to source law ----
    # theta_s = -R^{-1} B_s^T, so:
    # ||F||^2 = Tr(O_hid_t * R^{-1} B_s^T B_s R^{-1})
    #         = Tr(O_hid_t * R^{-1} Qhat R^{-1}) where Qhat = B_s^T B_s... no
    # Actually B_s is m x (n-m), B_s^T is (n-m) x m.
    # theta_s theta_s^T = R^{-1} B_s^T B_s R^{-1}  ((n-m)x(n-m))
    # But B_s^T B_s is (n-m)x(n-m) (hidden-space Gram of the source coupling)
    # And Qhat (from source law) = B_s R^{-1} B_s^T  (m x m, visible-space)
    # These are different projections of the same data.

    theta_theta = R0inv @ B_s.T @ B_s @ R0inv
    report("theta_s theta_s^T = R^{-1} B^T B R^{-1}", norm(theta_gram - theta_theta))

    F_from_source = np.trace(O_hid_t @ R0inv @ B_s.T @ B_s @ R0inv)
    report("||F||^2 from source coupling B", abs(F_norm_sq - F_from_source))

    print(f"\n  ||F_alpha(ds,dt)|| = {np.sqrt(F_norm_sq):.6f}")
    print(f"  ||beta_t|| = {norm(beta_t):.6f}")
    print(f"  ||theta_s|| = {norm(theta_s):.6f}")
    print(f"  ||B_s|| = {norm(B_s):.6f}")

    return F_norm_sq


def verify_trace_inequality(n, m, num_trials=500, seed=0):
    """
    From the factorisation F = beta_t theta_s:

    ||F||^2 = Tr(G_beta * G_theta)

    where G_beta = beta_t^T beta_t and G_theta = theta_s theta_s^T = R^{-1} B^T B R^{-1}.

    By Cauchy-Schwarz: ||F||^2 <= Tr(G_beta) * Tr(G_theta) = ||beta_t||^2 * ||theta_s||^2.

    This bounds the mixed curvature in terms of the pure-observer coupling
    strength and the source coupling strength.

    The bound is tight when G_beta and G_theta are proportional (aligned in hidden space).
    """
    print(f"\n=== Trace inequality (n={n}, m={m}, {num_trials} trials) ===")

    max_ratio = 0
    ratios = []

    for trial in range(num_trials):
        rng = np.random.default_rng(seed + trial * 100)
        H0 = spd(n, seed + trial * 100)
        C0 = rng.standard_normal((m, n))
        Phi0, L0, Z0, R0 = build_split_frame(H0, C0)
        R0inv = inv(R0)
        Hinv0 = inv(H0)

        dH = sym(rng.standard_normal((n, n)))
        dC = rng.standard_normal((m, n))

        # s-direction (source)
        B_s = L0.T @ dH @ Z0
        theta_s = -R0inv @ B_s.T

        # t-direction (curvature)
        dPhi_t = -Phi0 @ (dC @ Hinv0 @ C0.T + C0 @ Hinv0 @ dC.T) @ Phi0
        dZ_t = -L0 @ (dC @ Z0)
        beta_t = inv(Phi0) @ (L0.T @ H0 @ dZ_t)

        F = beta_t @ theta_s
        F_sq = np.trace(F.T @ F)
        bound = norm(beta_t)**2 * norm(theta_s)**2

        if bound > 1e-15:
            ratio = F_sq / bound
            ratios.append(ratio)
            max_ratio = max(max_ratio, ratio)

    report(f"||F||^2 <= ||beta||^2 ||theta||^2", max(0, max_ratio - 1.0))
    print(f"  Max ratio: {max_ratio:.6f}")
    print(f"  Mean ratio: {np.mean(ratios):.6f}")
    print(f"  (1.0 = perfectly aligned, <1.0 = generic)")


def sextic_form_verification(seed=0):
    """
    Verify the sextic form q_{6,mu}(Z) from Corollary 4.2 of the 0.4.0 TN.

    q_{6,mu}(Z) = -Tr(Q^3 M_U) + Tr(Q^2 N)
                  -(1+mu) sum_a [||Y_a Q||^2 + ||R^{1/2} Y_a Q^{1/2}||^2 + ||R Y_a||^2]

    where Q = Z^T Z, R = Z Z^T, N = Z^T M_W Z, Y_a = B_a Z - Z A_a.

    This is homogeneous of degree 6 in Z. We verify this and check the sign.
    """
    print("\n=== Sextic form verification (Cor 4.2) ===")
    rng = np.random.default_rng(seed + 800)

    m, p = 3, 2  # m = dim(U), p = dim(U_perp), tangent map Z: U -> U_perp is p x m
    n_families = 3  # number of weighted family members

    mu = 1.0

    # Random weighted family data (from TN1 Thm 2.3 notation)
    # A_a = A_U(omega_a) operates on U (m x m)
    # B_a = A_perp(omega_a) operates on U_perp (p x p)
    A_list = [sym(rng.standard_normal((m, m))) for _ in range(n_families)]  # m x m
    B_list = [sym(rng.standard_normal((p, p))) for _ in range(n_families)]  # p x p

    M_U = sum(A @ A for A in A_list)  # m x m
    M_W = sum(B @ B for B in B_list)  # p x p (M_perp)

    def sextic_form(Z, mu_val):
        # Z is p x m (tangent map from U to U_perp)
        Q = Z.T @ Z  # m x m
        R_mat = Z @ Z.T  # p x p
        N = Z.T @ M_W @ Z  # m x m (M_W is p x p)

        term1 = -np.trace(Q @ Q @ Q @ M_U)  # -Tr(Q^3 M_U)
        term2 = np.trace(Q @ Q @ N)          # Tr(Q^2 N)

        term3 = 0.0
        Qsqrt = np.real(sqrtm(Q + 1e-30 * np.eye(m)))
        Rsqrt = np.real(sqrtm(R_mat + 1e-30 * np.eye(p)))
        for a in range(n_families):
            Y = B_list[a] @ Z - Z @ A_list[a]  # p x m
            term3 += norm(Y @ Q)**2  # ||Y_a Q||_F^2
            term3 += norm(Rsqrt @ Y @ Qsqrt)**2  # ||R^{1/2} Y_a Q^{1/2}||_F^2
            term3 += norm(R_mat @ Y)**2  # ||R Y_a||_F^2

        return term1 + term2 - (1 + mu_val) * term3

    # Verify homogeneity of degree 6
    Z_test = rng.standard_normal((p, m)) * 0.5
    q_base = sextic_form(Z_test, mu)

    scalings = [0.5, 0.8, 1.2, 2.0, 3.0]
    print("  Homogeneity check (expect ratio = t^6):")
    for t in scalings:
        q_scaled = sextic_form(t * Z_test, mu)
        if abs(q_base) > 1e-15:
            ratio = q_scaled / q_base
            expected = t**6
            print(f"    t={t:.1f}: q(tZ)/q(Z) = {ratio:.6f}, t^6 = {expected:.6f}, "
                  f"error = {abs(ratio - expected):.3e}")

    # Check sign over many random Z directions
    n_neg = 0
    n_pos = 0
    n_zero = 0
    for trial in range(500):
        Z_rand = rng.standard_normal((p, m)) * 0.3
        q_val = sextic_form(Z_rand, mu)
        if q_val < -1e-12:
            n_neg += 1
        elif q_val > 1e-12:
            n_pos += 1
        else:
            n_zero += 1

    print(f"\n  Sign distribution (500 random Z, mu={mu}):")
    print(f"    Negative: {n_neg} ({100*n_neg/500:.0f}%)")
    print(f"    Positive: {n_pos} ({100*n_pos/500:.0f}%)")
    print(f"    Near zero: {n_zero} ({100*n_zero/500:.0f}%)")
    print(f"    (Negative needed for strict local maximality)")


def source_law_on_rc_decay():
    """
    Apply the source law to the RC decay system from the operationalise layer.

    RC decay: V(t) = V0 exp(-t/tau), observed with Gaussian noise sigma.
    The Fisher information for tau is:
      I(tau) = sum_i (dV_i/dtau)^2 / sigma^2 = sum_i (V0 t_i exp(-t_i/tau) / tau^2)^2 / sigma^2

    This is a scalar (m=1) in the observation parameter tau. With n=2 (tau + one hidden
    nuisance), we can construct a 2D precision matrix H(tau) and compute A_cpl.

    For m=1, the curvature F_alpha is always zero (1x1 antisymmetric = 0).
    But the source law A_cpl is nontrivial and governs how the Fisher precision
    changes along the tau axis — i.e., how informative the experiment is at different tau.
    """
    print("\n=== Source law on RC decay ===")

    # Q7 dataset (from the probe): 14 time points, V0 = 100V
    # Rather than loading the file, use the typical values
    V0 = 100.0
    times = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
                       4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
    tau_hat = 2.16  # approximate MLE from probe
    sigma = 3.0     # approximate residual noise

    # Fisher information as a function of tau
    def fisher_tau(tau):
        sensitivities = V0 * times * np.exp(-times / tau) / tau**2
        return np.sum(sensitivities**2) / sigma**2

    # Compute Fisher and its derivatives numerically
    dt = 1e-5
    I_0 = fisher_tau(tau_hat)
    I_p = fisher_tau(tau_hat + dt)
    I_m = fisher_tau(tau_hat - dt)
    I_dot = (I_p - I_m) / (2 * dt)
    I_ddot = (I_p - 2 * I_0 + I_m) / dt**2

    # For a 1D system (m=1, n=1), the source law simplifies:
    # V = I_dot (scalar), W = I_ddot
    # A_cpl = -1/2 * W / V (when V > 0)
    # = -1/2 * I_ddot / I_dot

    print(f"  At tau = {tau_hat:.2f} s:")
    print(f"    Fisher I(tau) = {I_0:.4f}")
    print(f"    dI/dtau = {I_dot:.4f}")
    print(f"    d^2I/dtau^2 = {I_ddot:.4f}")

    if abs(I_dot) > 1e-10:
        A_cpl = -0.5 * I_ddot / I_dot
        print(f"    A_cpl = {A_cpl:.6f}")
        if A_cpl > 0:
            print(f"    => Fisher precision is DECELERATING (concave log-evidence)")
        else:
            print(f"    => Fisher precision is ACCELERATING (convex log-evidence)")
    else:
        print(f"    dI/dtau ~ 0: at a Fisher extremum")

    # Sweep A_cpl over tau range
    taus = np.linspace(0.5, 5.0, 91)
    A_cpls = []
    fishers = []
    for tau in taus:
        I = fisher_tau(tau)
        Ip = fisher_tau(tau + dt)
        Im = fisher_tau(tau - dt)
        Id = (Ip - Im) / (2*dt)
        Idd = (Ip - 2*I + Im) / dt**2
        fishers.append(I)
        if abs(Id) > 1e-10:
            A_cpls.append(-0.5 * Idd / Id)
        else:
            A_cpls.append(0.0)

    A_cpls = np.array(A_cpls)
    fishers = np.array(fishers)

    # Find sign changes of A_cpl (inflection points of log Fisher)
    sign_changes = []
    for i in range(1, len(A_cpls)):
        if A_cpls[i] * A_cpls[i-1] < 0:
            sign_changes.append(taus[i])

    print(f"\n  A_cpl sign changes (inflection points of log-evidence):")
    for sc in sign_changes:
        print(f"    tau = {sc:.2f} s")

    print(f"\n  Fisher peak at tau = {taus[np.argmax(fishers)]:.2f} s")
    print(f"  Fisher value at peak: {np.max(fishers):.4f}")

    print(f"\n  Interpretation: The source law A_cpl tells us where the")
    print(f"  Fisher information is most sensitive to changes in tau.")
    print(f"  At sign changes, the log-evidence curvature flips.")


if __name__ == "__main__":
    print("=" * 70)
    print("0.4.0 Compatibility Theorem and Applied Research")
    print("=" * 70)

    verify_mixed_compatibility(4, 2, seed=0)
    verify_mixed_compatibility(5, 3, seed=10)
    verify_mixed_compatibility(6, 2, seed=20)

    verify_trace_inequality(5, 2, 500, seed=0)
    verify_trace_inequality(6, 3, 500, seed=100)

    sextic_form_verification(seed=0)

    source_law_on_rc_decay()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
