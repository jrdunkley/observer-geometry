"""
Goal A: Evidence second derivative from A_cpl.

The log-determinant of visible precision along a path H(t) with fixed C:
  f(t) = log det Phi(t)

First derivative (proved in research report):
  f'(t) = Tr(Phi^{-1} dPhi/dt) = Tr(Phi^{-1} V) + 2 Tr(alpha)

Second derivative (to be derived and verified):
  f''(t) = d/dt [Tr(Phi^{-1} dPhi/dt)]

This should be expressible in terms of A_cpl, V, W, alpha.

DERIVATION:
  f'(t) = Tr(Phi^{-1} dPhi)

  f''(t) = d/dt Tr(Phi^{-1} dPhi)
         = Tr(-Phi^{-1} dPhi Phi^{-1} dPhi + Phi^{-1} d^2Phi/dt^2)
         = -Tr((Phi^{-1} dPhi)^2) + Tr(Phi^{-1} d^2Phi/dt^2)

Now dPhi = V + alpha^T Phi + Phi alpha, so:
  Phi^{-1} dPhi = Phi^{-1} V + Phi^{-1} alpha^T Phi + alpha

And d^2Phi/dt^2 needs: Vdot + d/dt(alpha^T Phi + Phi alpha).

From W = Vdot - alpha^T V - V alpha (covariant visible second jet):
  Vdot = W + alpha^T V + V alpha

So d^2Phi/dt^2 involves W, alpha, dalpha/dt, Phi, dPhi.

Rather than chase the full algebra, verify numerically that f'' can be computed
from the split-frame quantities, and find the exact formula.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det
from scipy.linalg import sqrtm

TOL = 1e-6  # finite difference tolerance
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


def full_evidence_second_derivative(n, m, seed=0):
    """
    Compute f''(0) both numerically and from split-frame data.
    Derive the exact formula relating f'' to A_cpl.
    """
    print(f"\n=== Evidence second derivative (n={n}, m={m}) ===")
    rng = np.random.default_rng(seed)

    H0 = spd(n, seed)
    H1 = sym(rng.standard_normal((n, n)) * 2.0)
    H2 = sym(rng.standard_normal((n, n)))
    C = rng.standard_normal((m, n))

    # H(t) = H0 + t H1 + t^2 H2
    # At t=0: H = H0, Hdot = H1, Hddot = 2 H2

    # ---- Numerical f, f', f'' ----
    def log_det_phi(t):
        H_t = H0 + t * H1 + t**2 * H2
        ev = eigh(H_t)[0]
        if np.min(ev) < 0.01:
            return np.nan
        return np.log(det(inv(C @ inv(H_t) @ C.T)))

    dt = 1e-5
    f0 = log_det_phi(0)
    fp = log_det_phi(dt)
    fm = log_det_phi(-dt)
    f_prime_num = (fp - fm) / (2 * dt)
    f_pprime_num = (fp - 2 * f0 + fm) / dt**2

    # ---- Split-frame exact computation ----
    Hinv = inv(H0)
    Phi = inv(C @ Hinv @ C.T)
    Phi_inv = inv(Phi)
    L = Hinv @ C.T @ Phi
    U_, s_, Vt_ = svd(C)
    Z = Vt_[m:].T
    R = Z.T @ H0 @ Z
    Rinv = inv(R)

    # Connection forms and jets
    dHinv = -Hinv @ H1 @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi

    alpha_t = Phi_inv @ (L.T @ H0 @ dL)
    theta_t = Rinv @ (Z.T @ H0 @ dL)

    V = L.T @ H1 @ L
    B = L.T @ H1 @ Z

    # f'(0) exact
    f_prime_exact = np.trace(Phi_inv @ dPhi)
    report("f'(0) numerical vs exact", abs(f_prime_num - f_prime_exact))

    # ---- Derive f''(0) ----
    # f''(0) = -Tr((Phi^{-1} dPhi)^2) + Tr(Phi^{-1} d^2Phi/dt^2)
    #
    # Need d^2Phi/dt^2. From Phi(t) = (C H(t)^{-1} C^T)^{-1}:
    # dPhi/dt = -Phi (C dH^{-1} C^T) Phi
    # d^2Phi/dt^2 = -dPhi (C dH^{-1} C^T) Phi - Phi (C d^2H^{-1} C^T) Phi - Phi (C dH^{-1} C^T) dPhi

    # d^2(H^{-1})/dt^2 at t=0:
    # d(H^{-1})/dt = -H^{-1} Hdot H^{-1}
    # d^2(H^{-1})/dt^2 = 2 H^{-1} Hdot H^{-1} Hdot H^{-1} - H^{-1} Hddot H^{-1}
    Hddot = 2.0 * H2
    d2Hinv = 2.0 * Hinv @ H1 @ Hinv @ H1 @ Hinv - Hinv @ Hddot @ Hinv

    term_CdHinvC = C @ dHinv @ C.T
    term_Cd2HinvC = C @ d2Hinv @ C.T

    d2Phi = (-dPhi @ term_CdHinvC @ Phi
             - Phi @ term_Cd2HinvC @ Phi
             - Phi @ term_CdHinvC @ dPhi)

    # f''(0) from exact formula
    P = Phi_inv @ dPhi  # m x m
    f_pprime_exact = -np.trace(P @ P) + np.trace(Phi_inv @ d2Phi)
    report("f''(0) numerical vs exact", abs(f_pprime_num - f_pprime_exact))

    # ---- Now express f''(0) in terms of split-frame quantities ----
    # P = Phi^{-1} dPhi = Phi^{-1} V + Phi^{-1} alpha^T Phi + alpha
    #   = Phi^{-1} V + alpha^T + alpha  (wait, Phi^{-1} alpha^T Phi != alpha^T in general)

    # Let's define: S = Phi^{-1} V (the "source contribution" to the log-derivative)
    S = Phi_inv @ V

    # P = S + Phi^{-1} alpha^T Phi + alpha
    P_from_split = S + Phi_inv @ alpha_t.T @ Phi + alpha_t
    report("P = Phi^{-1}V + Phi^{-1}alpha^T Phi + alpha", norm(P - P_from_split))

    # f'(0) = Tr(P) = Tr(S) + Tr(Phi^{-1} alpha^T Phi) + Tr(alpha)
    # Tr(Phi^{-1} alpha^T Phi) = Tr(alpha^T) = Tr(alpha)  (cyclic)
    # So f'(0) = Tr(S) + 2 Tr(alpha). Confirmed.

    # f''(0) = -Tr(P^2) + Tr(Phi^{-1} d^2Phi)
    # The key question: can Tr(Phi^{-1} d^2Phi) be expressed via W and A_cpl?

    # From the calligraphic V identity differentiated:
    # d^2Phi = dV + d(alpha^T Phi + Phi alpha)
    #        = dV + dalpha^T Phi + alpha^T dPhi + dPhi alpha + Phi dalpha
    #        = (W + alpha^T V + V alpha) + dalpha^T Phi + alpha^T dPhi + dPhi alpha + Phi dalpha

    # This is getting complex. Let me compute Tr(Phi^{-1} d^2Phi) directly
    # and compare with a formula involving W.

    # Compute W
    Vdot = dL.T @ H1 @ L + L.T @ Hddot @ L + L.T @ H1 @ dL
    W = Vdot - alpha_t.T @ V - V @ alpha_t

    # Compute Tr(Phi^{-1} W)
    tr_phi_inv_W = np.trace(Phi_inv @ W)

    # Also need the connection terms
    # dalpha/dt at t=0... this requires second derivatives of L
    # which are complex. Let me instead verify the following conjecture:

    # f''(0) = Tr(Phi^{-1} W) + connection_correction

    # The connection correction involves alpha and its derivatives.
    # For a simple check: when alpha = 0 (a specific frame choice),
    # f''(0) should equal -Tr((Phi^{-1} V)^2) + Tr(Phi^{-1} (W + some correction))

    # Actually, let me take a different approach. Use:
    # d/dt [log det Phi] = Tr(Phi^{-1} dPhi) = Tr(Phi^{-1}(V + alpha^T Phi + Phi alpha))
    # = Tr(Phi^{-1} V) + 2 Tr(alpha) = Tr(S) + 2 Tr(alpha)

    # d^2/dt^2 [log det Phi] = d/dt[Tr(S)] + 2 d/dt[Tr(alpha)]

    # d/dt [Tr(S)] = d/dt [Tr(Phi^{-1} V)]
    # = Tr(-Phi^{-1} dPhi Phi^{-1} V + Phi^{-1} Vdot)
    # = -Tr(P S) + Tr(Phi^{-1} Vdot)
    # = -Tr(P S) + Tr(Phi^{-1}(W + alpha^T V + V alpha))
    # = -Tr(P S) + Tr(Phi^{-1} W) + Tr(Phi^{-1} alpha^T V) + Tr(Phi^{-1} V alpha)
    # = -Tr(P S) + Tr(Phi^{-1} W) + Tr(alpha^T V Phi^{-1}) + Tr(S alpha)
    # = -Tr(P S) + Tr(Phi^{-1} W) + Tr(S alpha^T)^T ... hmm
    # Let me just use Tr(AB) = Tr(BA):
    # Tr(Phi^{-1} alpha^T V) = Tr(V Phi^{-1} alpha^T) = Tr(alpha^T V Phi^{-1})
    # Hmm, this gets messy. Let me just verify numerically.

    # d/dt [Tr(Phi^{-1} V)] numerically:
    def phi_inv_V_trace(t):
        H_t = H0 + t*H1 + t**2*H2
        Hdot_t = H1 + 2*t*H2
        ev = eigh(H_t)[0]
        if np.min(ev) < 0.01:
            return np.nan
        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        L_t = Hinv_t @ C.T @ Phi_t
        V_t = L_t.T @ Hdot_t @ L_t
        return np.trace(inv(Phi_t) @ V_t)

    trS_0 = phi_inv_V_trace(0)
    trS_p = phi_inv_V_trace(dt)
    trS_m = phi_inv_V_trace(-dt)
    ddt_trS = (trS_p - trS_m) / (2*dt)

    # Check: d/dt[Tr(S)] + 2 d/dt[Tr(alpha)] = f''(0)
    # We need d/dt[Tr(alpha)] numerically too
    def alpha_trace(t):
        H_t = H0 + t*H1 + t**2*H2
        Hdot_t = H1 + 2*t*H2
        ev = eigh(H_t)[0]
        if np.min(ev) < 0.01:
            return np.nan
        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        L_t = Hinv_t @ C.T @ Phi_t
        dHinv_t = -Hinv_t @ Hdot_t @ Hinv_t
        dPhi_t = -Phi_t @ (C @ dHinv_t @ C.T) @ Phi_t
        dL_t = dHinv_t @ C.T @ Phi_t + Hinv_t @ C.T @ dPhi_t
        alpha_val = inv(Phi_t) @ (L_t.T @ H_t @ dL_t)
        return np.trace(alpha_val)

    trA_0 = alpha_trace(0)
    trA_p = alpha_trace(dt)
    trA_m = alpha_trace(-dt)
    ddt_trA = (trA_p - trA_m) / (2*dt)

    f_pprime_from_parts = ddt_trS + 2 * ddt_trA
    report("f'' = d/dt[Tr(S)] + 2 d/dt[Tr(alpha)]", abs(f_pprime_num - f_pprime_from_parts))

    # ---- Key result: express d/dt[Tr(S)] via W ----
    # d/dt[Tr(Phi^{-1} V)] = -Tr(P S) + Tr(Phi^{-1} W) + Tr(Phi^{-1} alpha^T V) + Tr(S alpha)
    # = -Tr(P S) + Tr(Phi^{-1} W) + 2 Re(Tr(S alpha))  [since Tr(Phi^{-1} alpha^T V) = Tr(S^T alpha^T)^T... ]

    # Actually: Tr(Phi^{-1} alpha^T V) = Tr(alpha^T V Phi^{-1}) = Tr(alpha S^T)^T
    # For real matrices: Tr(alpha S^T) is real.
    # And Tr(S alpha) is also Tr(alpha S).
    # So Tr(Phi^{-1} alpha^T V) + Tr(S alpha) = Tr(alpha S^T) + Tr(alpha S)
    # Hmm wait: Tr(Phi^{-1} alpha^T V) = Tr(V Phi^{-1} alpha^T) = Tr(S alpha^T)^T... let me just compute numerically.

    term_PS = np.trace(P @ S)
    term_PhiinvW = np.trace(Phi_inv @ W)
    term_cross1 = np.trace(Phi_inv @ alpha_t.T @ V)
    term_cross2 = np.trace(S @ alpha_t)

    ddt_trS_exact = -term_PS + term_PhiinvW + term_cross1 + term_cross2
    report("d/dt[Tr(S)] from W and alpha", abs(ddt_trS - ddt_trS_exact))

    # ---- Connection to A_cpl ----
    # When V > 0: A_cpl = -1/2 V^{-1/2} W V^{-1/2}
    # So W = -2 V^{1/2} A_cpl V^{1/2}
    # And Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})

    eigvals_V = eigh(V)[0]
    if np.min(eigvals_V) > 0.01:
        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        A_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)

        # Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})
        tr_from_Acpl = -2.0 * np.trace(Phi_inv @ Vsqrt @ A_cpl @ Vsqrt)
        report("Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})", abs(term_PhiinvW - tr_from_Acpl))

        # For the scalar case (m=1): Phi, V, A_cpl are all scalars.
        # Tr(Phi^{-1} W) = W/Phi = -2 V A_cpl / Phi
        # d/dt[Tr(S)] = -(V/Phi)^2 - 2VA_cpl/Phi + 2(V/Phi)*alpha + ... hmm still messy.

        # Clean summary formula:
        # f''(0) = d/dt[Tr(S)] + 2 d/dt[Tr(alpha)]
        # where d/dt[Tr(S)] = -Tr(P S) - 2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2}) + cross_terms
        # The A_cpl contribution is the genuinely new part.

        print(f"\n  Summary at (n={n}, m={m}):")
        print(f"    f'(0)  = {f_prime_exact:.6f}")
        print(f"    f''(0) = {f_pprime_exact:.6f}")
        print(f"    Tr(Phi^-1 W) = {term_PhiinvW:.6f}")
        print(f"    -2 Tr(Phi^-1 V^1/2 A_cpl V^1/2) = {tr_from_Acpl:.6f}")
        print(f"    A_cpl eigenvalues: {eigh(A_cpl)[0]}")
        print(f"    V eigenvalues: {eigh(V)[0]}")

        # The sign of Tr(Phi^{-1} W) is governed by A_cpl:
        # When A_cpl > 0 (positive definite), Tr(Phi^{-1} W) < 0 (decelerating)
        # When A_cpl < 0 (negative definite), Tr(Phi^{-1} W) > 0 (accelerating)
    else:
        print(f"  V not positive definite, skipping A_cpl analysis")

    return f_pprime_num, f_pprime_exact


def scalar_evidence_formula(seed=0):
    """
    For the scalar case (m=1), derive the complete evidence second derivative.

    When m=1, all quantities are scalars:
      Phi = phi (positive scalar)
      V = v = L^T Hdot L (scalar, sign depends on Hdot)
      W = w (scalar)
      alpha = a (scalar)
      S = v/phi

    f'(0) = v/phi + 2a
    f''(0) = d/dt[v/phi] + 2 d/dt[a]

    d/dt[v/phi] = (phi vdot - v dphi)/ phi^2
                = vdot/phi - v dphi/phi^2
                = (w + a v + v a)/phi - v(v + 2a phi)/phi^2
                = w/phi + 2av/phi - v^2/phi^2 - 2av/phi
                = w/phi - v^2/phi^2
                = w/phi - (v/phi)^2

    When v > 0: A_cpl = -w/(2v), so w = -2v A_cpl.
    d/dt[v/phi] = -2v A_cpl / phi - v^2/phi^2

    f''(0) = -2v A_cpl / phi - v^2/phi^2 + 2 da/dt

    The A_cpl term is the genuinely geometric contribution.
    The v^2/phi^2 term is the "kinematic" quadratic correction.
    The da/dt term is the connection acceleration.
    """
    print("\n=== Scalar evidence formula (m=1) ===")
    rng = np.random.default_rng(seed + 900)

    n = 4
    m = 1

    for attempt in range(50):
        H0 = spd(n, seed + attempt * 300)
        H1 = sym(rng.standard_normal((n, n)) * 3.0)
        H2 = sym(rng.standard_normal((n, n)))
        C = rng.standard_normal((m, n))

        Hinv = inv(H0)
        Phi = inv(C @ Hinv @ C.T)  # 1x1
        L = Hinv @ C.T @ Phi
        U_, s_, Vt_ = svd(C)
        Z = Vt_[m:].T
        R = Z.T @ H0 @ Z
        Rinv = inv(R)

        V = L.T @ H1 @ L  # 1x1
        if V.item() > 0.1:
            break
    else:
        print("  Skipped: could not find V > 0")
        return

    phi = Phi.item()
    v = V.item()

    dHinv = -Hinv @ H1 @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi
    alpha_t = inv(Phi) @ (L.T @ H0 @ dL)
    a = alpha_t.item()

    Hddot = 2.0 * H2
    B = L.T @ H1 @ Z
    Vdot = dL.T @ H1 @ L + L.T @ Hddot @ L + L.T @ H1 @ dL
    W = Vdot - alpha_t.T @ V - V @ alpha_t
    w = W.item()

    A_cpl_scalar = -0.5 * w / v

    # Numerical check
    dt = 1e-5
    def log_det_phi(t):
        H_t = H0 + t*H1 + t**2*H2
        return np.log(det(inv(C @ inv(H_t) @ C.T)))

    f0 = log_det_phi(0)
    fp = log_det_phi(dt)
    fm = log_det_phi(-dt)
    f_pprime_num = (fp - 2*f0 + fm) / dt**2

    # Formula: f''(0) = -2v A_cpl / phi - v^2/phi^2 + 2 da/dt
    # Need da/dt numerically
    def alpha_val(t):
        H_t = H0 + t*H1 + t**2*H2
        Hdot_t = H1 + 2*t*H2
        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        L_t = Hinv_t @ C.T @ Phi_t
        dHinv_t = -Hinv_t @ Hdot_t @ Hinv_t
        dPhi_t = -Phi_t @ (C @ dHinv_t @ C.T) @ Phi_t
        dL_t = dHinv_t @ C.T @ Phi_t + Hinv_t @ C.T @ dPhi_t
        return (inv(Phi_t) @ (L_t.T @ H_t @ dL_t)).item()

    da_dt = (alpha_val(dt) - alpha_val(-dt)) / (2*dt)

    f_pprime_formula = -2*v*A_cpl_scalar/phi - (v/phi)**2 + 2*da_dt
    report("Scalar f'' formula", abs(f_pprime_num - f_pprime_formula))

    # Also check: w/phi - (v/phi)^2 = d/dt[v/phi]
    ddt_S = w/phi - (v/phi)**2
    def S_val(t):
        H_t = H0 + t*H1 + t**2*H2
        Hdot_t = H1 + 2*t*H2
        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        L_t = Hinv_t @ C.T @ Phi_t
        V_t = L_t.T @ Hdot_t @ L_t
        return (inv(Phi_t) @ V_t).item()
    ddt_S_num = (S_val(dt) - S_val(-dt)) / (2*dt)
    report("d/dt[v/phi] = w/phi - (v/phi)^2", abs(ddt_S - ddt_S_num))

    print(f"\n  Scalar evidence decomposition:")
    print(f"    phi = {phi:.6f}")
    print(f"    v = {v:.6f}")
    print(f"    w = {w:.6f}")
    print(f"    A_cpl = {A_cpl_scalar:.6f}")
    print(f"    alpha = {a:.6f}")
    print(f"    da/dt = {da_dt:.6f}")
    print(f"    f''(0) = {f_pprime_num:.6f}")
    print(f"")
    print(f"    f''(0) = -2vA/phi  -  (v/phi)^2  +  2 da/dt")
    print(f"           = {-2*v*A_cpl_scalar/phi:.4f}   {-(v/phi)**2:.4f}   {2*da_dt:.4f}")
    print(f"             source     kinematic    connection")
    print(f"")
    print(f"    The source term (-2vA/phi) is the A_cpl contribution.")
    print(f"    When A_cpl > 0: source term < 0 (decelerates f, concave).")
    print(f"    When A_cpl < 0: source term > 0 (accelerates f, convex).")


if __name__ == "__main__":
    print("=" * 70)
    print("Goal A: Evidence Second Derivative from A_cpl")
    print("=" * 70)

    full_evidence_second_derivative(3, 1, seed=0)
    full_evidence_second_derivative(4, 2, seed=10)
    full_evidence_second_derivative(5, 3, seed=20)

    scalar_evidence_formula(seed=0)
    scalar_evidence_formula(seed=42)

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
