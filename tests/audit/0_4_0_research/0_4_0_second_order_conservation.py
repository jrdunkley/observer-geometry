"""
Second-order conservation law.

First order (proved in R2):
  d/dt[log det Phi] + d/dt[log det R] = d/dt[log det H]

Second order (to prove):
  d^2/dt^2[log det Phi] + d^2/dt^2[log det R] = d^2/dt^2[log det H]

The visible acceleration involves A_cpl.
The hidden acceleration involves W_h = Z^T Hddot Z (direct, no observer structure).
The ambient acceleration is Tr(H^{-1} Hddot) - Tr((H^{-1} Hdot)^2).

CONSEQUENCE: A_cpl is constrained by the ambient geometry.
The observer can redistribute information between visible and hidden,
but the total acceleration is fixed. This is the information-conservation
constraint on observer design.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det, slogdet
from scipy.linalg import sqrtm

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


def second_order_conservation(n, m, seed=0):
    """
    Verify: f''_vis + f''_hid = f''_amb

    where:
    f''_vis = d^2/dt^2 [log det Phi(t)]
    f''_hid = d^2/dt^2 [log det R(t)]
    f''_amb = d^2/dt^2 [log det H(t)]

    All evaluated at t=0 for H(t) = H0 + t*H1 + t^2*H2, fixed C.
    """
    print(f"\n=== Second-order conservation (n={n}, m={m}) ===")
    rng = np.random.default_rng(seed)

    H0 = spd(n, seed)
    H1 = sym(rng.standard_normal((n, n)) * 2.0)
    H2 = sym(rng.standard_normal((n, n)))
    C = rng.standard_normal((m, n))
    _, _, Vt = svd(C)
    Z = Vt[m:].T

    # Numerical computation via finite differences
    dt = 1e-5

    def log_dets(t):
        H_t = H0 + t * H1 + t**2 * H2
        ev = eigh(H_t)[0]
        if np.min(ev) < 0.01:
            return np.nan, np.nan, np.nan
        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        R_t = Z.T @ H_t @ Z
        _, ld_phi = slogdet(Phi_t)
        _, ld_R = slogdet(R_t)
        _, ld_H = slogdet(H_t)
        return ld_phi, ld_R, ld_H

    ld_phi_0, ld_R_0, ld_H_0 = log_dets(0)
    ld_phi_p, ld_R_p, ld_H_p = log_dets(dt)
    ld_phi_m, ld_R_m, ld_H_m = log_dets(-dt)

    # Second derivatives
    f2_vis = (ld_phi_p - 2*ld_phi_0 + ld_phi_m) / dt**2
    f2_hid = (ld_R_p - 2*ld_R_0 + ld_R_m) / dt**2
    f2_amb = (ld_H_p - 2*ld_H_0 + ld_H_m) / dt**2

    # Conservation check
    conservation_2nd = f2_vis + f2_hid - f2_amb
    report(f"f''_vis + f''_hid = f''_amb", abs(conservation_2nd), tol=1e-3)

    # Analytical computation of each term
    Hinv = inv(H0)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    R = Z.T @ H0 @ Z
    Rinv = inv(R)
    Hdot = H1
    Hddot = 2.0 * H2

    # --- Ambient second derivative ---
    # d^2/dt^2 [log det H] = -Tr((H^{-1} Hdot)^2) + Tr(H^{-1} Hddot)
    P_H = Hinv @ Hdot
    f2_amb_exact = -np.trace(P_H @ P_H) + np.trace(Hinv @ Hddot)
    report("f''_amb exact vs numerical", abs(f2_amb - f2_amb_exact), tol=1e-3)

    # --- Hidden second derivative ---
    # d^2/dt^2 [log det R] where R(t) = Z^T H(t) Z, dR/dt = Z^T Hdot Z = U_h
    U_h = Z.T @ Hdot @ Z
    U_h2 = Z.T @ Hddot @ Z  # d^2R/dt^2
    P_R = Rinv @ U_h
    f2_hid_exact = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)
    report("f''_hid exact vs numerical", abs(f2_hid - f2_hid_exact), tol=1e-3)

    # --- Visible second derivative ---
    # From the conservation: f2_vis = f2_amb - f2_hid
    f2_vis_from_conservation = f2_amb_exact - f2_hid_exact
    report("f''_vis from conservation", abs(f2_vis - f2_vis_from_conservation), tol=1e-3)

    # Also compute f2_vis directly (involves A_cpl)
    V = L.T @ Hdot @ L
    dHinv = -Hinv @ Hdot @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    P_Phi = inv(Phi) @ dPhi
    d2Hinv = 2.0 * Hinv @ Hdot @ Hinv @ Hdot @ Hinv - Hinv @ Hddot @ Hinv
    d2Phi = (-dPhi @ (C @ dHinv @ C.T) @ Phi
             - Phi @ (C @ d2Hinv @ C.T) @ Phi
             - Phi @ (C @ dHinv @ C.T) @ dPhi)
    f2_vis_exact = -np.trace(P_Phi @ P_Phi) + np.trace(inv(Phi) @ d2Phi)
    report("f''_vis exact vs numerical", abs(f2_vis - f2_vis_exact), tol=1e-3)

    # --- THE KEY RESULT: decomposition ---
    print(f"\n  Second-order conservation at (n={n}, m={m}):")
    print(f"    f''_vis (visible acceleration) = {f2_vis_exact:.6f}")
    print(f"    f''_hid (hidden acceleration)  = {f2_hid_exact:.6f}")
    print(f"    f''_amb (ambient acceleration)  = {f2_amb_exact:.6f}")
    print(f"    f''_vis + f''_hid              = {f2_vis_exact + f2_hid_exact:.6f}")
    print(f"    Conservation residual:          {abs(f2_vis_exact + f2_hid_exact - f2_amb_exact):.3e}")

    # Fractional split
    if abs(f2_amb_exact) > 1e-10:
        vis_frac = f2_vis_exact / f2_amb_exact
        hid_frac = f2_hid_exact / f2_amb_exact
        print(f"\n    Visible fraction of acceleration: {vis_frac:.4f}")
        print(f"    Hidden fraction of acceleration:  {hid_frac:.4f}")
        print(f"    Sum:                              {vis_frac + hid_frac:.4f}")

    # --- Connection to A_cpl ---
    eigvals_V = eigh(V)[0]
    if np.min(eigvals_V) > 0.01:
        dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi
        alpha_t = inv(Phi) @ (L.T @ H0 @ dL)
        Vdot = dL.T @ Hdot @ L + L.T @ Hddot @ L + L.T @ Hdot @ dL
        W = Vdot - alpha_t.T @ V - V @ alpha_t

        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        A_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)

        # The A_cpl contribution to f''_vis:
        # Tr(Phi^{-1} W) = -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})
        Acpl_contribution = -2.0 * np.trace(inv(Phi) @ Vsqrt @ A_cpl @ Vsqrt)

        # The kinematic contribution: -Tr(P_Phi^2) + other terms
        # f''_vis = -Tr(P_Phi^2) + Tr(Phi^{-1} d^2Phi)
        kinematic = -np.trace(P_Phi @ P_Phi)

        print(f"\n    A_cpl eigenvalues: {sorted(eigh(A_cpl)[0])}")
        print(f"    A_cpl contribution to f''_vis: {Acpl_contribution:.6f}")
        print(f"    Kinematic contribution:        {kinematic:.6f}")
        print(f"    Remainder:                     {f2_vis_exact - Acpl_contribution - kinematic:.6f}")

    return f2_vis_exact, f2_hid_exact, f2_amb_exact


def conservation_constraint_on_acpl():
    """
    THEOREM (A_cpl constraint from ambient geometry):

    Since f''_vis + f''_hid = f''_amb, and f''_hid depends only on (H, Z)
    (not on C), the visible acceleration f''_vis is CONSTRAINED:

      f''_vis = f''_amb - f''_hid

    The A_cpl contribution to f''_vis is:
      -2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2})

    This means A_cpl cannot be chosen freely. Given the ambient geometry (H, Hdot, Hddot)
    and the hidden frame (Z), the visible acceleration is determined.

    The only freedom is in C (the observer choice), which affects:
    1. The split between visible and hidden (Phi vs R)
    2. The visible jet V = L^T Hdot L
    3. The connection forms alpha, theta
    4. A_cpl itself

    But all of these are constrained by the conservation law. An observer that
    maximises the visible acceleration must minimise the hidden acceleration,
    and vice versa — within the constraint imposed by the ambient geometry.

    This is the information-conservation principle for observer design.
    """
    print("\n=== A_cpl constraint from ambient geometry ===")

    # Demonstrate with coupled spring: varying the observer C while keeping H fixed
    n = 3
    H = np.array([
        [4.0, 1.0, 0.5],
        [1.0, 3.0, 0.2],
        [0.5, 0.2, 2.0],
    ])
    Hdot = np.array([
        [1.0, 0.3, 0.1],
        [0.3, 0.5, -0.1],
        [0.1, -0.1, 0.2],
    ])
    Hddot = np.array([
        [0.2, 0.1, 0.0],
        [0.1, 0.3, 0.0],
        [0.0, 0.0, 0.1],
    ])

    Hinv = inv(H)

    # Ambient acceleration (FIXED, independent of observer)
    P_H = Hinv @ Hdot
    f2_amb = -np.trace(P_H @ P_H) + np.trace(Hinv @ Hddot)

    print(f"  Ambient acceleration (FIXED): {f2_amb:.6f}")
    print(f"  Now varying the observer C (m=1) while H is fixed...\n")

    # Try different observers
    rng = np.random.default_rng(42)
    observers = []
    for i in range(20):
        C = rng.standard_normal((1, n))
        C = C / norm(C)  # normalise

        _, _, Vt = svd(C)
        Z = Vt[1:].T

        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        R = Z.T @ H @ Z

        # Visible acceleration
        dHinv = -Hinv @ Hdot @ Hinv
        dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
        P_Phi = inv(Phi) @ dPhi
        d2Hinv = 2 * Hinv @ Hdot @ Hinv @ Hdot @ Hinv - Hinv @ Hddot @ Hinv
        d2Phi = (-dPhi @ (C @ dHinv @ C.T) @ Phi
                 - Phi @ (C @ d2Hinv @ C.T) @ Phi
                 - Phi @ (C @ dHinv @ C.T) @ dPhi)
        f2_vis = (-np.trace(P_Phi @ P_Phi) + np.trace(inv(Phi) @ d2Phi))

        # Hidden acceleration
        U_h = Z.T @ Hdot @ Z
        U_h2 = Z.T @ Hddot @ Z
        Rinv = inv(R)
        P_R = Rinv @ U_h
        f2_hid = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)

        # Conservation check
        residual = abs(f2_vis + f2_hid - f2_amb)

        observers.append({
            'C': C.flatten(),
            'phi': Phi.item(),
            'f2_vis': f2_vis,
            'f2_hid': f2_hid,
            'residual': residual,
        })

    # Sort by visible acceleration
    observers.sort(key=lambda x: x['f2_vis'])

    print(f"  {'Phi':>8s}  {'f2_vis':>10s}  {'f2_hid':>10s}  {'sum':>10s}  {'residual':>10s}")
    for o in observers:
        print(f"  {o['phi']:8.4f}  {o['f2_vis']:10.4f}  {o['f2_hid']:10.4f}  "
              f"{o['f2_vis']+o['f2_hid']:10.4f}  {o['residual']:10.2e}")

    max_residual = max(o['residual'] for o in observers)
    report("Conservation holds for all observers", max_residual, tol=1e-6)

    # The KEY observation:
    vis_range = max(o['f2_vis'] for o in observers) - min(o['f2_vis'] for o in observers)
    print(f"\n  Visible acceleration range: {vis_range:.4f}")
    print(f"  Ambient acceleration (fixed): {f2_amb:.4f}")
    print(f"\n  The observer choice redistributes acceleration between vis and hid,")
    print(f"  but the total is ALWAYS {f2_amb:.4f}. This is the conservation constraint.")
    print(f"  The observer is not a source of information — only a splitter.")


if __name__ == "__main__":
    print("=" * 70)
    print("Second-Order Conservation Law")
    print("=" * 70)

    second_order_conservation(3, 1, seed=0)
    second_order_conservation(4, 2, seed=10)
    second_order_conservation(5, 3, seed=20)
    second_order_conservation(6, 2, seed=30)

    conservation_constraint_on_acpl()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
