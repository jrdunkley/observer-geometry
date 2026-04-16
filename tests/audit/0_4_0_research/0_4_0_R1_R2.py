"""
R1: Information-destruction locus.
R2: Split determinant conservation law.

R1 maps the codimension-1 surface in parameter space where visible information
is destroyed (V drops rank). R2 derives and verifies the pathwise determinant
decomposition into visible + hidden + coupling contributions.
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


# ============================================================
# R1: Information-destruction locus
# ============================================================

def information_destruction_coupled_spring():
    """
    Coupled spring: H(kc) = [[k1+kc, -kc], [-kc, k2+kc]], C = [1,0].

    V = L^T dH L where dH = [[1,-1],[-1,1]] (d/dkc).

    For n=2, m=1: V is scalar. V = 0 is the information-destruction event.

    ANALYTICAL:
    L = H^{-1} C^T Phi = [1, 0]^T (constant for this diagonal observer).
    Wait: L = H^{-1} [1;0] * Phi.
    H^{-1} = (1/det(H)) * [[k2+kc, kc], [kc, k1+kc]]
    H^{-1} [1;0] = (1/det(H)) * [k2+kc, kc]^T
    Phi = det(H)/(k2+kc) = (k1*k2 + k1*kc + k2*kc)/(k2+kc)
    L = [1, kc/(k2+kc)]^T ... wait.

    H^{-1} C^T = H^{-1} [1;0] = [H^{-1}_{11}, H^{-1}_{21}]^T
    = [(k2+kc)/det(H), kc/det(H)]^T

    Phi = 1/H^{-1}_{11} = det(H)/(k2+kc)
    L = H^{-1} C^T Phi = [(k2+kc)/det(H), kc/det(H)]^T * det(H)/(k2+kc)
    = [1, kc/(k2+kc)]^T

    V = L^T dH L = [1, kc/(k2+kc)] [[1,-1],[-1,1]] [1, kc/(k2+kc)]^T
    = [1, kc/(k2+kc)] [1 - kc/(k2+kc), -1 + kc/(k2+kc)]^T
    = [1, kc/(k2+kc)] [k2/(k2+kc), -(k2)/(k2+kc)]^T
    = k2/(k2+kc) - kc*k2/(k2+kc)^2
    = k2(k2+kc-kc)/(k2+kc)^2
    = k2^2/(k2+kc)^2

    V = k2^2/(k2+kc)^2 > 0 for all kc > 0 and k2 > 0.

    The information-destruction locus is EMPTY for the coupled spring system!
    V never reaches zero — the observer always receives information about the
    coupling change. This is because the observer (mass 1) is directly coupled
    to mass 2 through the coupling spring.

    V -> 0 only as kc -> infinity, which is the limit where the system becomes
    rigid and there's no internal dynamics to observe.
    """
    print("\n=== R1: Information-destruction locus (coupled spring) ===")

    k1, k2 = 3.0, 2.0
    C = np.array([[1.0, 0.0]])
    dH = np.array([[1.0, -1.0], [-1.0, 1.0]])

    kcs = np.logspace(-2, 3, 200)
    vs = []
    for kc in kcs:
        H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        V = (L.T @ dH @ L).item()
        vs.append(V)

        # Verify analytical formula
        V_analytical = k2**2 / (k2 + kc)**2
        assert abs(V - V_analytical) < 1e-10

    vs = np.array(vs)
    print(f"  V = k2^2 / (k2+kc)^2 > 0 always (k2={k2})")
    print(f"  V range: [{vs.min():.6e}, {vs.max():.6e}]")
    print(f"  Information-destruction locus: EMPTY")
    print(f"  V -> 0 only as kc -> infinity (rigid limit)")
    print(f"  Observer always sees information about coupling change.")


def information_destruction_entanglement():
    """
    Two-mode squeezed thermal state:
    H(r) = inv(sigma(r)), C = [[1,0,0,0],[0,1,0,0]].

    The parameter is squeezing r. dH/dr at various r values
    determines V = L^T dH L (2x2 matrix). V drops rank when
    one eigenvalue crosses zero.

    Use finite differences for dH/dr.
    """
    print("\n=== R1: Information-destruction locus (entanglement) ===")

    def sigma(r, n_a=0.0, n_b=0.0):
        c = np.cosh(r); s = np.sinh(r)
        return np.array([
            [(2*n_a+1)*c**2 + (2*n_b+1)*s**2, 0, (2*n_a+2*n_b+2)*c*s, 0],
            [0, (2*n_a+1)*c**2 + (2*n_b+1)*s**2, 0, -(2*n_a+2*n_b+2)*c*s],
            [(2*n_a+2*n_b+2)*c*s, 0, (2*n_b+1)*c**2 + (2*n_a+1)*s**2, 0],
            [0, -(2*n_a+2*n_b+2)*c*s, 0, (2*n_b+1)*c**2 + (2*n_a+1)*s**2],
        ])

    C = np.array([[1,0,0,0],[0,1,0,0]], dtype=float)
    n, m = 4, 2

    rs = np.linspace(0.05, 2.0, 100)
    dr = 1e-6
    v_min_eigvals = []
    v_max_eigvals = []
    phi_min_eigvals = []

    for r in rs:
        H = inv(sigma(r))
        H_plus = inv(sigma(r + dr))
        H_minus = inv(sigma(r - dr))
        dH = (H_plus - H_minus) / (2 * dr)

        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi

        V = L.T @ dH @ L
        eigvals_V = sorted(eigh(V)[0])
        v_min_eigvals.append(eigvals_V[0])
        v_max_eigvals.append(eigvals_V[-1])
        phi_min_eigvals.append(min(eigh(Phi)[0]))

    v_min = np.array(v_min_eigvals)
    v_max = np.array(v_max_eigvals)

    # Find where V eigenvalues cross zero
    sign_changes_min = []
    sign_changes_max = []
    for i in range(1, len(rs)):
        if v_min[i] * v_min[i-1] < 0:
            sign_changes_min.append(rs[i])
        if v_max[i] * v_max[i-1] < 0:
            sign_changes_max.append(rs[i])

    print(f"  Parameter: squeezing r in [{rs[0]:.2f}, {rs[-1]:.2f}]")
    print(f"  V eigenvalue 1 range: [{v_min.min():.4f}, {v_min.max():.4f}]")
    print(f"  V eigenvalue 2 range: [{v_max.min():.4f}, {v_max.max():.4f}]")
    print(f"  V eigenvalue 1 zero-crossings: {sign_changes_min}")
    print(f"  V eigenvalue 2 zero-crossings: {sign_changes_max}")

    if len(sign_changes_min) > 0 or len(sign_changes_max) > 0:
        print(f"  INFORMATION-DESTRUCTION LOCUS FOUND!")
        for sc in sign_changes_min + sign_changes_max:
            print(f"    r = {sc:.4f}: a visible direction goes dark")
    else:
        # Check: are both eigenvalues always the same sign?
        if np.all(v_min > 0) and np.all(v_max > 0):
            print(f"  Both V eigenvalues positive throughout: no information destruction")
        elif np.all(v_min < 0) and np.all(v_max < 0):
            print(f"  Both V eigenvalues negative: Phi DECREASING with squeezing")
        else:
            print(f"  V has mixed sign eigenvalues: one direction sees information, other loses it")

    # Report the eigenvalue structure
    print(f"\n  V eigenvalue profile:")
    print(f"  {'r':>6s}  {'V_min':>10s}  {'V_max':>10s}  {'Phi_min':>10s}")
    for i in [0, 10, 25, 50, 75, 99]:
        print(f"  {rs[i]:6.3f}  {v_min[i]:10.4f}  {v_max[i]:10.4f}  {phi_min_eigvals[i]:10.4f}")

    print(f"\n  Interpretation:")
    if np.all(v_min < 0):
        print(f"  V_min < 0 throughout: increasing squeezing DECREASES precision in")
        print(f"  at least one visible direction. This is the price of entanglement:")
        print(f"  to entangle the modes, local precision must be sacrificed.")
    if np.all(v_max < 0):
        print(f"  V_max < 0 too: ALL visible precision decreases with squeezing.")
        print(f"  The information-destruction locus is at r=0 (no squeezing).")
        print(f"  Beyond that point, every increase in squeezing destroys visible information.")


# ============================================================
# R2: Split determinant conservation law
# ============================================================

def split_determinant_conservation(n, m, seed=0):
    """
    THEOREM (Split determinant identity):
    For M = [L Z] with M^T H M = diag(Phi, R):

      det(H) = det(Phi) * det(R) / det(M^T M)   ... wait, this isn't right.

    Actually: M^T H M = diag(Phi, R) implies:
      det(M)^2 det(H) = det(M^T H M) = det(Phi) det(R)

    But M is n x n (not necessarily orthogonal), so det(M) != 1 in general.

    The correct identity is:
      det(Phi) * det(R) = det(H) * det(M)^2

    Taking logs and differentiating:
      d/dt [log det Phi] + d/dt [log det R] = d/dt [log det H] + 2 d/dt [log det M]

    The LHS splits into visible and hidden evidence contributions.
    The RHS has the ambient contribution plus the frame-transport contribution.

    Rearranging:
      f_vis(t) + f_hid(t) = f_amb(t) + 2 f_frame(t)

    where f_vis = log det Phi, f_hid = log det R, f_amb = log det H, f_frame = log det M.

    This IS the conservation law: the total evidence (visible + hidden)
    equals the ambient evidence plus the frame contribution.
    """
    print(f"\n=== R2: Split determinant conservation (n={n}, m={m}) ===")
    rng = np.random.default_rng(seed)

    H0 = spd(n, seed)
    H1 = sym(rng.standard_normal((n, n)) * 2.0)
    H2 = sym(rng.standard_normal((n, n)))
    C = rng.standard_normal((m, n))

    # Verify the static identity: det(Phi) det(R) = det(H) det(M)^2
    Hinv = inv(H0)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C); Z = Vt[m:].T
    R = Z.T @ H0 @ Z
    M_frame = np.hstack([L, Z])

    lhs = det(Phi) * det(R)
    rhs = det(H0) * det(M_frame)**2
    report(f"det(Phi)*det(R) = det(H)*det(M)^2 (static)", abs(lhs - rhs) / max(abs(lhs), 1e-15))

    # Now verify the pathwise version
    ts = np.linspace(-0.15, 0.15, 201)
    log_det_phi = []
    log_det_R = []
    log_det_H = []
    log_det_M = []

    for t in ts:
        H_t = H0 + t * H1 + t**2 * H2
        ev = eigh(H_t)[0]
        if np.min(ev) < 0.01:
            log_det_phi.append(np.nan)
            log_det_R.append(np.nan)
            log_det_H.append(np.nan)
            log_det_M.append(np.nan)
            continue

        Hinv_t = inv(H_t)
        Phi_t = inv(C @ Hinv_t @ C.T)
        L_t = Hinv_t @ C.T @ Phi_t
        R_t = Z.T @ H_t @ Z  # Z is fixed (from SVD of fixed C)
        M_t = np.hstack([L_t, Z])

        _, ld_phi = slogdet(Phi_t)
        _, ld_R = slogdet(R_t)
        _, ld_H = slogdet(H_t)
        _, ld_M = slogdet(M_t)

        log_det_phi.append(ld_phi)
        log_det_R.append(ld_R)
        log_det_H.append(ld_H)
        log_det_M.append(ld_M)

    log_det_phi = np.array(log_det_phi)
    log_det_R = np.array(log_det_R)
    log_det_H = np.array(log_det_H)
    log_det_M = np.array(log_det_M)

    # Verify conservation: log det Phi + log det R = log det H + 2 log det M
    conservation = log_det_phi + log_det_R - log_det_H - 2 * log_det_M
    valid = ~np.isnan(conservation)
    max_violation = np.max(np.abs(conservation[valid]))
    report("Conservation along path", max_violation)

    # Derivatives at t=0
    dt = ts[1] - ts[0]
    mid = len(ts) // 2

    def deriv(arr):
        return (arr[mid+1] - arr[mid-1]) / (2*dt)

    d_log_phi = deriv(log_det_phi)
    d_log_R = deriv(log_det_R)
    d_log_H = deriv(log_det_H)
    d_log_M = deriv(log_det_M)

    conservation_deriv = d_log_phi + d_log_R - d_log_H - 2 * d_log_M
    report("Conservation derivative at t=0", abs(conservation_deriv), tol=1e-4)

    print(f"\n  Split determinant at t=0:")
    print(f"    d/dt [log det Phi] = {d_log_phi:.6f}  (visible evidence rate)")
    print(f"    d/dt [log det R]   = {d_log_R:.6f}  (hidden evidence rate)")
    print(f"    d/dt [log det H]   = {d_log_H:.6f}  (ambient rate)")
    print(f"    2 d/dt [log det M] = {2*d_log_M:.6f}  (frame transport rate)")
    print(f"    Conservation check: {conservation_deriv:.6e}")

    # Connection to source law: d/dt [log det Phi] = Tr(Phi^{-1} dPhi) = f'(0)
    # which we already know = Tr(Phi^{-1} V) + 2 Tr(alpha)
    # And d/dt [log det R] = Tr(R^{-1} dR) involves the hidden jet U_h = Z^T Hdot Z
    # and the hidden connection omega.

    # The conservation law says:
    #   Tr(Phi^{-1} V) + 2Tr(alpha) + Tr(R^{-1} U_h) + 2Tr(omega) = Tr(H^{-1} Hdot) + 2 d/dt[log det M]
    # Since Tr(H^{-1} Hdot) = d/dt[log det H] (the ambient rate).

    # Verify this decomposition
    V_mat = L.T @ H1 @ L
    U_h = Z.T @ H1 @ Z

    dHinv = -Hinv @ H1 @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi
    alpha_t = inv(Phi) @ (L.T @ H0 @ dL)
    omega_t = inv(R) @ (Z.T @ H0 @ np.zeros_like(Z))  # dZ = 0

    vis_rate = np.trace(inv(Phi) @ V_mat) + 2 * np.trace(alpha_t)
    hid_rate = np.trace(inv(R) @ U_h) + 2 * np.trace(omega_t)
    amb_rate = np.trace(Hinv @ H1)

    report("Visible rate matches d/dt[log det Phi]", abs(vis_rate - d_log_phi), tol=1e-4)
    report("Ambient rate matches d/dt[log det H]", abs(amb_rate - d_log_H), tol=1e-4)

    print(f"\n  Rate decomposition:")
    print(f"    Visible rate (Phi): {vis_rate:.6f}")
    print(f"    Hidden rate (R):    {hid_rate:.6f}")
    print(f"    Sum:                {vis_rate + hid_rate:.6f}")
    print(f"    Ambient rate:       {amb_rate:.6f}")
    print(f"    Frame rate:         {2*d_log_M:.6f}")
    print(f"    Ambient + frame:    {amb_rate + 2*d_log_M:.6f}")

    # The "coupling" is the frame transport contribution
    coupling_rate = 2 * d_log_M
    print(f"\n  CONSERVATION LAW:")
    print(f"    (visible rate) + (hidden rate) = (ambient rate) + (frame transport rate)")
    print(f"    {vis_rate:.4f} + {hid_rate:.4f} = {amb_rate:.4f} + {coupling_rate:.4f}")
    print(f"    {vis_rate + hid_rate:.4f} = {amb_rate + coupling_rate:.4f}")

    return vis_rate, hid_rate, amb_rate, coupling_rate


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("R1 + R2: Information-Destruction Locus + Determinant Conservation")
    print("=" * 70)

    # R1
    information_destruction_coupled_spring()
    information_destruction_entanglement()

    # R2
    split_determinant_conservation(3, 1, seed=0)
    split_determinant_conservation(4, 2, seed=10)
    split_determinant_conservation(5, 3, seed=20)

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
