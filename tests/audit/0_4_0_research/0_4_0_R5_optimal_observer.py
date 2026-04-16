"""
R5: Optimal observer alignment with ambient geometry.

The second-order conservation constrains f''_vis:
  f''_vis = f''_amb - f''_hid

The observer C determines f''_hid (through Z = ker C and R = Z^T H Z).
Maximising f''_vis means minimising f''_hid.

For m=1 (scalar observer), the observer is a unit vector c in R^n.
Z = ker(c) is the (n-1)-dimensional orthogonal complement.
f''_hid depends on Z through R = Z^T H Z and U_h = Z^T Hdot Z.

The optimal observer is the c that minimises:
  f''_hid(c) = -Tr((R^{-1} U_h)^2) + Tr(R^{-1} U_h2)

where R = Z^T H Z, U_h = Z^T Hdot Z, U_h2 = Z^T Hddot Z, and Z = ker(c^T).

This is a smooth function on the projective space RP^{n-1} (or equivalently
on the Grassmannian Gr(1,n)). Its critical points are the observers that
are stationary for the evidence acceleration.

QUESTION: Do these critical points coincide with the weighted-family frontier
observers from TN1? They should, at least in a specific limit.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det, slogdet
from scipy.linalg import sqrtm
from scipy.optimize import minimize

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok

def sym(M):
    return 0.5 * (M + M.T)


def optimal_observer_scalar(n, seed=0):
    """
    For m=1, find the observer c that maximises f''_vis.

    Parameterise: c = (cos theta_1, sin theta_1 cos theta_2, ...) on S^{n-1}.
    Use scipy.optimize on the sphere.
    """
    print(f"\n=== R5: Optimal observer (n={n}, m=1) ===")
    rng = np.random.default_rng(seed)

    H = np.array([
        [4.0, 1.0, 0.5],
        [1.0, 3.0, 0.2],
        [0.5, 0.2, 2.0],
    ]) if n == 3 else (rng.standard_normal((n,n)) @ rng.standard_normal((n,n)).T + np.eye(n))

    Hdot = sym(rng.standard_normal((n, n)))
    Hddot = sym(rng.standard_normal((n, n)))
    Hinv = inv(H)

    # Ambient acceleration (fixed)
    P_H = Hinv @ Hdot
    f2_amb = -np.trace(P_H @ P_H) + np.trace(Hinv @ Hddot)

    def f2_vis_from_c(c_vec):
        """Compute f''_vis for observer c."""
        c = c_vec.reshape(1, n)
        _, _, Vt = svd(c)
        Z = Vt[1:].T

        Phi = inv(c @ Hinv @ c.T)
        R = Z.T @ H @ Z
        Rinv = inv(R)
        U_h = Z.T @ Hdot @ Z
        U_h2 = Z.T @ Hddot @ Z
        P_R = Rinv @ U_h

        f2_hid = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)
        return f2_amb - f2_hid  # = f''_vis

    def neg_f2_vis(angles):
        """Negative f''_vis as function of spherical angles."""
        c = np.zeros(n)
        c[0] = np.cos(angles[0])
        for i in range(1, n - 1):
            c[i] = np.prod([np.sin(angles[j]) for j in range(i)]) * np.cos(angles[i])
        c[n-1] = np.prod([np.sin(angles[j]) for j in range(n-1)])
        return -f2_vis_from_c(c)

    # Grid search + refinement
    best_val = float('inf')
    best_angles = None
    best_c = None

    for trial in range(200):
        angles0 = rng.uniform(0, np.pi, n - 1)
        try:
            result = minimize(neg_f2_vis, angles0, method='Nelder-Mead',
                            options={'maxiter': 1000, 'xatol': 1e-10, 'fatol': 1e-12})
            if result.fun < best_val:
                best_val = result.fun
                best_angles = result.x
                # Reconstruct c
                c = np.zeros(n)
                c[0] = np.cos(best_angles[0])
                for i in range(1, n - 1):
                    c[i] = np.prod([np.sin(best_angles[j]) for j in range(i)]) * np.cos(best_angles[i])
                c[n-1] = np.prod([np.sin(best_angles[j]) for j in range(n-1)])
                best_c = c
        except:
            pass

    # Search for minimum f''_vis (= maximise f''_hid = maximise neg_f2_vis)
    worst_val = float('inf')
    worst_c = None
    for trial in range(200):
        angles0 = rng.uniform(0, np.pi, n - 1)
        try:
            # Minimise f''_vis = maximise -f''_vis = minimise neg(neg_f2_vis) ...
            # Actually: neg_f2_vis = -f2_vis. To minimise f2_vis, maximise neg_f2_vis.
            # So minimise -neg_f2_vis = minimise f2_vis.
            result = minimize(lambda a: -neg_f2_vis(a), angles0, method='Nelder-Mead',
                            options={'maxiter': 1000, 'xatol': 1e-10, 'fatol': 1e-12})
            if result.fun < worst_val:
                worst_val = result.fun
                angles = result.x
                c = np.zeros(n)
                c[0] = np.cos(angles[0])
                for i in range(1, n - 1):
                    c[i] = np.prod([np.sin(angles[j]) for j in range(i)]) * np.cos(angles[i])
                c[n-1] = np.prod([np.sin(angles[j]) for j in range(n-1)])
                worst_c = c
        except:
            pass

    max_f2_vis = -best_val
    min_f2_vis = worst_val

    print(f"  Ambient acceleration (FIXED): {f2_amb:.6f}")
    print(f"  Max f''_vis: {max_f2_vis:.6f}  (optimal observer)")
    print(f"  Min f''_vis: {min_f2_vis:.6f}  (worst observer)")
    print(f"  Range: {max_f2_vis - min_f2_vis:.6f}")

    if best_c is not None:
        print(f"\n  Optimal observer direction: c = [{', '.join(f'{x:.4f}' for x in best_c)}]")
        phi_opt = (best_c.reshape(1,n) @ Hinv @ best_c.reshape(n,1)).item()**(-1)
        print(f"  Visible precision at optimum: Phi = {phi_opt:.4f}")

    if worst_c is not None:
        print(f"  Worst observer direction:    c = [{', '.join(f'{x:.4f}' for x in worst_c)}]")
        phi_worst = (worst_c.reshape(1,n) @ Hinv @ worst_c.reshape(n,1)).item()**(-1)
        print(f"  Visible precision at worst:  Phi = {phi_worst:.4f}")

    # Check: does the optimal observer align with an eigenvector of H or Hdot?
    eigvals_H, eigvecs_H = eigh(H)
    eigvals_Hdot, eigvecs_Hdot = eigh(Hdot)

    if best_c is not None:
        print(f"\n  Alignment with H eigenvectors:")
        for i in range(n):
            alignment = abs(np.dot(best_c, eigvecs_H[:, i]))
            print(f"    eigval {eigvals_H[i]:.2f}: |<c*, e_i>| = {alignment:.4f}")

        print(f"  Alignment with Hdot eigenvectors:")
        for i in range(n):
            alignment = abs(np.dot(best_c, eigvecs_Hdot[:, i]))
            print(f"    eigval {eigvals_Hdot[i]:.2f}: |<c*, e_i>| = {alignment:.4f}")

    return max_f2_vis, min_f2_vis, best_c


def observer_landscape(seed=0):
    """
    For n=3, m=1, plot the full landscape of f''_vis on S^2.
    Identify all critical points and classify them.
    """
    print("\n=== R5: Observer landscape (n=3, m=1) ===")

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
    P_H = Hinv @ Hdot
    f2_amb = -np.trace(P_H @ P_H) + np.trace(Hinv @ Hddot)

    def f2_vis(theta, phi_angle):
        c = np.array([
            np.sin(theta) * np.cos(phi_angle),
            np.sin(theta) * np.sin(phi_angle),
            np.cos(theta),
        ])
        c_row = c.reshape(1, 3)
        _, _, Vt = svd(c_row)
        Z = Vt[1:].T
        R = Z.T @ H @ Z
        U_h = Z.T @ Hdot @ Z
        U_h2 = Z.T @ Hddot @ Z
        Rinv = inv(R)
        P_R = Rinv @ U_h
        f2_hid = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)
        return f2_amb - f2_hid

    # Sample the landscape
    n_theta = 50
    n_phi = 100
    thetas = np.linspace(0.01, np.pi - 0.01, n_theta)
    phis = np.linspace(0, 2*np.pi, n_phi)

    values = np.zeros((n_theta, n_phi))
    for i, th in enumerate(thetas):
        for j, ph in enumerate(phis):
            values[i, j] = f2_vis(th, ph)

    print(f"  f''_vis landscape on S^2:")
    print(f"    Min: {values.min():.6f}")
    print(f"    Max: {values.max():.6f}")
    print(f"    Mean: {values.mean():.6f}")
    print(f"    Std: {values.std():.6f}")
    print(f"    Ambient: {f2_amb:.6f}")

    # Find the location of max and min
    i_max, j_max = np.unravel_index(values.argmax(), values.shape)
    i_min, j_min = np.unravel_index(values.argmin(), values.shape)

    th_max, ph_max = thetas[i_max], phis[j_max]
    th_min, ph_min = thetas[i_min], phis[j_min]

    c_max = np.array([np.sin(th_max)*np.cos(ph_max), np.sin(th_max)*np.sin(ph_max), np.cos(th_max)])
    c_min = np.array([np.sin(th_min)*np.cos(ph_min), np.sin(th_min)*np.sin(ph_min), np.cos(th_min)])

    print(f"\n    Best direction:  theta={th_max:.2f}, phi={ph_max:.2f}")
    print(f"    Worst direction: theta={th_min:.2f}, phi={ph_min:.2f}")

    # Compute first-order quantities at the optimal observer
    c_opt = c_max.reshape(1, 3)
    Phi_opt = inv(c_opt @ Hinv @ c_opt.T).item()
    _, _, Vt = svd(c_opt)
    Z_opt = Vt[1:].T
    L_opt = Hinv @ c_opt.T * Phi_opt
    V_opt = (L_opt.T @ Hdot @ L_opt).item()

    print(f"\n    At optimal observer:")
    print(f"      Phi = {Phi_opt:.4f}")
    print(f"      V = {V_opt:.6f}")
    print(f"      f''_vis = {values.max():.6f}")
    print(f"      f''_hid = {f2_amb - values.max():.6f}")
    print(f"      Visible fraction of acceleration: {values.max()/f2_amb:.4f}")

    # Conservation check at every point
    conservation_errors = []
    for i, th in enumerate(thetas):
        for j, ph in enumerate(phis):
            c = np.array([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)])
            c_row = c.reshape(1, 3)
            _, _, Vt = svd(c_row)
            Z = Vt[1:].T
            R = Z.T @ H @ Z
            U_h = Z.T @ Hdot @ Z
            U_h2 = Z.T @ Hddot @ Z
            Rinv = inv(R)
            P_R = Rinv @ U_h
            f2_hid = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)
            conservation_errors.append(abs(values[i,j] + f2_hid - f2_amb))

    max_cons_err = max(conservation_errors)
    report(f"Conservation on full S^2 ({n_theta*n_phi} points)", max_cons_err, tol=1e-8)

    return values


def alignment_principle():
    """
    Test the alignment principle: the optimal observer aligns with the
    direction where H changes fastest (largest eigenvalue of Hdot weighted by H^{-1}).

    HYPOTHESIS: the optimal c* for f''_vis is related to the eigenvectors of
    the "effective acceleration" matrix H^{-1} Hddot - (H^{-1} Hdot)^2.

    If this is true, then the optimal observer can be found from a single
    eigendecomposition, not from optimisation on the Grassmannian.
    """
    print("\n=== R5: Alignment principle ===")

    n = 3
    rng = np.random.default_rng(42)

    for trial in range(5):
        H = np.eye(n) + rng.standard_normal((n,n)) @ rng.standard_normal((n,n)).T
        Hdot = sym(rng.standard_normal((n, n)))
        Hddot = sym(rng.standard_normal((n, n)))
        Hinv = inv(H)

        P_H = Hinv @ Hdot
        f2_amb = -np.trace(P_H @ P_H) + np.trace(Hinv @ Hddot)

        # The "effective second jet" of H in the H^{-1} metric
        # d^2/dt^2 [log det H] = Tr(H^{-1} Hddot - (H^{-1} Hdot)^2)
        # The per-direction version: for c, the visible part captures
        # the c-projection of this.
        # Effective acceleration matrix:
        A_eff = Hinv @ Hddot - P_H @ P_H  # n x n
        A_eff_sym = sym(A_eff)
        eigvals_A, eigvecs_A = eigh(A_eff_sym)

        # Find optimal observer by optimisation
        def neg_f2_vis(angles):
            c = np.zeros(n)
            c[0] = np.cos(angles[0])
            c[1] = np.sin(angles[0]) * np.cos(angles[1])
            c[2] = np.sin(angles[0]) * np.sin(angles[1])
            c_row = c.reshape(1, n)
            _, _, Vt = svd(c_row)
            Z = Vt[1:].T
            R = Z.T @ H @ Z; Rinv = inv(R)
            U_h = Z.T @ Hdot @ Z; U_h2 = Z.T @ Hddot @ Z
            P_R = Rinv @ U_h
            f2_hid = -np.trace(P_R @ P_R) + np.trace(Rinv @ U_h2)
            return -(f2_amb - f2_hid)

        best = float('inf')
        best_c = None
        for _ in range(100):
            a0 = rng.uniform(0, np.pi, 2)
            try:
                res = minimize(neg_f2_vis, a0, method='Nelder-Mead', options={'maxiter': 500})
                if res.fun < best:
                    best = res.fun
                    c = np.zeros(n)
                    c[0] = np.cos(res.x[0])
                    c[1] = np.sin(res.x[0]) * np.cos(res.x[1])
                    c[2] = np.sin(res.x[0]) * np.sin(res.x[1])
                    best_c = c / norm(c)
            except:
                pass

        # Check alignment with A_eff eigenvectors
        alignments = [abs(np.dot(best_c, eigvecs_A[:, i])) for i in range(n)]
        best_align = max(alignments)
        best_eigvec_idx = alignments.index(best_align)

        print(f"  Trial {trial}: A_eff eigvals = [{', '.join(f'{e:.2f}' for e in eigvals_A)}], "
              f"best alignment = {best_align:.4f} with eigvec {best_eigvec_idx} "
              f"(eigval {eigvals_A[best_eigvec_idx]:.2f})")

    print(f"\n  If alignment is consistently high with the largest eigenvector of A_eff,")
    print(f"  then the optimal observer aligns with the direction of maximum")
    print(f"  effective acceleration — confirming the alignment principle.")


if __name__ == "__main__":
    print("=" * 70)
    print("R5: Optimal Observer Alignment")
    print("=" * 70)

    optimal_observer_scalar(3, seed=0)
    observer_landscape(seed=0)
    alignment_principle()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
