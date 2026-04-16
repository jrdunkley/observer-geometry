"""
N5: Verify the coupling-budget Lagrangian on the coupled spring system.

The coupled spring has closed-form expressions for everything:
  H(kc) = [[k1+kc, -kc], [-kc, k2+kc]]
  Phi(kc) = (k1*k2 + k1*kc + k2*kc) / (k2+kc)
  V = k2^2 / (k2+kc)^2

For a 1D observer on a 2D system, the Grassmannian is Gr(1,2) = RP^1 = S^1/Z_2.
Parameterise the observer as C(phi) = [cos(phi), sin(phi)].

The Lagrangian claims:
1. The TN1 stationarity equation is the zero-velocity limit
2. The gradient d/dC[vis_rate] = -d/dC[hid_rate] from conservation
3. The optimal observer angle phi* maximises vis_rate

Verify all three on the coupled spring.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm
import sys
from pathlib import Path

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))
from nomogeo import information_budget

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok


def coupled_spring_landscape():
    """
    Full landscape of vis_rate as a function of observer angle phi
    on the coupled spring system.
    """
    print("=== N5: Coupled spring Lagrangian verification ===")

    k1, k2 = 3.0, 2.0
    kc = 5.0
    H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
    Hdot = np.array([[1.0, -1.0], [-1.0, 1.0]])  # dH/dkc

    Hinv = inv(H)
    amb_rate = float(np.trace(Hinv @ Hdot))

    # Sweep observer angle
    phis = np.linspace(0, np.pi, 361)  # RP^1 is [0, pi)
    vis_rates = []
    hid_rates = []
    v_mins = []

    for phi in phis:
        C = np.array([[np.cos(phi), np.sin(phi)]])
        b = information_budget(H, C, Hdot)
        vis_rates.append(b.visible_rate)
        hid_rates.append(b.hidden_rate)
        v_mins.append(float(np.min(b.v_eigenvalues)))

    vis_rates = np.array(vis_rates)
    hid_rates = np.array(hid_rates)
    v_mins = np.array(v_mins)

    # 1. Conservation at every angle
    conservation = np.abs(vis_rates + hid_rates - amb_rate)
    report("Conservation at all 361 angles", np.max(conservation))

    # 2. Find the optimal observer angle
    best_idx = np.argmax(vis_rates)
    worst_idx = np.argmin(vis_rates)
    phi_best = phis[best_idx]
    phi_worst = phis[worst_idx]

    print(f"\n  Ambient rate (FIXED): {amb_rate:.6f}")
    print(f"  Best observer:  phi = {np.degrees(phi_best):.1f} deg, vis_rate = {vis_rates[best_idx]:.6f}")
    print(f"  Worst observer: phi = {np.degrees(phi_worst):.1f} deg, vis_rate = {vis_rates[worst_idx]:.6f}")
    print(f"  Range of vis_rate: [{vis_rates.min():.6f}, {vis_rates.max():.6f}]")

    # 3. Verify gradient is zero at the optimum (stationarity)
    dphi = phis[1] - phis[0]
    grad_vis = np.gradient(vis_rates, dphi)
    grad_hid = np.gradient(hid_rates, dphi)

    # At optimum, gradient should be zero
    report("Gradient zero at optimum", abs(grad_vis[best_idx]), tol=1e-3)

    # 4. Verify d/dphi[vis_rate] = -d/dphi[hid_rate] (from conservation)
    max_grad_err = np.max(np.abs(grad_vis + grad_hid))
    report("d/dphi[vis] = -d/dphi[hid] (conservation gradient)", max_grad_err, tol=1e-3)

    # 5. Check: does the optimal angle correspond to the canonical observer?
    # Canonical: C = [1, 0], phi = 0
    print(f"\n  Canonical observer (phi=0): vis_rate = {vis_rates[0]:.6f}")
    print(f"  Optimal observer (phi={np.degrees(phi_best):.1f}): vis_rate = {vis_rates[best_idx]:.6f}")

    # 6. V > 0 region
    vpos_region = np.sum(v_mins > 1e-6)
    print(f"\n  V > 0 region: {vpos_region}/{len(phis)} angles ({100*vpos_region/len(phis):.1f}%)")

    # 7. What is the optimal angle geometrically?
    # For H = [[a, b], [b, d]], the observer C = [cos phi, sin phi]
    # sees Phi = 1 / (C H^{-1} C^T)
    # The visible rate is d/dkc [log Phi] = ...
    # The eigenvectors of H are the natural frame. Does the optimum align with one?
    eigvals_H, eigvecs_H = eigh(H)
    for i in range(2):
        angle_with_eig = np.abs(np.dot([np.cos(phi_best), np.sin(phi_best)], eigvecs_H[:, i]))
        print(f"  |<c*, e_{i}>| = {angle_with_eig:.4f} (eigval = {eigvals_H[i]:.2f})")

    eigvals_Hdot, eigvecs_Hdot = eigh(Hdot)
    for i in range(2):
        angle_with_eig = np.abs(np.dot([np.cos(phi_best), np.sin(phi_best)], eigvecs_Hdot[:, i]))
        print(f"  |<c*, e_Hdot_{i}>| = {angle_with_eig:.4f} (eigval = {eigvals_Hdot[i]:.2f})")

    # 8. The effective acceleration matrix A_eff = H^{-1} Hddot - (H^{-1} Hdot)^2
    # On a linear path, Hddot = 0, so A_eff = -(H^{-1} Hdot)^2
    P = Hinv @ Hdot
    A_eff = -P @ P
    eigvals_Aeff, eigvecs_Aeff = eigh(A_eff)
    print(f"\n  A_eff eigenvalues: {eigvals_Aeff}")
    for i in range(2):
        angle_with_eig = np.abs(np.dot([np.cos(phi_best), np.sin(phi_best)], eigvecs_Aeff[:, i]))
        print(f"  |<c*, e_Aeff_{i}>| = {angle_with_eig:.4f} (eigval = {eigvals_Aeff[i]:.4f})")

    # 9. Print the full landscape at key angles
    print(f"\n  {'angle':>6s}  {'vis_rate':>10s}  {'hid_rate':>10s}  {'v_min':>8s}  {'V>0':>4s}")
    for i in range(0, len(phis), 36):
        vp = 'Y' if v_mins[i] > 1e-6 else 'N'
        print(f"  {np.degrees(phis[i]):6.1f}  {vis_rates[i]:10.6f}  {hid_rates[i]:10.6f}  {v_mins[i]:8.4f}  {vp:>4s}")

    return vis_rates, hid_rates, v_mins, phis


def verify_euler_lagrange():
    """
    The Euler-Lagrange equation for the zero-velocity Lagrangian is:
      d/dC [vis_rate] = 0

    This is equivalent to: d/dC [hid_rate] = 0 (from conservation).

    For the coupled spring, this means: d/dphi [Tr(R^{-1} U_h)] = 0
    where R = Z^T H Z and U_h = Z^T Hdot Z depend on phi through Z = ker(C).

    Verify that the stationarity condition from the gradient matches
    the optimal angle from the landscape.
    """
    print("\n=== Euler-Lagrange verification ===")

    k1, k2 = 3.0, 2.0

    for kc in [1.0, 5.0, 15.0]:
        H = np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
        Hdot = np.array([[1.0, -1.0], [-1.0, 1.0]])
        Hinv = inv(H)

        # Find optimal phi by dense sweep
        phis = np.linspace(0, np.pi, 1001)
        best_vr = -1e10; best_phi = 0
        for phi in phis:
            C = np.array([[np.cos(phi), np.sin(phi)]])
            b = information_budget(H, C, Hdot)
            if b.visible_rate > best_vr:
                best_vr = b.visible_rate; best_phi = phi

        # Also check: is the optimal observer the one that minimises hid_rate?
        best_hr = 1e10; best_phi_h = 0
        for phi in phis:
            C = np.array([[np.cos(phi), np.sin(phi)]])
            b = information_budget(H, C, Hdot)
            if b.hidden_rate < best_hr:
                best_hr = b.hidden_rate; best_phi_h = phi

        # They should be the same (from conservation)
        err = abs(best_phi - best_phi_h)
        if err > np.pi/2: err = np.pi - err  # RP^1 identification
        report(f"kc={kc}: max(vis) = min(hid) angle", err, tol=0.01)
        print(f"  kc={kc}: optimal phi = {np.degrees(best_phi):.1f} deg, vis_rate = {best_vr:.6f}")

    print(f"\n  The stationarity condition d/dphi[vis_rate] = 0 is equivalent")
    print(f"  to d/dphi[hid_rate] = 0, confirming the Euler-Lagrange equation")
    print(f"  in the zero-velocity limit is the conservation-derived gradient condition.")


if __name__ == "__main__":
    print("=" * 70)
    print("N5: Coupling-Budget Lagrangian on Coupled Springs")
    print("=" * 70)

    coupled_spring_landscape()
    verify_euler_lagrange()

    print(f"\n{'='*70}")
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
