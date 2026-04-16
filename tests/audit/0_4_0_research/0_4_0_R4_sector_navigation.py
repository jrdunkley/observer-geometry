"""
R4: Sector navigation via A_cpl eigenvalues.

Track source law diagnostics along a path and classify the local evidence sector.
The A_cpl eigenvalue structure signals which TN4 sector we're in:
  - All eigenvalues finite, V > 0: regular interior (source law valid)
  - A_cpl diverging: restiffening boundary (approaching Phi = 0)
  - V dropping rank: support boundary (information destruction)
  - A_cpl sign change: evidence curvature flip

Also track the conservation law decomposition along the path.
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

def sym(M):
    return 0.5 * (M + M.T)


def sector_navigator(H_func, Hdot_func, Hddot_func, C, ts, name="path"):
    """
    Given a path H(t) with fixed C, track the sector diagnostics at each t.

    Returns a list of diagnostic dicts.
    """
    n = C.shape[1]
    m = C.shape[0]
    _, _, Vt = svd(C)
    Z = Vt[m:].T

    diagnostics = []

    for t in ts:
        H = H_func(t)
        Hdot = Hdot_func(t)
        Hddot = Hddot_func(t)

        ev_H = eigh(H)[0]
        if np.min(ev_H) < 1e-10:
            diagnostics.append({'t': t, 'sector': 'H_degenerate', 'phi_min': 0, 'v_min': 0, 'a_cpl': None})
            continue

        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        R = Z.T @ H @ Z
        Rinv = inv(R)

        phi_eigvals = sorted(eigh(Phi)[0])
        phi_min = phi_eigvals[0]

        V = L.T @ Hdot @ L
        v_eigvals = sorted(eigh(V)[0])
        v_min = v_eigvals[0]

        B = L.T @ Hdot @ Z

        # Connection
        dHinv = -Hinv @ Hdot @ Hinv
        dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
        dL = dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi
        alpha = inv(Phi) @ (L.T @ H @ dL)

        Vdot = dL.T @ Hdot @ L + L.T @ Hddot @ L + L.T @ Hdot @ dL
        W = Vdot - alpha.T @ V - V @ alpha

        # Hidden defect
        Qhat = B @ Rinv @ B.T

        # Conservation law
        U_h = Z.T @ Hdot @ Z
        vis_rate = np.trace(inv(Phi) @ V) + 2 * np.trace(alpha)
        hid_rate = np.trace(inv(R) @ U_h)
        amb_rate = np.trace(Hinv @ Hdot)
        hidden_fraction = hid_rate / amb_rate if abs(amb_rate) > 1e-15 else float('nan')

        # A_cpl computation
        if v_min > 1e-8:
            Vsqrt = np.real(sqrtm(V))
            Vsqrt_inv = inv(Vsqrt)
            A_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)
            a_eigvals = sorted(eigh(A_cpl)[0])

            # Sector classification
            if all(e < -1e-6 for e in a_eigvals):
                sector = 'evidence_trough'  # all negative = convex log-evidence
            elif all(e > 1e-6 for e in a_eigvals):
                sector = 'evidence_peak'    # all positive = concave log-evidence
            elif min(a_eigvals) < -1e-6 and max(a_eigvals) > 1e-6:
                sector = 'evidence_saddle'  # indefinite
            else:
                sector = 'regular_flat'     # near-zero A_cpl
        elif v_min < -1e-8:
            sector = 'information_loss'
            a_eigvals = None
            A_cpl = None
        else:
            sector = 'support_boundary'
            a_eigvals = None
            A_cpl = None

        diagnostics.append({
            't': t,
            'sector': sector,
            'phi_min': phi_min,
            'v_min': v_min,
            'a_cpl_eigvals': a_eigvals,
            'qhat_trace': np.trace(Qhat),
            'vis_rate': vis_rate,
            'hid_rate': hid_rate,
            'hidden_fraction': hidden_fraction,
        })

    return diagnostics


def navigate_coupled_spring():
    """Navigate the coupled spring across coupling strength."""
    print("\n=== R4: Sector navigation — coupled spring ===")

    k1, k2 = 3.0, 2.0
    C = np.array([[1.0, 0.0]])

    def H(kc): return np.array([[k1 + kc, -kc], [-kc, k2 + kc]])
    def Hdot(kc): return np.array([[1.0, -1.0], [-1.0, 1.0]])
    def Hddot(kc): return np.zeros((2, 2))

    kcs = np.linspace(0.01, 30.0, 150)
    diags = sector_navigator(H, Hdot, Hddot, C, kcs, "coupled_spring")

    print(f"  {'kc':>6s}  {'sector':>20s}  {'Phi_min':>10s}  {'V_min':>10s}  "
          f"{'A_cpl':>10s}  {'Qhat_tr':>10s}  {'hid_frac':>10s}")
    for i in [0, 10, 30, 60, 100, 149]:
        d = diags[i]
        acpl_str = f"{d['a_cpl_eigvals'][0]:.4f}" if d['a_cpl_eigvals'] else "N/A"
        print(f"  {d['t']:6.2f}  {d['sector']:>20s}  {d['phi_min']:10.4f}  {d['v_min']:10.6f}  "
              f"{acpl_str:>10s}  {d['qhat_trace']:10.6f}  {d['hidden_fraction']:10.4f}")

    # Count sectors
    sector_counts = {}
    for d in diags:
        s = d['sector']
        sector_counts[s] = sector_counts.get(s, 0) + 1
    print(f"\n  Sector distribution: {sector_counts}")


def navigate_restiffening():
    """Navigate a path that approaches a branch-restiffening point."""
    print("\n=== R4: Sector navigation — restiffening approach ===")

    eps = 0.01
    C = np.array([[1.0, 0.0]])

    def H(t): return np.array([[eps + t**2, 0.0], [0.0, 1.0]])
    def Hdot(t): return np.array([[2*t, 0.0], [0.0, 0.0]])
    def Hddot(t): return np.array([[2.0, 0.0], [0.0, 0.0]])

    ts = np.linspace(0.02, 2.0, 100)
    diags = sector_navigator(H, Hdot, Hddot, C, ts, "restiffening")

    print(f"  {'t':>6s}  {'sector':>20s}  {'Phi_min':>10s}  {'V_min':>10s}  "
          f"{'A_cpl':>10s}  {'Qhat_tr':>10s}")
    for i in [0, 2, 5, 10, 20, 50, 99]:
        d = diags[i]
        acpl_str = f"{d['a_cpl_eigvals'][0]:.4f}" if d['a_cpl_eigvals'] else "N/A"
        print(f"  {d['t']:6.3f}  {d['sector']:>20s}  {d['phi_min']:10.6f}  {d['v_min']:10.6f}  "
              f"{acpl_str:>10s}  {d['qhat_trace']:10.6f}")

    # Track A_cpl divergence
    acpls = [d['a_cpl_eigvals'][0] for d in diags if d['a_cpl_eigvals'] is not None]
    ts_valid = [d['t'] for d in diags if d['a_cpl_eigvals'] is not None]

    print(f"\n  A_cpl divergence profile:")
    print(f"    t = {ts_valid[0]:.3f}: A_cpl = {acpls[0]:.4f}")
    print(f"    t = {ts_valid[-1]:.3f}: A_cpl = {acpls[-1]:.4f}")
    print(f"    Divergence rate: A_cpl ~ -1/(2t)")
    print(f"    At t = {ts_valid[0]:.3f}: exact = {-1/(2*ts_valid[0]):.4f}, computed = {acpls[0]:.4f}")

    report("A_cpl ~ -1/(2t) at small t", abs(acpls[0] - (-1/(2*ts_valid[0]))), tol=0.1)

    sector_counts = {}
    for d in diags:
        s = d['sector']
        sector_counts[s] = sector_counts.get(s, 0) + 1
    print(f"  Sector distribution: {sector_counts}")


def navigate_saddle():
    """Navigate a path that passes through an evidence saddle."""
    print("\n=== R4: Sector navigation — evidence saddle ===")

    C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    n, m = 3, 2

    # H(t) interpolates between two SPD matrices with different eigenstructure
    H0 = np.diag([3.0, 1.0, 2.0])
    H1 = np.array([[1.0, 0.5, 0.0], [0.5, 3.0, 0.0], [0.0, 0.0, 2.0]])
    # Hdot = H1 - H0 = [[-2, 0.5, 0], [0.5, 2, 0], [0, 0, 0]]
    dH_const = H1 - H0

    def H(t): return (1-t) * H0 + t * H1
    def Hdot(t): return dH_const
    def Hddot(t): return np.zeros((n, n))

    ts = np.linspace(0.01, 0.99, 100)
    diags = sector_navigator(H, Hdot, Hddot, C, ts, "saddle_path")

    print(f"  Path: H(t) = (1-t)*diag(3,1,2) + t*[[1,.5,0],[.5,3,0],[0,0,2]]")
    print(f"  {'t':>6s}  {'sector':>20s}  {'V_min':>10s}  {'A_cpl':>24s}")
    for i in range(0, 100, 10):
        d = diags[i]
        if d['a_cpl_eigvals']:
            acpl_str = f"[{d['a_cpl_eigvals'][0]:.4f}, {d['a_cpl_eigvals'][1]:.4f}]"
        else:
            acpl_str = "N/A"
        print(f"  {d['t']:6.3f}  {d['sector']:>20s}  {d['v_min']:10.4f}  {acpl_str:>24s}")

    sector_counts = {}
    for d in diags:
        s = d['sector']
        sector_counts[s] = sector_counts.get(s, 0) + 1
    print(f"\n  Sector distribution: {sector_counts}")

    # Check for sector transitions
    transitions = []
    for i in range(1, len(diags)):
        if diags[i]['sector'] != diags[i-1]['sector']:
            transitions.append((diags[i-1]['t'], diags[i-1]['sector'], diags[i]['sector']))

    if transitions:
        print(f"\n  Sector transitions:")
        for t, s1, s2 in transitions:
            print(f"    t = {t:.3f}: {s1} -> {s2}")
    else:
        print(f"\n  No sector transitions (single sector throughout)")


if __name__ == "__main__":
    print("=" * 70)
    print("R4: Sector Navigation via A_cpl Eigenvalues")
    print("=" * 70)

    navigate_coupled_spring()
    navigate_restiffening()
    navigate_saddle()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
