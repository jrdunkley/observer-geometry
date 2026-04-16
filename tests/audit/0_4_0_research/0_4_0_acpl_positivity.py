"""
Target 1: Is A_cpl >= 0 a theorem on linear paths?

For a linear path H(t) = (1-t)H_0 + t H_1 with fixed C:
  Hdot = H_1 - H_0 (constant)
  Hddot = 0

From the completed square (eq 2.8 of 0.4.0 TN):
  W = L^T Hddot L - 2 Qhat + 1/2 O
  With Hddot = 0 and beta = 0 (fixed C): W = -2 Qhat, O = 0

  A_cpl = -1/2 V^{-1/2} W V^{-1/2} = V^{-1/2} Qhat V^{-1/2}

  where Qhat = B R^{-1} B^T (PSD by construction).

THEOREM: On any linear path through SPD(n) with fixed observer C,
if V > 0 on the active support, then A_cpl is PSD.

PROOF: A_cpl = V^{-1/2} (B R^{-1} B^T) V^{-1/2}.
Since R is SPD, R^{-1} is SPD. Therefore B R^{-1} B^T is PSD (it's X D X^T
with D SPD). Conjugation by V^{-1/2} preserves positive semidefiniteness.

COROLLARY: On linear paths, the log-evidence curvature source term
-2 Tr(Phi^{-1} V^{1/2} A_cpl V^{1/2}) = -2 Tr(Phi^{-1} Qhat) <= 0.
The evidence is always concave (decelerating) in the source-term contribution.

This is NOT true for nonlinear paths (Hddot != 0), where A_direct =
-1/2 V^{-1/2} (L^T Hddot L) V^{-1/2} can be negative and can dominate
the positive hidden defect.

Let's verify this theorem comprehensively and explore its boundary.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd
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


def verify_positivity_theorem(num_trials=5000):
    """
    Verify: A_cpl >= 0 on linear paths (Hddot = 0, fixed C).
    """
    print("=== Theorem: A_cpl >= 0 on linear paths ===")

    violations = 0
    total_tested = 0
    min_acpl_eigval = float('inf')

    for trial in range(num_trials):
        rng = np.random.default_rng(trial)
        n = rng.integers(3, 8)
        m = rng.integers(1, n)

        H0 = spd(n, trial * 10)
        H1 = spd(n, trial * 10 + 1)
        C = rng.standard_normal((m, n))

        Hdot = H1 - H0  # linear path: Hddot = 0

        Hinv = inv(H0)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C)
        Z = Vt[m:].T
        R = Z.T @ H0 @ Z
        Rinv = inv(R)

        V = L.T @ Hdot @ L
        B = L.T @ Hdot @ Z

        v_eigvals = eigh(V)[0]
        if np.min(v_eigvals) < 0.01:
            continue  # V not positive, skip

        total_tested += 1

        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)

        # Qhat = B R^{-1} B^T (PSD by construction)
        Qhat = B @ Rinv @ B.T

        # A_cpl = V^{-1/2} Qhat V^{-1/2} (should be PSD)
        A_cpl = sym(Vsqrt_inv @ Qhat @ Vsqrt_inv)
        acpl_eigvals = eigh(A_cpl)[0]
        min_eigval = float(np.min(acpl_eigvals))

        if min_eigval < -1e-10:
            violations += 1

        min_acpl_eigval = min(min_acpl_eigval, min_eigval)

    report(f"A_cpl >= 0 on linear paths ({total_tested} trials, 0 violations expected)",
           float(violations))
    print(f"  Tested: {total_tested}/{num_trials} (rest had V not positive)")
    print(f"  Min A_cpl eigenvalue across all trials: {min_acpl_eigval:.6e}")
    print(f"  Violations: {violations}")

    return violations == 0


def verify_algebraic_proof():
    """
    Verify the algebraic chain: A_cpl = V^{-1/2} Qhat V^{-1/2}
    where Qhat = B R^{-1} B^T is PSD.
    """
    print("\n=== Algebraic proof verification ===")

    rng = np.random.default_rng(42)
    for n, m in [(3,1), (4,2), (5,3), (6,2), (8,4)]:
        H0 = spd(n, 42 + n)
        H1 = spd(n, 43 + n)
        C = rng.standard_normal((m, n))
        Hdot = H1 - H0

        Hinv = inv(H0)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C)
        Z = Vt[m:].T
        R = Z.T @ H0 @ Z
        Rinv = inv(R)

        V = L.T @ Hdot @ L
        B = L.T @ Hdot @ Z

        if np.min(eigh(V)[0]) < 0.01:
            continue

        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)

        # Hddot = 0 path: W = theta^T B^T + B theta = -2 B R^{-1} B^T
        theta = -Rinv @ B.T
        W_from_theta = theta.T @ B.T + B @ theta
        W_from_Qhat = -2.0 * B @ Rinv @ B.T

        report(f"(n={n},m={m}) W = -2 Qhat", norm(W_from_theta - W_from_Qhat))

        # A_cpl from W vs from Qhat
        A_from_W = sym(-0.5 * Vsqrt_inv @ W_from_theta @ Vsqrt_inv)
        A_from_Qhat = sym(Vsqrt_inv @ (B @ Rinv @ B.T) @ Vsqrt_inv)

        report(f"(n={n},m={m}) A_cpl = V^{{-1/2}} Qhat V^{{-1/2}}", norm(A_from_W - A_from_Qhat))

        # Verify Qhat is PSD
        qhat_eigvals = eigh(B @ Rinv @ B.T)[0]
        report(f"(n={n},m={m}) Qhat PSD", max(0, -np.min(qhat_eigvals)))

        # Verify A_cpl is PSD
        acpl_eigvals = eigh(A_from_Qhat)[0]
        report(f"(n={n},m={m}) A_cpl PSD", max(0, -np.min(acpl_eigvals)))


def boundary_analysis():
    """
    Explore the boundary: how much Hddot is needed to make A_cpl indefinite?

    A_cpl = A_direct + hidden_defect
    A_direct = -1/2 V^{-1/2} (L^T Hddot L) V^{-1/2}  (can be negative)
    hidden_defect = V^{-1/2} Qhat V^{-1/2}  (always PSD)

    A_cpl indefinite iff A_direct has a negative eigenvalue that exceeds
    the hidden defect in magnitude.
    """
    print("\n=== Boundary: Hddot needed to make A_cpl indefinite ===")

    rng = np.random.default_rng(100)
    n, m = 4, 2

    H0 = spd(n, 100)
    H1 = spd(n, 101)
    C = rng.standard_normal((m, n))
    Hdot = H1 - H0

    Hinv = inv(H0)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C)
    Z = Vt[m:].T
    R = Z.T @ H0 @ Z
    Rinv = inv(R)

    V = L.T @ Hdot @ L
    B = L.T @ Hdot @ Z

    if np.min(eigh(V)[0]) < 0.01:
        print("  Skipped: V not positive")
        return

    Vsqrt = np.real(sqrtm(V))
    Vsqrt_inv = inv(Vsqrt)
    Qhat = B @ Rinv @ B.T
    hidden_defect = sym(Vsqrt_inv @ Qhat @ Vsqrt_inv)
    hd_min_eigval = np.min(eigh(hidden_defect)[0])

    print(f"  Hidden defect min eigenvalue: {hd_min_eigval:.6f}")
    print(f"  Hidden defect trace: {np.trace(hidden_defect):.6f}")

    # Sweep Hddot magnitude
    # Choose Hddot that pushes A_direct negative
    Hddot_dir = sym(rng.standard_normal((n, n)))
    # A_direct = -1/2 V^{-1/2} (L^T Hddot L) V^{-1/2}
    A_dir_unit = sym(-0.5 * Vsqrt_inv @ (L.T @ Hddot_dir @ L) @ Vsqrt_inv)
    # Find the most negative eigenvalue of A_dir_unit
    ad_eigvals = eigh(A_dir_unit)[0]
    if np.min(ad_eigvals) >= 0:
        # A_dir is PSD for this Hddot direction; try negative
        Hddot_dir = -Hddot_dir
        A_dir_unit = sym(-0.5 * Vsqrt_inv @ (L.T @ Hddot_dir @ L) @ Vsqrt_inv)
        ad_eigvals = eigh(A_dir_unit)[0]

    if np.min(ad_eigvals) >= 0:
        print("  Could not find Hddot that makes A_direct negative")
        return

    ad_min = np.min(ad_eigvals)
    print(f"  A_direct min eigenvalue (unit Hddot): {ad_min:.6f}")

    # Scale: A_cpl = scale * A_dir_unit + hidden_defect
    # Indefinite when scale * ad_min + hd_min_eigval < 0
    # => scale > -hd_min_eigval / ad_min (since ad_min < 0)
    critical_scale = -hd_min_eigval / ad_min if ad_min < 0 else float('inf')
    print(f"  Critical Hddot scale: {critical_scale:.6f}")
    print(f"  Below this: A_cpl PSD (evidence peak)")
    print(f"  Above this: A_cpl indefinite (evidence saddle)")

    # Verify
    scales = [0, 0.5 * critical_scale, 0.9 * critical_scale, critical_scale,
              1.1 * critical_scale, 2.0 * critical_scale, 5.0 * critical_scale]
    print(f"\n  {'scale':>10s}  {'A_cpl_min':>12s}  {'regime':>20s}")
    for scale in scales:
        Hddot = scale * Hddot_dir
        A_direct = sym(-0.5 * Vsqrt_inv @ (L.T @ Hddot @ L) @ Vsqrt_inv)
        A_cpl = A_direct + hidden_defect
        acpl_min = float(np.min(eigh(A_cpl)[0]))
        if acpl_min > 1e-8:
            regime = "evidence_peak (PSD)"
        elif acpl_min < -1e-8:
            regime = "evidence_saddle (indef)"
        else:
            regime = "critical boundary"
        print(f"  {scale:10.4f}  {acpl_min:12.6f}  {regime:>20s}")


def acpl_equals_zero_characterisation():
    """
    A_cpl = 0 iff Qhat = B R^{-1} B^T = 0 iff B = 0.

    B = L^T Hdot Z = 0 means the perturbation Hdot does not couple
    visible and hidden sectors through the observer.

    This happens when:
    1. Hdot is block-diagonal in the (visible, hidden) decomposition, or
    2. The perturbation is entirely in the visible or entirely in the hidden sector.

    Verify and characterise.
    """
    print("\n=== A_cpl = 0 characterisation ===")

    n, m = 4, 2
    rng = np.random.default_rng(200)

    H = spd(n, 200)
    C = rng.standard_normal((m, n))

    Hinv = inv(H)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C)
    Z = Vt[m:].T

    # Construct Hdot that is block-diagonal in (L, Z) frame
    # M = [L Z], M^T H M = diag(Phi, R)
    # A perturbation Hdot that is block-diagonal means:
    # L^T Hdot Z = 0, i.e., B = 0
    # In the M basis: Hdot_M = diag(dPhi_block, dR_block)
    # Hdot = M^{-T} diag(dPhi_block, dR_block) M^{-1}

    M = np.hstack([L, Z])
    M_inv = inv(M)
    dPhi_block = sym(rng.standard_normal((m, m)))
    dR_block = sym(rng.standard_normal((n-m, n-m)))
    Hdot_block = M_inv.T @ np.block([
        [dPhi_block, np.zeros((m, n-m))],
        [np.zeros((n-m, m)), dR_block]
    ]) @ M_inv

    V = L.T @ Hdot_block @ L
    B = L.T @ Hdot_block @ Z

    report("B = 0 for block-diagonal Hdot", norm(B))

    if np.min(eigh(V)[0]) > 0.01:
        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        Qhat = B @ inv(Z.T @ H @ Z) @ B.T
        A_cpl = sym(Vsqrt_inv @ Qhat @ Vsqrt_inv)
        report("A_cpl = 0 for block-diagonal Hdot", norm(A_cpl))
        print(f"  When the perturbation doesn't couple visible and hidden,")
        print(f"  the hidden defect vanishes and A_cpl = 0.")
        print(f"  The visible precision changes, but without hidden mediation.")


def trace_identity():
    """
    Verify: Tr(A_cpl) = Tr(V^{-1} Qhat) = Tr(V^{-1} B R^{-1} B^T).

    This scalar quantity is the total "hidden correction rate" divided by V.
    On linear paths it equals Tr(A_cpl) >= 0.
    """
    print("\n=== Trace identity ===")

    for trial in range(20):
        rng = np.random.default_rng(300 + trial)
        n = rng.integers(3, 7)
        m = rng.integers(1, n)

        H0 = spd(n, 300 + trial)
        H1 = spd(n, 301 + trial)
        C = rng.standard_normal((m, n))
        Hdot = H1 - H0

        Hinv = inv(H0)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C)
        Z = Vt[m:].T
        R = Z.T @ H0 @ Z
        Rinv = inv(R)

        V = L.T @ Hdot @ L
        B = L.T @ Hdot @ Z

        if np.min(eigh(V)[0]) < 0.01:
            continue

        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        Qhat = B @ Rinv @ B.T
        A_cpl = sym(Vsqrt_inv @ Qhat @ Vsqrt_inv)

        # Tr(A_cpl) = Tr(V^{-1} Qhat)
        V_inv = inv(V)
        tr_acpl = np.trace(A_cpl)
        tr_vinv_qhat = np.trace(V_inv @ Qhat)
        report(f"Tr(A_cpl) = Tr(V^{{-1}} Qhat)", abs(tr_acpl - tr_vinv_qhat))
        break  # one is enough for the identity


if __name__ == "__main__":
    print("=" * 70)
    print("Target 1: A_cpl Positivity Theorem")
    print("=" * 70)

    verify_positivity_theorem(5000)
    verify_algebraic_proof()
    boundary_analysis()
    acpl_equals_zero_characterisation()
    trace_identity()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
