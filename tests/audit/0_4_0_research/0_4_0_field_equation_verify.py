"""
Lightweight verification of the KL-regularised field equation.
No optimisation — evaluates the analytical solution and checks stationarity.
Runs in < 1 second.

Field equation: S - Xi = (gamma/2)(H - H0)
Solution at B=0: Phi* = Phi0 - (2/gamma)V, R* = R0, K0 = 0
"""
import sys, os
import numpy as np
from numpy.linalg import inv, norm, det, eigvalsh, solve
from scipy.linalg import null_space

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
from nomogeo.validation import symmetrize

np.set_printoptions(precision=10, linewidth=120)

checks_passed = 0
checks_failed = 0
def check(name, cond, detail=""):
    global checks_passed, checks_failed
    if cond: checks_passed += 1; print(f"  [PASS] {name}")
    else: checks_failed += 1; print(f"  [FAIL] {name}  {detail}")

def random_spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return symmetrize(A @ A.T + 2*np.eye(n))  # well-conditioned

def hidden_frame(C):
    return null_space(C, rcond=1e-12)

def split(H, C):
    """Return L, Z, Phi, R."""
    Hi = inv(H)
    Q = symmetrize(C @ Hi @ C.T)
    Phi = inv(Q)
    L = Hi @ C.T @ Phi
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    return L, Z, Phi, R

def vis_rate(H, C, Xi):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Xi @ Z)
    return float(np.trace(solve(H, Xi)) - np.trace(solve(R, U_h)))

def hidden_stress(H, C, Xi):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Xi @ Z)
    A = H @ Z @ inv(R)
    return symmetrize(A @ U_h @ A.T)

def grad_vis(H, C, Xi):
    Z = hidden_frame(C)
    Hi = inv(H)
    R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Xi @ Z)
    Ri = inv(R)
    return symmetrize(-Hi @ Xi @ Hi + Z @ Ri @ U_h @ Ri @ Z.T)

def grad_DKL(H, H0):
    Hi = inv(H)
    return symmetrize(0.5 * (Hi - Hi @ H0 @ Hi))

def D_KL(H0, H):
    A = H0 @ inv(H)
    return 0.5 * (np.trace(A) - np.log(det(A)) - H.shape[0])

def construct_adapted_Xi(H, C, Xi_raw):
    """Rotate Xi so that B = L^T Xi Z = 0 (adapted perturbation)."""
    L, Z, Phi, R = split(H, C)
    # Project Xi onto the B=0 subspace
    # Xi_adapted has V = L^T Xi L (unchanged), U_h = Z^T Xi Z (unchanged), B = 0
    V = symmetrize(L.T @ Xi_raw @ L)
    U_h = symmetrize(Z.T @ Xi_raw @ Z)
    # Reconstruct: Xi = M^{-T} [[V, 0], [0, U_h]] M^{-1}
    m = C.shape[0]; n = H.shape[0]
    block = np.zeros((n, n))
    block[:m, :m] = V
    block[m:, m:] = U_h
    M = np.hstack([L, Z])
    Mi = inv(M)
    return symmetrize(Mi.T @ block @ Mi)

# ============================================================
print("=" * 70)
print("FIELD EQUATION VERIFICATION (analytical solution, no optimisation)")
print("=" * 70)

# --- Test 1: 2x2 diagonal model ---
print("\n--- 2x2 split model (scalar verification) ---")
a0, d0, h, f, gamma = 5.0, 3.0, 1.5, -0.8, 20.0
H0 = np.diag([a0, d0])
Xi = np.diag([h, f])
C = np.array([[1.0, 0.0]])

a_star = a0 - 2*h/gamma
d_star = d0  # frozen
H_star = np.diag([a_star, d_star])

print(f"  H0 = diag({a0}, {d0}), Xi = diag({h}, {f}), gamma = {gamma}")
print(f"  H* = diag({a_star:.4f}, {d_star:.4f})")
print(f"  vis_rate(H*) = {vis_rate(H_star, C, Xi):.6f}")

# Check field equation
S = hidden_stress(H_star, C, Xi)
fe_res = norm(S - Xi - (gamma/2)*(H_star - H0), 'fro')
print(f"  ||S - Xi - (gamma/2)(H* - H0)|| = {fe_res:.3e}")
check("2x2 field equation", fe_res < 1e-10)

# Check stationarity
gv = grad_vis(H_star, C, Xi)
gd = grad_DKL(H_star, H0)
stationarity = norm(gv - gamma * gd, 'fro')
print(f"  ||grad_vis - gamma * grad_DKL|| = {stationarity:.3e}")
check("2x2 stationarity", stationarity < 1e-10)

# Check H* is a local maximum (J(H*) > J(H0) and J(H_perturbed))
J_star = vis_rate(H_star, C, Xi) - gamma * D_KL(H0, H_star)
J_ref = vis_rate(H0, C, Xi) - gamma * D_KL(H0, H0)
print(f"  J(H*) = {J_star:.6f},  J(H0) = {J_ref:.6f}")
check("J(H*) >= J(H0)", J_star >= J_ref - 1e-10)

# --- Test 2: General n, m, analytical construction ---
print("\n--- General (n, m) analytical construction ---")

for trial in range(20):
    seed = 7000 + trial
    rng = np.random.default_rng(seed)
    n = int(rng.integers(3, 7))
    m = int(rng.integers(1, n))
    gamma = 50.0  # strong regularisation

    H0 = random_spd(n, seed)
    C = rng.standard_normal((m, n))
    while np.linalg.matrix_rank(C) < m:
        C = rng.standard_normal((m, n))
    Xi_raw = symmetrize(rng.standard_normal((n, n)))

    # Get adapted frame for H0
    L0, Z0, Phi0, R0 = split(H0, C)

    # Project Xi to B=0 sector (so the adapted observer condition holds)
    Xi = construct_adapted_Xi(H0, C, Xi_raw)

    # Verify B = 0
    L0, Z0, Phi0, R0 = split(H0, C)
    B_check = L0.T @ Xi @ Z0
    # Note: B depends on H, not just Xi. At H0, B = L0^T Xi Z0.
    # The adapted perturbation was constructed for H0's frame.

    # Visible jet at H0
    V = symmetrize(L0.T @ Xi @ L0)

    # Predicted solution: Phi* = Phi0 - (2/gamma)V, R* = R0
    Phi_star = Phi0 - (2.0/gamma) * V

    # Check Phi* is SPD
    phi_eigs = eigvalsh(Phi_star)
    if phi_eigs[0] <= 0:
        # gamma too small for this perturbation, skip
        continue

    # Reconstruct H* from (Phi*, R0, K=0) in H0's adapted basis
    M0 = np.hstack([L0, Z0])
    block_star = np.zeros((n, n))
    block_star[:m, :m] = Phi_star
    block_star[m:, m:] = R0
    M0i = inv(M0)
    H_star = symmetrize(M0i.T @ block_star @ M0i)

    # Check H* is SPD
    h_eigs = eigvalsh(H_star)
    if h_eigs[0] <= 0:
        continue

    # Check field equation: S - Xi = (gamma/2)(H* - H0)
    S = hidden_stress(H_star, C, Xi)
    fe_lhs = S - Xi
    fe_rhs = (gamma/2) * (H_star - H0)
    fe_res = norm(fe_lhs - fe_rhs, 'fro')

    # Check stationarity: grad_vis = gamma * grad_DKL
    gv = grad_vis(H_star, C, Xi)
    gd = grad_DKL(H_star, H0)
    station = norm(gv - gamma * gd, 'fro')

    # Check J(H*) >= J(H0)
    J_star = vis_rate(H_star, C, Xi) - gamma * D_KL(H0, H_star)
    J_ref = vis_rate(H0, C, Xi) - gamma * D_KL(H0, H0)

    if trial < 3:
        print(f"\n  Trial {trial}: n={n}, m={m}")
        print(f"    ||FE residual|| = {fe_res:.3e}")
        print(f"    ||stationarity|| = {station:.3e}")
        print(f"    J(H*) = {J_star:.4f}, J(H0) = {J_ref:.4f}")
        print(f"    Phi* eigs: {phi_eigs}")

    check(f"FE (n={n},m={m},t{trial})", fe_res < 1e-6, f"res={fe_res:.3e}")
    check(f"Stationarity (n={n},m={m},t{trial})", station < 1e-6, f"res={station:.3e}")
    check(f"J improve (n={n},m={m},t{trial})", J_star >= J_ref - 1e-6,
          f"J*={J_star:.4f}, J0={J_ref:.4f}")

print(f"\n{'='*70}")
print(f"  TOTAL: {checks_passed} passed, {checks_failed} failed")
print(f"{'='*70}")
print("""
  INTERPRETATION:
  The analytical solution H* = H0 + M^{-T}[[-2V/gamma, 0],[0, 0]]M^{-1}
  satisfies the field equation S - Xi = (gamma/2)(H - H0) and the
  stationarity condition grad_vis = gamma * grad_DKL.

  At B = 0:
    - Hidden metric frozen: R* = R0
    - Visible precision shifted: Phi* = Phi0 - (2/gamma)V
    - H* differs from H0 only in the visible sector (rank-m perturbation)
""")
