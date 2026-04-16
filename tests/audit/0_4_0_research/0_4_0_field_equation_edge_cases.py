"""
Edge-case stress tests for the KL field equation.
All instant — no optimisation, just evaluates analytical formulas.

Run: cd C:\\observer_geometry_workspace && set PYTHONPATH=nomogeo\\src && python nomogeo\\audit\\0_4_0_research\\0_4_0_field_equation_edge_cases.py
"""
import sys, os
import numpy as np
from numpy.linalg import inv, norm, det, eigvalsh
from scipy.linalg import null_space

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
from nomogeo.validation import symmetrize

np.set_printoptions(precision=8, linewidth=120)

P, F = 0, 0
def check(name, cond, detail=""):
    global P, F
    if cond: P += 1; print(f"  [PASS] {name}")
    else: F += 1; print(f"  [FAIL] {name}  {detail}")

def random_spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return symmetrize(A @ A.T + 2*np.eye(n))

def split(H, C):
    Hi = inv(H); Q = symmetrize(C @ Hi @ C.T); Phi = inv(Q)
    L = Hi @ C.T @ Phi; Z = null_space(C, rcond=1e-12)
    R = symmetrize(Z.T @ H @ Z)
    return L, Z, Phi, R

def vis_rate(H, C, Xi):
    Z = null_space(C, rcond=1e-12)
    R = symmetrize(Z.T @ H @ Z); U = symmetrize(Z.T @ Xi @ Z)
    return float(np.trace(np.linalg.solve(H, Xi)) - np.trace(np.linalg.solve(R, U)))

def hidden_stress(H, C, Xi):
    Z = null_space(C, rcond=1e-12)
    R = symmetrize(Z.T @ H @ Z); U = symmetrize(Z.T @ Xi @ Z)
    A = H @ Z @ inv(R); return symmetrize(A @ U @ A.T)

def D_KL(H0, H):
    A = H0 @ inv(H); return 0.5*(np.trace(A) - np.log(det(A)) - H.shape[0])

def build_solution(H0, C, Xi, gamma):
    """Build H* from the closed-form formula."""
    L, Z, Phi0, R0 = split(H0, C)
    V = symmetrize(L.T @ Xi @ L)
    Phi_star = Phi0 - (2.0/gamma)*V
    n = H0.shape[0]; m = C.shape[0]
    M = np.hstack([L, Z]); Mi = inv(M)
    block = np.zeros((n,n)); block[:m,:m] = Phi_star; block[m:,m:] = R0
    return symmetrize(Mi.T @ block @ Mi), Phi_star, V

def check_fe(H_star, H0, C, Xi, gamma, label):
    S = hidden_stress(H_star, C, Xi)
    res = norm(S - Xi - (gamma/2)*(H_star - H0), 'fro')
    check(f"FE {label}", res < 1e-9, f"res={res:.3e}")
    return res

# ============================================================
print("=" * 60)
print("EDGE CASE 1: Xi = 0 (no perturbation)")
print("=" * 60)
for n in [2,3,5]:
    m = max(1, n//2)
    H0 = random_spd(n, 100+n)
    C = np.random.default_rng(200+n).standard_normal((m,n))
    Xi = np.zeros((n,n))
    gamma = 50.0
    H_star, Phi_star, V = build_solution(H0, C, Xi, gamma)
    check(f"H*=H0 (n={n})", norm(H_star - H0, 'fro') < 1e-12)
    check(f"V=0 (n={n})", norm(V, 'fro') < 1e-12)
    check_fe(H_star, H0, C, Xi, gamma, f"Xi=0,n={n}")

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 2: Xi = alpha*H0 (proportional)")
print("=" * 60)
for alpha in [0.1, 1.0, -0.5]:
    n, m, gamma = 4, 2, 50.0
    H0 = random_spd(n, 300)
    C = np.random.default_rng(400).standard_normal((m,n))
    Xi = alpha * H0
    H_star, Phi_star, V = build_solution(H0, C, Xi, gamma)
    L, Z, Phi0, R0 = split(H0, C)
    # V should be alpha*Phi0
    check(f"V=alpha*Phi0 (a={alpha})", norm(V - alpha*Phi0, 'fro') < 1e-10)
    # Phi* should be Phi0*(1 - 2*alpha/gamma)
    expected = Phi0 * (1 - 2*alpha/gamma)
    check(f"Phi*=scaled (a={alpha})", norm(Phi_star - expected, 'fro') < 1e-10)
    check_fe(H_star, H0, C, Xi, gamma, f"Xi=aH0,a={alpha}")

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 3: m = n-1 (almost full observation)")
print("=" * 60)
for n in [3, 5, 7]:
    m = n - 1
    H0 = random_spd(n, 500+n)
    C = np.random.default_rng(600+n).standard_normal((m,n))
    while np.linalg.matrix_rank(C) < m:
        C = np.random.default_rng(600+n+100).standard_normal((m,n))
    # Construct adapted Xi
    L, Z, Phi0, R0 = split(H0, C)
    Xi_raw = symmetrize(np.random.default_rng(700+n).standard_normal((n,n)))
    V = symmetrize(L.T @ Xi_raw @ L); U = symmetrize(Z.T @ Xi_raw @ Z)
    M = np.hstack([L, Z]); Mi = inv(M)
    block = np.zeros((n,n)); block[:m,:m] = V; block[m:,m:] = U
    Xi = symmetrize(Mi.T @ block @ Mi)
    gamma = 50.0
    H_star, Phi_star, _ = build_solution(H0, C, Xi, gamma)
    if eigvalsh(Phi_star)[0] > 0:
        check_fe(H_star, H0, C, Xi, gamma, f"m=n-1,n={n}")
        # R0 is scalar here
        check(f"R0 scalar (n={n})", R0.shape == (1,1))

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 4: m = 1 (minimal observation)")
print("=" * 60)
for n in [3, 5, 8]:
    m = 1
    H0 = random_spd(n, 800+n)
    C = np.random.default_rng(900+n).standard_normal((m,n))
    L, Z, Phi0, R0 = split(H0, C)
    Xi_raw = symmetrize(np.random.default_rng(1000+n).standard_normal((n,n)))
    V = symmetrize(L.T @ Xi_raw @ L); U = symmetrize(Z.T @ Xi_raw @ Z)
    M = np.hstack([L, Z]); Mi = inv(M)
    block = np.zeros((n,n)); block[:m,:m] = V; block[m:,m:] = U
    Xi = symmetrize(Mi.T @ block @ Mi)
    gamma = 50.0
    H_star, Phi_star, _ = build_solution(H0, C, Xi, gamma)
    if eigvalsh(Phi_star)[0] > 0:
        check_fe(H_star, H0, C, Xi, gamma, f"m=1,n={n}")
        # Phi0 is scalar
        check(f"Phi0 scalar (n={n})", Phi0.shape == (1,1))

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 5: Rank-1 Xi")
print("=" * 60)
for n in [3, 4, 6]:
    m = max(1, n//2)
    H0 = random_spd(n, 1100+n)
    C = np.random.default_rng(1200+n).standard_normal((m,n))
    xi = np.random.default_rng(1300+n).standard_normal(n)
    Xi_raw = np.outer(xi, xi)  # rank 1, PSD
    L, Z, Phi0, R0 = split(H0, C)
    V = symmetrize(L.T @ Xi_raw @ L); U = symmetrize(Z.T @ Xi_raw @ Z)
    M = np.hstack([L, Z]); Mi = inv(M)
    block = np.zeros((n,n)); block[:m,:m] = V; block[m:,m:] = U
    Xi = symmetrize(Mi.T @ block @ Mi)
    gamma = 100.0  # strong reg for rank-1
    H_star, Phi_star, Vout = build_solution(H0, C, Xi, gamma)
    if eigvalsh(Phi_star)[0] > 0:
        check_fe(H_star, H0, C, Xi, gamma, f"rank1,n={n}")
        # V should have rank <= 1
        v_eigs = eigvalsh(Vout)
        r = np.sum(np.abs(v_eigs) > 1e-10)
        check(f"V rank<=1 (n={n})", r <= 1, f"rank={r}")

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 6: gamma near SPD threshold")
print("=" * 60)
n, m = 3, 1
H0 = random_spd(n, 1400)
C = np.random.default_rng(1500).standard_normal((m,n))
L, Z, Phi0, R0 = split(H0, C)
Xi_raw = symmetrize(np.random.default_rng(1600).standard_normal((n,n)))
V = symmetrize(L.T @ Xi_raw @ L); U = symmetrize(Z.T @ Xi_raw @ Z)
M = np.hstack([L, Z]); Mi = inv(M)
block = np.zeros((n,n)); block[:m,:m] = V; block[m:,m:] = U
Xi = symmetrize(Mi.T @ block @ Mi)

# Find the critical gamma where Phi* becomes singular
# Phi* = Phi0 - (2/gamma)*V. For m=1, Phi0 and V are scalars.
phi0_val = float(Phi0); v_val = float(V)
if v_val > 0:
    gamma_crit = 2*v_val/phi0_val
    print(f"  Phi0 = {phi0_val:.4f}, V = {v_val:.4f}, gamma_crit = {gamma_crit:.4f}")

    for factor in [1.01, 1.1, 2.0, 10.0]:
        gamma = gamma_crit * factor
        Phi_star_val = phi0_val - 2*v_val/gamma
        print(f"  gamma={gamma:.4f} ({factor}x crit): Phi*={Phi_star_val:.6f}")
        if Phi_star_val > 0:
            H_star, _, _ = build_solution(H0, C, Xi, gamma)
            check_fe(H_star, H0, C, Xi, gamma, f"gamma={factor}x_crit")
        else:
            print(f"    Phi* <= 0, closed form invalid (expected)")

# ============================================================
print("\n" + "=" * 60)
print("EDGE CASE 7: H0 ill-conditioned")
print("=" * 60)
for cond in [10, 100, 1000, 10000]:
    n, m = 4, 2
    rng = np.random.default_rng(1700 + cond)
    Q = np.linalg.qr(rng.standard_normal((n,n)))[0]
    eigs = np.array([cond, 1.0, 1.0, 1.0/cond])
    H0 = Q @ np.diag(eigs) @ Q.T
    H0 = symmetrize(H0)
    C = rng.standard_normal((m,n))
    L, Z, Phi0, R0 = split(H0, C)
    Xi_raw = symmetrize(rng.standard_normal((n,n)))
    V = symmetrize(L.T @ Xi_raw @ L); U = symmetrize(Z.T @ Xi_raw @ Z)
    M = np.hstack([L, Z]); Mi = inv(M)
    block_xi = np.zeros((n,n)); block_xi[:m,:m] = V; block_xi[m:,m:] = U
    Xi = symmetrize(Mi.T @ block_xi @ Mi)
    gamma = 100.0 * cond  # scale with condition number
    H_star, Phi_star, _ = build_solution(H0, C, Xi, gamma)
    phi_eigs = eigvalsh(Phi_star)
    if phi_eigs[0] > 0:
        res = check_fe(H_star, H0, C, Xi, gamma, f"cond={cond}")
    else:
        print(f"  cond={cond}: Phi* not SPD (gamma too small)")

# ============================================================
print(f"\n{'='*60}")
print(f"  TOTAL: {P} passed, {F} failed")
print(f"{'='*60}")
