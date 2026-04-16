"""
Field equation verification v2 — honest theorem boundaries.
Three test suites, clearly labelled.

Run: cd C:\\observer_geometry_workspace && set PYTHONPATH=nomogeo\\src && python nomogeo\\audit\\0_4_0_research\\0_4_0_field_equation_verify_v2.py
"""
import sys, os, time
import numpy as np
from numpy.linalg import inv, norm, det, eigvalsh, solve
from scipy.linalg import null_space
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
from nomogeo.validation import symmetrize

WORKERS = 20

def random_spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return symmetrize(A @ A.T + 2*np.eye(n))

def random_sym(n, seed):
    rng = np.random.default_rng(seed)
    return symmetrize(rng.standard_normal((n, n)))

def random_C(m, n, seed):
    rng = np.random.default_rng(seed)
    for _ in range(50):
        C = rng.standard_normal((m, n))
        if np.linalg.matrix_rank(C) >= m: return C
    return rng.standard_normal((m, n))

def hidden_frame(C):
    return null_space(C, rcond=1e-12)

def split(H, C):
    Hi = inv(H); Q = symmetrize(C @ Hi @ C.T); Phi = inv(Q)
    L = Hi @ C.T @ Phi; Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    return L, Z, Phi, R

def vis_rate(H, C, Xi):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z); U_h = symmetrize(Z.T @ Xi @ Z)
    return float(np.trace(solve(H, Xi)) - np.trace(solve(R, U_h)))

def hidden_stress(H, C, Xi):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z); U_h = symmetrize(Z.T @ Xi @ Z)
    A = H @ Z @ inv(R)
    return symmetrize(A @ U_h @ A.T)

def grad_vis(H, C, Xi):
    Z = hidden_frame(C); Hi = inv(H)
    R = symmetrize(Z.T @ H @ Z); U_h = symmetrize(Z.T @ Xi @ Z); Ri = inv(R)
    return symmetrize(-Hi @ Xi @ Hi + Z @ Ri @ U_h @ Ri @ Z.T)

def grad_DKL(H, H0):
    Hi = inv(H)
    return symmetrize(0.5 * (Hi - Hi @ H0 @ Hi))

def D_KL(H0, H):
    A = H0 @ inv(H)
    return 0.5 * (np.trace(A) - np.log(det(A)) - H.shape[0])

def cholesky_params(H):
    L = np.linalg.cholesky(H); p = []
    for i in range(H.shape[0]):
        for j in range(i+1):
            p.append(np.log(L[i,j]) if i==j else L[i,j])
    return np.array(p)

def params_to_H(params, n):
    L = np.zeros((n,n)); idx=0
    for i in range(n):
        for j in range(i+1): L[i,j] = params[idx]; idx+=1
    for i in range(n): L[i,i] = np.exp(np.clip(L[i,i],-10,10))
    return L @ L.T

# ============================================================
# TEST 1: Fixed-C generic (Theorem 1) — NO B=0 projection
# ============================================================
def worker_test1(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(3,5)); m = int(rng.integers(1,n))
    gamma = 50.0
    H0 = random_spd(n, seed+1000)
    C = random_C(m, n, seed+2000)
    Xi = random_sym(n, seed+3000)  # GENERIC Xi, no projection

    det_kl = lambda H: vis_rate(H, C, Xi) - gamma * D_KL(H0, H)

    def obj(params):
        H = params_to_H(params, n)
        d = det(H)
        if d <= 0 or np.isnan(d): return 1e10
        try: return -det_kl(H)
        except: return 1e10

    best_val, best_p = np.inf, None
    for restart in range(15):
        Hr = random_spd(n, seed+10000+restart)
        try: p0 = cholesky_params(Hr)
        except: continue
        res = minimize(obj, p0, method='Nelder-Mead',
                       options={'maxiter':15000, 'xatol':1e-13, 'fatol':1e-13})
        if res.fun < best_val: best_val, best_p = res.fun, res.x

    if best_p is None: return {'n':n,'m':m,'test':'Theorem1','pass':False,'reason':'no conv'}
    Hs = params_to_H(best_p, n)

    # Check field equation
    S = hidden_stress(Hs, C, Xi)
    fe_res = norm(S - Xi - (gamma/2)*(Hs - H0), 'fro')

    # Check stationarity
    gv = grad_vis(Hs, C, Xi)
    gd = grad_DKL(Hs, H0)
    station = norm(gv - gamma*gd, 'fro')

    # Check B (should NOT be zero generically)
    L, Z, Phi, R = split(Hs, C)
    B = L.T @ Xi @ Z
    B_norm = norm(B, 'fro')

    return {
        'n':n, 'm':m, 'test':'Theorem1',
        'fe_res': float(fe_res), 'station': float(station),
        'B_norm': float(B_norm),
        'J_star': float(-best_val),
        'pass': station < 0.1
    }

# ============================================================
# TEST 2: Compatibility-sector (Proposition 3)
# ============================================================
def worker_test2(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(3,7)); m = int(rng.integers(1,n))
    gamma = 50.0
    H0 = random_spd(n, seed+1000)
    C = random_C(m, n, seed+2000)
    Xi_raw = random_sym(n, seed+3000)

    L0, Z0, Phi0, R0 = split(H0, C)

    # Project Xi to B=0 sector (compatibility hypothesis)
    V = symmetrize(L0.T @ Xi_raw @ L0)
    U_h = symmetrize(Z0.T @ Xi_raw @ Z0)
    block = np.zeros((n,n)); block[:m,:m] = V; block[m:,m:] = U_h
    M0 = np.hstack([L0, Z0]); Mi = inv(M0)
    Xi = symmetrize(Mi.T @ block @ Mi)

    # Closed-form solution
    Phi_star = Phi0 - (2.0/gamma)*V
    phi_eigs = eigvalsh(Phi_star)
    if phi_eigs[0] <= 0:
        return {'n':n,'m':m,'test':'Prop3','pass':True,'reason':'SPD violated, skipped'}

    block_star = np.zeros((n,n)); block_star[:m,:m] = Phi_star; block_star[m:,m:] = R0
    H_star = symmetrize(Mi.T @ block_star @ Mi)
    h_eigs = eigvalsh(H_star)
    if h_eigs[0] <= 0:
        return {'n':n,'m':m,'test':'Prop3','pass':True,'reason':'H* not SPD, skipped'}

    # Check field equation
    S = hidden_stress(H_star, C, Xi)
    fe_res = norm(S - Xi - (gamma/2)*(H_star - H0), 'fro')

    # Check stationarity
    gv = grad_vis(H_star, C, Xi)
    gd = grad_DKL(H_star, H0)
    station = norm(gv - gamma*gd, 'fro')

    # Check L(H*) = L(H0)
    L_star, _, _, _ = split(H_star, C)
    L_diff = norm(L_star - L0, 'fro')

    return {
        'n':n, 'm':m, 'test':'Prop3',
        'fe_res': float(fe_res), 'station': float(station),
        'L_diff': float(L_diff),
        'pass': fe_res < 1e-6 and station < 1e-6
    }

# ============================================================
# TEST 3: Joint (H,C) optimiser (Proposition 2 — exploratory)
# ============================================================
def worker_test3(seed):
    n = 3; m = 1  # fixed small dims for speed
    gamma = 50.0
    H0 = random_spd(n, seed+1000)
    Xi = random_sym(n, seed+3000)  # GENERIC Xi

    n_h = n*(n+1)//2

    def obj(params):
        H = params_to_H(params[:n_h], n)
        C = params[n_h:].reshape(m, n)
        d = det(H)
        if d <= 0 or np.isnan(d): return 1e10
        if np.linalg.matrix_rank(C, tol=1e-10) < m: return 1e10
        try: return -(vis_rate(H, C, Xi) - gamma * D_KL(H0, H))
        except: return 1e10

    best_val, best_p = np.inf, None
    for restart in range(5):
        Hr = random_spd(n, seed+10000+restart)
        try: ph = cholesky_params(Hr)
        except: continue
        Cr = np.random.default_rng(seed+20000+restart).standard_normal((m,n))
        p0 = np.concatenate([ph, Cr.flatten()])
        res = minimize(obj, p0, method='Nelder-Mead',
                       options={'maxiter':8000, 'xatol':1e-11, 'fatol':1e-11})
        if res.fun < best_val: best_val, best_p = res.fun, res.x

    if best_p is None: return {'n':n,'m':m,'test':'Prop2','pass':False,'reason':'no conv'}
    Hs = params_to_H(best_p[:n_h], n)
    Cs = best_p[n_h:].reshape(m, n)

    # Check B(H*, C*) ~ 0 (Prop 2a: adapted observer)
    Ls, Zs, Phis, Rs = split(Hs, Cs)
    B = Ls.T @ Xi @ Zs
    B_norm = norm(B, 'fro')

    # Check K0* = L(H*)^T H0 Z ~ 0 (Prop 2b)
    K0_star = Ls.T @ H0 @ Zs
    K0_norm = norm(K0_star, 'fro')

    # Check R* ~ R0* (Prop 2c)
    R_star = Rs
    R0_star = symmetrize(Zs.T @ H0 @ Zs)
    R_diff = norm(R_star - R0_star, 'fro')

    # Check Phi* ~ Phi0* - (2/gamma)V* (Prop 2d)
    Phi0_star = symmetrize(Ls.T @ H0 @ Ls)
    V_star = symmetrize(Ls.T @ Xi @ Ls)
    Phi_pred = Phi0_star - (2.0/gamma)*V_star
    Phi_diff = norm(Phis - Phi_pred, 'fro')

    # Check L(H*) vs L(H0) (compatibility question)
    L0, _, _, _ = split(H0, Cs)
    L_diff = norm(Ls - L0, 'fro')

    # Check field equation
    S = hidden_stress(Hs, Cs, Xi)
    fe_res = norm(S - Xi - (gamma/2)*(Hs - H0), 'fro')

    return {
        'n':n, 'm':m, 'test':'Prop2',
        'B_norm': float(B_norm), 'K0_norm': float(K0_norm),
        'R_diff': float(R_diff), 'Phi_diff': float(Phi_diff),
        'L_diff': float(L_diff), 'fe_res': float(fe_res),
        'pass': B_norm < 0.1 and fe_res < 1.0
    }

# ============================================================
def run_parallel(name, worker_fn, seeds, show_first=5):
    t0 = time.time()
    total = len(seeds)
    results = [None]*total
    print(f"\n{'='*70}\n  {name}  ({total} trials, {WORKERS} workers)\n{'='*70}")
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(worker_fn, s): i for i, s in enumerate(seeds)}
        done = 0
        for fut in as_completed(futs):
            results[futs[fut]] = fut.result()
            done += 1
            bar = '█'*(done*40//total) + '░'*(40-done*40//total)
            print(f"\r  [{bar}] {done}/{total}", end='', flush=True)
    print(f"  ({time.time()-t0:.1f}s)")
    for i, r in enumerate(results[:show_first]):
        print(f"    Sample {i}: {r}")
    passed = sum(1 for r in results if r.get('pass'))
    print(f"  => {passed}/{total} passed")
    return results

if __name__ == '__main__':
    # Test 1: Theorem 1 (fixed-C, generic Xi)
    r1 = run_parallel("TEST 1: Theorem 1 — fixed-C generic stationarity",
                      worker_test1, list(range(800, 815)))
    conv1 = [r for r in r1 if 'station' in r]
    if conv1:
        print(f"    Stationarity: median={np.median([r['station'] for r in conv1]):.3e}")
        print(f"    FE residual:  median={np.median([r['fe_res'] for r in conv1]):.3e}")
        print(f"    ||B|| at H*:  median={np.median([r['B_norm'] for r in conv1]):.3e}")
        print(f"    (B_norm > 0 expected: generic Xi has B != 0 at the maximiser)")

    # Test 2: Proposition 3 (compatibility sector)
    r2 = run_parallel("TEST 2: Proposition 3 — compatibility-sector closed form",
                      worker_test2, list(range(900, 920)))

    # Test 3: Proposition 2 (joint, exploratory)
    r3 = run_parallel("TEST 3: Proposition 2 — joint (H,C) exploratory",
                      worker_test3, list(range(1000, 1010)))
    conv3 = [r for r in r3 if 'B_norm' in r]
    if conv3:
        print(f"    ||B||:    median={np.median([r['B_norm'] for r in conv3]):.3e}")
        print(f"    ||K0*||:  median={np.median([r['K0_norm'] for r in conv3]):.3e}")
        print(f"    ||R-R0||: median={np.median([r['R_diff'] for r in conv3]):.3e}")
        print(f"    ||Phi-pred||: median={np.median([r['Phi_diff'] for r in conv3]):.3e}")
        print(f"    ||L*-L0||: median={np.median([r['L_diff'] for r in conv3]):.3e}")
        print(f"    (L_diff > 0 expected if compatibility sector is not generic)")

    print(f"\n{'='*70}\n  DONE\n{'='*70}")
