"""
Part 6 only — lightweight stationarity verification.
10 trials, 10 restarts, n in {3,4}, 20 workers.
"""
import sys, os, time
import numpy as np
from numpy.linalg import det, inv, norm
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
from nomogeo.validation import symmetrize
from scipy.linalg import null_space

WORKERS = 20

def random_spd(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return symmetrize(A @ A.T + np.eye(n))

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

def vis_rate(H, C, Hdot):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Hdot @ Z)
    return float(np.trace(np.linalg.solve(H, Hdot)) - np.trace(np.linalg.solve(R, U_h)))

def grad_H(H, C, Hdot):
    Z = hidden_frame(C)
    Hi = inv(H); R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Hdot @ Z); Ri = inv(R)
    return symmetrize(-Hi @ Hdot @ Hi + Z @ Ri @ U_h @ Ri @ Z.T)

def hidden_stress(H, C, Hdot):
    Z = hidden_frame(C)
    R = symmetrize(Z.T @ H @ Z)
    U_h = symmetrize(Z.T @ Hdot @ Z)
    A = H @ Z @ inv(R)
    return symmetrize(A @ U_h @ A.T)

def cholesky_params(H):
    L = np.linalg.cholesky(H); p = []
    for i in range(H.shape[0]):
        for j in range(i+1):
            p.append(np.log(L[i,j]) if i==j else L[i,j])
    return np.array(p)

def params_to_H(params, n):
    L = np.zeros((n,n)); idx=0
    for i in range(n):
        for j in range(i+1):
            L[i,j] = params[idx]; idx+=1
    for i in range(n): L[i,i] = np.exp(np.clip(L[i,i],-10,10))
    return L @ L.T

def worker(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(3, 5))  # only n=3 or 4
    m = int(rng.integers(1, n))
    H_init = random_spd(n, seed+1000)
    C = random_C(m, n, seed+2000)
    Hdot = random_sym(n, seed+3000)
    det_tgt = det(H_init)
    if det_tgt <= 0: return {'n':n,'m':m,'pass':False,'reason':'bad det'}

    def obj(params):
        H = params_to_H(params, n)
        d = det(H)
        if d <= 0 or np.isnan(d): return 1e10
        pen = 500.0 * (np.log(d) - np.log(det_tgt))**2
        try: vr = vis_rate(H, C, Hdot)
        except: return 1e10
        return -vr + pen if np.isfinite(vr) else 1e10

    best_val, best_p = np.inf, None
    for restart in range(10):
        Hr = random_spd(n, seed+10000+restart)
        sc = (det_tgt / det(Hr))**(1.0/(2*n))
        Hr = Hr * sc
        try: p0 = cholesky_params(Hr)
        except: continue
        res = minimize(obj, p0, method='Nelder-Mead',
                       options={'maxiter':10000, 'xatol':1e-12, 'fatol':1e-12})
        if res.fun < best_val: best_val, best_p = res.fun, res.x

    if best_p is None: return {'n':n,'m':m,'pass':False,'reason':'no conv'}
    Hs = params_to_H(best_p, n)
    g = grad_H(Hs, C, Hdot)
    Hi = inv(Hs)
    lam = np.trace(Hs @ g) / n
    station = norm(g - lam * Hi, 'fro')
    vr0 = vis_rate(H_init, C, Hdot)
    vrs = vis_rate(Hs, C, Hdot)

    # Also check tensor field equation at optimum
    S = hidden_stress(Hs, C, Hdot)
    fe_res = norm(S - Hdot - lam * Hs, 'fro')

    return {
        'n':n, 'm':m,
        'vr_init': round(vr0, 4), 'vr_star': round(vrs, 4),
        'lambda': round(lam, 6),
        'stationarity': round(station, 6),
        'fe_tensor_res': round(fe_res, 6),
        'pass': station < 0.01,
    }

if __name__ == '__main__':
    seeds = list(range(500, 510))
    print(f"Part 6: Optimise H*, 10 trials, 10 restarts each, n∈{{3,4}}, {WORKERS} workers")
    t0 = time.time()
    results = [None]*len(seeds)
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(worker, s): i for i, s in enumerate(seeds)}
        done = 0
        for fut in as_completed(futs):
            results[futs[fut]] = fut.result()
            done += 1
            bar = '█'*(done*40//10) + '░'*(40-done*40//10)
            print(f"\r  [{bar}] {done}/10", end='', flush=True)
    print(f"  ({time.time()-t0:.1f}s)\n")

    for i, r in enumerate(results):
        tag = "PASS" if r.get('pass') else "FAIL"
        print(f"  [{tag}] Trial {i}: {r}")

    passed = sum(1 for r in results if r.get('pass'))
    print(f"\n  => {passed}/10 passed")

    conv = [r for r in results if 'stationarity' in r]
    if conv:
        stats = [r['stationarity'] for r in conv]
        fe_res = [r['fe_tensor_res'] for r in conv]
        print(f"  Stationarity: median={np.median(stats):.3e}, max={np.max(stats):.3e}")
        print(f"  FE tensor residual: median={np.median(fe_res):.3e}, max={np.max(fe_res):.3e}")
        print(f"  (FE tensor: if large, S=Hdot+lam*H is NOT an identity but holds at optima)")
