"""
Deep push: Grassmannian optimisation, molecule-by-molecule comparison,
Pareto frontier, m-dependence, capture curve universality, nomoselect robustness.

This is the comprehensive experiment that settles the unification question.
"""

import sys
import time
import csv
import json
import numpy as np
from numpy.linalg import inv, eigh, norm, svd, slogdet
from scipy.linalg import null_space, sqrtm, expm
from scipy.optimize import minimize
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))

N_WORKERS = 20
DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
OUTPUT_DIR = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/papers/qm9_deep_push")

ATOMIC_MASS = {1: 1.00782503223, 6: 12.0, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273}


def sym(M): return 0.5*(M+M.T)

def iter_arrow_records(split, shard, limit=None):
    import pyarrow.ipc as ipc
    path = DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow"
    if not path.exists(): return
    seen = 0
    with ipc.open_stream(path) as reader:
        for batch in reader:
            for row in batch.to_pylist():
                yield row
                seen += 1
                if limit and seen >= limit: return

def vibrational_hessian(row):
    an = np.asarray(row["atomic_numbers"], dtype=int)
    pos = np.asarray(row["positions"], dtype=float)
    na = len(an)
    masses = np.array([ATOMIC_MASS[int(z)] for z in an])
    h_raw = np.asarray(row["hessian"], dtype=float)
    if h_raw.shape == (na,3,na,3): hc = sym(h_raw.reshape(3*na,3*na))
    else: hc = sym(h_raw.transpose(0,2,1,3).reshape(3*na,3*na))
    dm = np.repeat(masses,3); ism = 1.0/np.sqrt(dm)
    hmw = sym((ism[:,None]*hc)*ism[None,:])
    w = np.sqrt(masses); ctr = np.average(pos,axis=0,weights=masses); ctrd = pos-ctr
    cols = []
    for ax in range(3):
        v = np.zeros((na,3)); v[:,ax]=w; cols.append(v.reshape(-1))
    for ax in np.eye(3):
        v = np.cross(ax[None,:],ctrd)*w[:,None]; cols.append(v.reshape(-1))
    raw = np.column_stack(cols)
    u,sv,_ = np.linalg.svd(raw,full_matrices=False)
    rk = int(np.sum(sv>1e-10*max(1,float(np.max(sv)))))
    vb = null_space(u[:,:rk].T)
    hv = sym(vb.T @ hmw @ vb)
    ev = np.linalg.eigvalsh(hv)
    counts = {}
    for z in an: counts[int(z)] = counts.get(int(z),0)+1
    sm = {1:"H",6:"C",7:"N",8:"O",9:"F"}
    formula = "".join(sm.get(z,"?")+str(c) for z,c in sorted(counts.items()))
    return {'label': row.get('label','?'), 'h_vib': hv, 'n_vib': hv.shape[0],
            'spd': float(np.min(ev))>1e-6, 'formula': formula, 'n_atoms': na,
            'min_eig': float(np.min(ev)), 'max_eig': float(np.max(ev)),
            'condition': float(np.max(ev)/max(np.min(ev),1e-15))}


# ---- Core computation functions ----

def compute_diagnostics(H, C, Hdot):
    """Full diagnostic for one (H, C, Hdot) triple."""
    n = H.shape[0]; m = C.shape[0]
    try:
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C); Z = Vt[m:].T
        R = Z.T @ H @ Z; Rinv = inv(R)

        V = L.T @ Hdot @ L
        B = L.T @ Hdot @ Z
        U_h = Z.T @ Hdot @ Z

        dHinv = -Hinv @ Hdot @ Hinv
        dPhi = sym(-Phi @ (C @ dHinv @ C.T) @ Phi)
        vis_rate = float(np.trace(inv(Phi) @ dPhi))
        hid_rate = float(np.trace(Rinv @ U_h))
        amb_rate = float(np.trace(Hinv @ Hdot))

        v_eigvals = eigh(V)[0]
        v_min = float(np.min(v_eigvals))

        # Closure leakage (commutator norm of Hdot with the projector)
        P = L @ C  # projector
        commutator = Hdot @ P - P @ Hdot
        leakage = float(norm(commutator, 'fro'))

        # Qhat trace (hidden coupling strength)
        Qhat = B @ Rinv @ B.T
        qhat_tr = float(np.trace(Qhat))

        vis_frac = vis_rate / amb_rate if abs(amb_rate) > 1e-10 else float('nan')

        return {
            'vis_rate': vis_rate, 'hid_rate': hid_rate, 'amb_rate': amb_rate,
            'vis_frac': vis_frac, 'v_min': v_min,
            'phi_min': float(np.min(eigh(Phi)[0])),
            'leakage': leakage, 'qhat_trace': qhat_tr,
            'cons_err': abs(vis_rate + hid_rate - amb_rate),
        }
    except:
        return None


def grassmann_optimise(H, Hdot, m, objective='vis_frac', n_restarts=10, max_iter=200):
    """
    Proper Riemannian optimisation on Gr(m,n).

    Parameterise the observer as C = Q[:m, :] where Q = exp(A) @ Q_init
    with A skew-symmetric. Optimise A via projected gradient.
    """
    n = H.shape[0]
    Hinv = inv(H)
    amb_rate = float(np.trace(Hinv @ Hdot))
    if abs(amb_rate) < 1e-10:
        return None, None

    def neg_objective(A_flat, Q_init):
        A = np.zeros((n, n))
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                A[i,j] = A_flat[idx]
                A[j,i] = -A_flat[idx]
                idx += 1
        Q = expm(A) @ Q_init
        C = Q[:m, :]
        diag = compute_diagnostics(H, C, Hdot)
        if diag is None: return 1e10
        if objective == 'vis_frac':
            vf = diag['vis_frac']
            return -vf if not np.isnan(vf) else 1e10
        elif objective == 'v_min':
            return -diag['v_min']
        elif objective == 'combined':
            # Weighted combination of vis_frac and v_min
            vf = diag['vis_frac'] if not np.isnan(diag['vis_frac']) else 0
            vm = diag['v_min']
            return -(vf + 0.5 * max(vm, 0))
        return 1e10

    n_params = n * (n - 1) // 2
    rng = np.random.default_rng(42)

    best_val = 1e10
    best_C = None
    best_diag = None

    for restart in range(n_restarts):
        Q_init = np.linalg.qr(rng.standard_normal((n, n)))[0]
        A0 = rng.standard_normal(n_params) * 0.1

        try:
            result = minimize(neg_objective, A0, args=(Q_init,),
                            method='L-BFGS-B', options={'maxiter': max_iter, 'ftol': 1e-12})
            if result.fun < best_val:
                best_val = result.fun
                # Reconstruct C
                A = np.zeros((n, n)); idx = 0
                for i in range(n):
                    for j in range(i+1, n):
                        A[i,j] = result.x[idx]; A[j,i] = -result.x[idx]; idx += 1
                Q = expm(A) @ Q_init
                best_C = Q[:m, :].copy()
                best_diag = compute_diagnostics(H, best_C, Hdot)
        except:
            pass

    return best_C, best_diag


def process_molecule_deep(args):
    """Full deep analysis for one molecule."""
    label, h_vac, h_sol, n_vib, formula, n_atoms, condition = args
    H = h_vac; Hdot = h_sol - h_vac

    try:
        Hinv = inv(H)
        amb_rate = float(np.trace(Hinv @ Hdot))
    except: return None
    if abs(amb_rate) < 0.01: return None

    results = {'label': label, 'formula': formula, 'n_vib': n_vib,
               'n_atoms': n_atoms, 'condition': condition, 'amb_rate': amb_rate}

    # ---- Test across m values (m = 1, 2, 3, 5, n//4, n//2) ----
    m_values = sorted(set([1, 2, 3, 5, max(1, n_vib//4), max(1, n_vib//2)]))
    m_values = [m for m in m_values if m < n_vib]

    m_results = []
    for m in m_values:
        # Canonical
        C_canon = np.zeros((m, n_vib))
        for i in range(m): C_canon[i,i] = 1.0
        d_canon = compute_diagnostics(H, C_canon, Hdot)

        # Adapted (static)
        d_adapted = None
        try:
            from nomogeo import closure_adapted_observer
            res = closure_adapted_observer(H, [Hdot], m)
            d_adapted = compute_diagnostics(H, res.C, Hdot)
        except: pass

        # Grassmannian optimised (proper, vis_frac objective)
        C_opt_vf, d_opt_vf = grassmann_optimise(H, Hdot, m, 'vis_frac', n_restarts=5, max_iter=100)

        # Grassmannian optimised (V > 0 objective)
        C_opt_vm, d_opt_vm = grassmann_optimise(H, Hdot, m, 'v_min', n_restarts=5, max_iter=100)

        m_results.append({
            'm': m,
            'canon': d_canon,
            'adapted': d_adapted,
            'opt_vf': d_opt_vf,
            'opt_vm': d_opt_vm,
        })

    results['m_results'] = m_results
    return results


def load_pairs(solvent, per_shard=500):
    vac_mols = {}; sol_mols = {}
    for shard in range(5):
        for row in iter_arrow_records("vacuum", shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: vac_mols[mol['label']] = mol
            except: pass
    for shard in range(5):
        for row in iter_arrow_records(solvent, shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: sol_mols[mol['label']] = mol
            except: pass
    matched = set(vac_mols.keys()) & set(sol_mols.keys())
    pairs = []
    for l in matched:
        v = vac_mols[l]; s = sol_mols[l]
        if v['n_vib'] == s['n_vib']:
            pairs.append((l, v['h_vib'], s['h_vib'], v['n_vib'], v['formula'],
                         v['n_atoms'], v['condition']))
    return pairs


def main():
    t0 = time.time()
    print("=" * 70)
    print("Deep Push: Grassmannian Optimisation + Full Analysis")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data for all three solvents
    all_results = {}
    for solvent in ["water", "thf", "toluene"]:
        print(f"\n  Loading {solvent} pairs...")
        pairs = load_pairs(solvent, per_shard=200)  # 200/shard for speed with Grassmann opt
        print(f"  Matched: {len(pairs)}")

        if not pairs:
            continue

        print(f"  Processing ({N_WORKERS} workers, Grassmannian opt)...")
        with Pool(N_WORKERS) as pool:
            results = [r for r in pool.imap_unordered(process_molecule_deep, pairs, chunksize=2)
                      if r is not None]
        print(f"  Processed: {len(results)}")
        all_results[solvent] = results

    # ---- Analysis ----
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    for solvent, results in all_results.items():
        print(f"\n  === {solvent.upper()} ({len(results)} molecules) ===")

        # Per-m analysis
        for target_m in [1, 2, 3, 5]:
            canon_vf = []; adapted_vf = []; opt_vf_vf = []; opt_vm_vf = []
            canon_vm = []; adapted_vm = []; opt_vf_vm = []; opt_vm_vm = []
            canon_leak = []; adapted_leak = []

            for r in results:
                for mr in r['m_results']:
                    if mr['m'] != target_m: continue
                    if mr['canon']:
                        vf = mr['canon']['vis_frac']
                        if not np.isnan(vf) and abs(vf) < 100:
                            canon_vf.append(vf)
                            canon_vm.append(mr['canon']['v_min'])
                            canon_leak.append(mr['canon']['leakage'])
                    if mr['adapted']:
                        vf = mr['adapted']['vis_frac']
                        if not np.isnan(vf) and abs(vf) < 100:
                            adapted_vf.append(vf)
                            adapted_vm.append(mr['adapted']['v_min'])
                            adapted_leak.append(mr['adapted']['leakage'])
                    if mr['opt_vf']:
                        vf = mr['opt_vf']['vis_frac']
                        if not np.isnan(vf) and abs(vf) < 100:
                            opt_vf_vf.append(vf)
                            opt_vf_vm.append(mr['opt_vf']['v_min'])
                    if mr['opt_vm']:
                        vf = mr['opt_vm']['vis_frac']
                        if not np.isnan(vf) and abs(vf) < 100:
                            opt_vm_vf.append(vf)
                            opt_vm_vm.append(mr['opt_vm']['v_min'])

            if not canon_vf: continue

            print(f"\n  m = {target_m}:")
            print(f"  {'Observer':>25s}  {'med_vf':>8s}  {'med_vm':>8s}  {'V>0':>6s}  {'leak':>8s}  {'n':>4s}")

            def row(name, vfs, vms, leaks=None):
                n = len(vfs)
                if n == 0: return
                mvf = np.median(vfs)
                mvm = np.median(vms) if vms else float('nan')
                vpos = sum(1 for v in vms if v > 1e-6) if vms else 0
                mlk = f"{np.median(leaks):.4f}" if leaks else "N/A"
                print(f"  {name:>25s}  {mvf:8.4f}  {mvm:8.4f}  {vpos:>4d}/{n:<1d}  {mlk:>8s}  {n:4d}")

            row("Canonical", canon_vf, canon_vm, canon_leak)
            row("Adapted (static)", adapted_vf, adapted_vm, adapted_leak)
            row("Grassmann opt (vis_frac)", opt_vf_vf, opt_vf_vm)
            row("Grassmann opt (V>0)", opt_vm_vf, opt_vm_vm)

            # Correlation between adapted and Grassmann-optimised
            if adapted_vf and opt_vf_vf and len(adapted_vf) == len(opt_vf_vf):
                from scipy.stats import spearmanr
                sr, sp = spearmanr(adapted_vf, opt_vf_vf)
                print(f"    Spearman(adapted, opt_vf): {sr:.4f} (p={sp:.2e})")

        # ---- Molecule-by-molecule divergence analysis (at m=3) ----
        print(f"\n  Molecule-by-molecule divergence (m=3):")
        divergences = []
        for r in results:
            for mr in r['m_results']:
                if mr['m'] != 3: continue
                if mr['canon'] and mr['adapted'] and mr['opt_vf']:
                    vf_c = mr['canon']['vis_frac']
                    vf_a = mr['adapted']['vis_frac']
                    vf_o = mr['opt_vf']['vis_frac']
                    if all(not np.isnan(x) and abs(x) < 100 for x in [vf_c, vf_a, vf_o]):
                        divergences.append({
                            'label': r['label'], 'formula': r['formula'],
                            'n_vib': r['n_vib'], 'n_atoms': r['n_atoms'],
                            'condition': r['condition'],
                            'vf_canon': vf_c, 'vf_adapted': vf_a, 'vf_optimised': vf_o,
                            'delta_adapted': vf_a - vf_c,
                            'delta_opt': vf_o - vf_c,
                            'adapted_vs_opt': vf_a - vf_o,
                            'vm_canon': mr['canon']['v_min'],
                            'vm_adapted': mr['adapted']['v_min'],
                            'vm_opt': mr['opt_vf']['v_min'],
                        })

        if divergences:
            # Where does adapted diverge most from optimised?
            divergences.sort(key=lambda d: -abs(d['adapted_vs_opt']))
            print(f"  Largest adapted-vs-optimised gaps:")
            print(f"  {'label':>18s}  {'n':>3s}  {'cond':>8s}  {'canon':>7s}  {'adapt':>7s}  {'opt':>7s}  {'gap':>7s}")
            for d in divergences[:10]:
                print(f"  {d['label']:>18s}  {d['n_vib']:3d}  {d['condition']:8.1f}  "
                      f"{d['vf_canon']:7.3f}  {d['vf_adapted']:7.3f}  {d['vf_optimised']:7.3f}  "
                      f"{d['adapted_vs_opt']:7.3f}")

            # What molecular features predict divergence?
            gaps = np.array([d['adapted_vs_opt'] for d in divergences])
            n_vibs = np.array([d['n_vib'] for d in divergences])
            conds = np.array([d['condition'] for d in divergences])
            n_atoms_arr = np.array([d['n_atoms'] for d in divergences])

            from scipy.stats import spearmanr
            sr_nv, _ = spearmanr(np.abs(gaps), n_vibs)
            sr_cond, _ = spearmanr(np.abs(gaps), conds)
            sr_na, _ = spearmanr(np.abs(gaps), n_atoms_arr)
            print(f"\n  Divergence predictors:")
            print(f"    Spearman(|gap|, n_vib): {sr_nv:.3f}")
            print(f"    Spearman(|gap|, condition): {sr_cond:.3f}")
            print(f"    Spearman(|gap|, n_atoms): {sr_na:.3f}")

            # Write divergence CSV
            with open(OUTPUT_DIR / f"divergence_{solvent}.csv", 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=list(divergences[0].keys()))
                w.writeheader()
                for d in divergences: w.writerow(d)

    # ---- Capture curve across solvents ----
    print(f"\n  === Capture Curve Universality ===")
    for solvent, results in all_results.items():
        ratio_bins = defaultdict(list)
        for r in results:
            n = r['n_vib']
            for mr in r['m_results']:
                if mr['canon'] and not np.isnan(mr['canon']['vis_frac']) and abs(mr['canon']['vis_frac']) < 100:
                    ratio = mr['m'] / n
                    bin_idx = min(int(ratio * 10), 9)
                    ratio_bins[bin_idx].append(mr['canon']['vis_frac'])

        half_point = None
        print(f"\n  {solvent}:")
        print(f"  {'m/n':>6s}  {'n':>5s}  {'median':>8s}")
        for i in range(10):
            if ratio_bins[i]:
                mid = (i+0.5)/10
                med = np.median(ratio_bins[i])
                print(f"  {mid:6.2f}  {len(ratio_bins[i]):5d}  {med:8.4f}")
                if med >= 0.5 and half_point is None:
                    half_point = mid
        if half_point:
            print(f"  50% capture at m/n ~ {half_point:.2f}")

    t1 = time.time()
    print(f"\n  Total elapsed: {t1-t0:.1f}s")


if __name__ == "__main__":
    main()
