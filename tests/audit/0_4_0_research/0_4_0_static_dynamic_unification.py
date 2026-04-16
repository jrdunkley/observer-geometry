"""
The unification question: does the static optimum (closure_adapted_observer)
coincide with the dynamic optimum (maximum visible rate)?

For each molecule-solvent pair:
1. Compute the canonical observer (soft modes) vis_frac
2. Compute the closure-adapted observer (static optimum) vis_frac
3. Compute the random-search optimised observer (dynamic optimum) vis_frac
4. Compare: if adapted ~ optimised, the static and dynamic geometries agree.
"""

import sys
import time
import csv
import numpy as np
from numpy.linalg import inv, eigh, svd
from scipy.linalg import null_space
from pathlib import Path
from multiprocessing import Pool

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))

N_WORKERS = 20
DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
OUTPUT_DIR = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/papers/qm9_sweep_outputs")

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
    return {'label': row.get('label','?'), 'h_vib': hv, 'n_vib': hv.shape[0],
            'spd': float(np.min(ev))>1e-6}


def compute_vis_frac(H, C, Hdot):
    try:
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        dHinv = -Hinv @ Hdot @ Hinv
        dPhi = sym(-Phi @ (C @ dHinv @ C.T) @ Phi)
        vis_rate = float(np.trace(inv(Phi) @ dPhi))
        amb_rate = float(np.trace(Hinv @ Hdot))
        if abs(amb_rate) < 1e-10: return None, None, None
        V = (Hinv @ C.T @ Phi).T @ Hdot @ (Hinv @ C.T @ Phi)
        vm = float(np.min(eigh(V)[0]))
        return vis_rate/amb_rate, vis_rate, vm
    except:
        return None, None, None


def process_molecule(args):
    label, h_vac, h_sol, n_vib = args
    H = h_vac; Hdot = h_sol - h_vac
    m = min(3, n_vib - 1)
    if m < 1: return None

    try:
        Hinv = inv(H)
        amb_rate = float(np.trace(Hinv @ Hdot))
    except: return None
    if abs(amb_rate) < 0.01: return None

    # 1. Canonical observer
    C_canon = np.zeros((m, n_vib))
    for i in range(m): C_canon[i,i] = 1.0
    vf_canon, _, vm_canon = compute_vis_frac(H, C_canon, Hdot)

    # 2. Closure-adapted observer (static optimum: maximise visibility of Hdot through H)
    vf_adapted = None; vm_adapted = None
    try:
        from nomogeo import closure_adapted_observer
        result = closure_adapted_observer(H, [Hdot], m)
        C_adapted = result.C
        vf_adapted, _, vm_adapted = compute_vis_frac(H, C_adapted, Hdot)
    except:
        pass

    # 3. Dynamic optimum (random search + eigenvector observer)
    best_vf = vf_canon if vf_canon is not None else -999
    best_vm = vm_canon
    rng = np.random.default_rng(hash(label) % (2**31))

    for _ in range(50):
        Q = np.linalg.qr(rng.standard_normal((n_vib, n_vib)))[0]
        C_trial = Q[:m, :]
        vf, _, vm = compute_vis_frac(H, C_trial, Hdot)
        if vf is not None and vf > best_vf:
            best_vf = vf; best_vm = vm

    # Eigenvector observer: top m eigenvectors of H^{-1} Hdot
    try:
        P = Hinv @ Hdot
        _, evecs = eigh(sym(P))
        C_eig = evecs[:, -m:].T
        vf, _, vm = compute_vis_frac(H, C_eig, Hdot)
        if vf is not None and vf > best_vf:
            best_vf = vf; best_vm = vm
    except: pass

    return {
        'label': label, 'n_vib': n_vib, 'amb_rate': amb_rate,
        'vf_canon': vf_canon, 'vm_canon': vm_canon,
        'vf_adapted': vf_adapted, 'vm_adapted': vm_adapted,
        'vf_optimised': best_vf, 'vm_optimised': best_vm,
    }


def main():
    t0 = time.time()
    print("=" * 70)
    print("Static-Dynamic Unification: Does the adapted observer match the optimum?")
    print("=" * 70)

    vac_mols = {}; sol_mols = {}
    for shard in range(5):
        for row in iter_arrow_records("vacuum", shard, 500):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: vac_mols[mol['label']] = mol
            except: pass
    for shard in range(5):
        for row in iter_arrow_records("water", shard, 500):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: sol_mols[mol['label']] = mol
            except: pass

    matched = set(vac_mols.keys()) & set(sol_mols.keys())
    pairs = [(l, vac_mols[l]['h_vib'], sol_mols[l]['h_vib'], vac_mols[l]['n_vib'])
             for l in matched if vac_mols[l]['n_vib'] == sol_mols[l]['n_vib']]
    print(f"  Matched pairs: {len(pairs)}")

    with Pool(N_WORKERS) as pool:
        results = [r for r in pool.imap_unordered(process_molecule, pairs, chunksize=5) if r is not None]
    print(f"  Processed: {len(results)}")

    # Analyse
    has_all = [r for r in results
               if r['vf_canon'] is not None and r['vf_adapted'] is not None and r['vf_optimised'] is not None]

    print(f"  With all three observers: {len(has_all)}")

    if has_all:
        canon = np.array([r['vf_canon'] for r in has_all])
        adapted = np.array([r['vf_adapted'] for r in has_all])
        optimised = np.array([r['vf_optimised'] for r in has_all])

        print(f"\n  === Visible Fraction Comparison ===")
        print(f"  {'Observer':>20s}  {'Median':>8s}  {'Mean':>8s}  {'Std':>8s}")
        print(f"  {'Canonical':>20s}  {np.median(canon):8.4f}  {np.mean(canon):8.4f}  {np.std(canon):8.4f}")
        print(f"  {'Adapted (static)':>20s}  {np.median(adapted):8.4f}  {np.mean(adapted):8.4f}  {np.std(adapted):8.4f}")
        print(f"  {'Optimised (dynamic)':>20s}  {np.median(optimised):8.4f}  {np.mean(optimised):8.4f}  {np.std(optimised):8.4f}")

        # Correlation between adapted and optimised
        corr = np.corrcoef(adapted, optimised)[0,1]
        print(f"\n  Correlation(adapted, optimised): {corr:.4f}")

        # How often does the adapted observer beat canonical?
        adapted_beats_canon = np.sum(adapted > canon)
        optimised_beats_canon = np.sum(optimised > canon)
        adapted_near_optimised = np.sum(np.abs(adapted - optimised) < 0.1 * np.abs(optimised))

        print(f"  Adapted beats canonical: {adapted_beats_canon}/{len(has_all)} ({100*adapted_beats_canon/len(has_all):.1f}%)")
        print(f"  Optimised beats canonical: {optimised_beats_canon}/{len(has_all)} ({100*optimised_beats_canon/len(has_all):.1f}%)")
        print(f"  Adapted within 10% of optimised: {adapted_near_optimised}/{len(has_all)} ({100*adapted_near_optimised/len(has_all):.1f}%)")

        # V > 0 rates
        vm_canon = [r['vm_canon'] for r in has_all if r['vm_canon'] is not None]
        vm_adapted = [r['vm_adapted'] for r in has_all if r['vm_adapted'] is not None]
        vm_optimised = [r['vm_optimised'] for r in has_all if r['vm_optimised'] is not None]
        print(f"\n  V > 0 rates:")
        print(f"    Canonical: {sum(1 for v in vm_canon if v > 1e-6)}/{len(vm_canon)}")
        print(f"    Adapted:   {sum(1 for v in vm_adapted if v > 1e-6)}/{len(vm_adapted)}")
        print(f"    Optimised: {sum(1 for v in vm_optimised if v > 1e-6)}/{len(vm_optimised)}")

        # Rank correlation (are the same molecules best for both?)
        from scipy.stats import spearmanr
        sr, sp = spearmanr(adapted, optimised)
        print(f"\n  Spearman rank correlation(adapted, optimised): {sr:.4f} (p={sp:.2e})")

    t1 = time.time()
    print(f"\n  Elapsed: {t1-t0:.1f}s")

    # Write CSV
    with open(OUTPUT_DIR / "static_dynamic_comparison.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['label','n_vib','amb_rate','vf_canon','vf_adapted','vf_optimised',
                     'vm_canon','vm_adapted','vm_optimised'])
        for r in results:
            w.writerow([r['label'], r['n_vib'], r['amb_rate'],
                       r['vf_canon'], r['vf_adapted'], r['vf_optimised'],
                       r['vm_canon'], r['vm_adapted'], r['vm_optimised']])


if __name__ == "__main__":
    main()
