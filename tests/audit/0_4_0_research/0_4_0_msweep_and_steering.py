"""
Targets 2+3: Full m-sweep capture curve + observer optimisation on real molecules.

Target 2: For each matched pair, sweep m from 1 to n-1 with the canonical
          (soft-mode) observer C = [I_m | 0].

Target 3: For each matched pair at fixed m, optimise the observer on the
          Grassmannian to maximise the visible fraction. Compare with the
          static closure-adapted observer from nomogeo.

This is the "freezing point" experiment: can we steer the observer so the
problem enters a more exact sector?
"""

import sys
import time
import csv
import json
import numpy as np
from numpy.linalg import inv, eigh, norm, svd, slogdet
from scipy.linalg import null_space, sqrtm
from scipy.optimize import minimize
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

N_WORKERS = 20
DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
OUTPUT_DIR = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/papers/qm9_sweep_outputs")

ATOMIC_MASS = {1: 1.00782503223, 6: 12.0, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273}

# Add nomogeo to path
sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))


def sym(M):
    return 0.5 * (M + M.T)


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
    if h_raw.shape == (na, 3, na, 3):
        hc = sym(h_raw.reshape(3*na, 3*na))
    else:
        hc = sym(h_raw.transpose(0,2,1,3).reshape(3*na, 3*na))
    dm = np.repeat(masses, 3)
    ism = 1.0/np.sqrt(dm)
    hmw = sym((ism[:,None]*hc)*ism[None,:])
    w = np.sqrt(masses)
    ctr = np.average(pos, axis=0, weights=masses)
    ctrd = pos - ctr
    cols = []
    for ax in range(3):
        v = np.zeros((na,3)); v[:,ax] = w; cols.append(v.reshape(-1))
    for ax in np.eye(3):
        v = np.cross(ax[None,:], ctrd)*w[:,None]; cols.append(v.reshape(-1))
    raw = np.column_stack(cols)
    u, sv, _ = np.linalg.svd(raw, full_matrices=False)
    rk = int(np.sum(sv > 1e-10*max(1,float(np.max(sv)))))
    vb = null_space(u[:,:rk].T)
    hv = sym(vb.T @ hmw @ vb)
    ev = np.linalg.eigvalsh(hv)
    counts = {}
    for z in an: counts[int(z)] = counts.get(int(z),0)+1
    sm = {1:"H",6:"C",7:"N",8:"O",9:"F"}
    formula = "".join(sm[z]+(str(c) if c>1 else "") for z,c in sorted(counts.items()) if z in sm)
    return {'label': row.get('label','?'), 'h_vib': hv, 'n_vib': hv.shape[0],
            'spd': float(np.min(ev))>1e-6, 'formula': formula, 'n_atoms': na}


def compute_vis_frac(H, C, Hdot):
    """Compute visible fraction for given (H, C, Hdot)."""
    n = H.shape[0]
    try:
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        dHinv = -Hinv @ Hdot @ Hinv
        dPhi = sym(-Phi @ (C @ dHinv @ C.T) @ Phi)
        vis_rate = float(np.trace(inv(Phi) @ dPhi))
        amb_rate = float(np.trace(Hinv @ Hdot))
        if abs(amb_rate) < 1e-10:
            return None
        return vis_rate / amb_rate
    except:
        return None


def compute_v_min(H, C, Hdot):
    """Compute min eigenvalue of V = L^T Hdot L."""
    try:
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        V = L.T @ Hdot @ L
        return float(np.min(eigh(V)[0]))
    except:
        return -999


def process_molecule(args):
    """Full analysis for one molecule: m-sweep + observer optimisation."""
    label, h_vac, h_sol, n_vib, formula = args
    Hdot = h_sol - h_vac
    H = h_vac

    try:
        Hinv = inv(H)
        amb_rate = float(np.trace(Hinv @ Hdot))
    except:
        return None

    if abs(amb_rate) < 0.01:
        return None

    # --- m-sweep with canonical observer ---
    m_sweep = []
    for m in range(1, min(n_vib, 30)):  # cap at 30 for speed
        C = np.zeros((m, n_vib))
        for i in range(m): C[i,i] = 1.0
        vf = compute_vis_frac(H, C, Hdot)
        vm = compute_v_min(H, C, Hdot)
        if vf is not None:
            m_sweep.append({'m': m, 'vis_frac': vf, 'v_min': vm})

    # --- Observer optimisation at m=3 ---
    m_opt = min(3, n_vib - 1)
    if m_opt < 1 or n_vib < 4:
        return {'label': label, 'formula': formula, 'n_vib': n_vib,
                'amb_rate': amb_rate, 'm_sweep': m_sweep, 'opt': None}

    # Canonical observer (identity block)
    C_canon = np.zeros((m_opt, n_vib))
    for i in range(m_opt): C_canon[i,i] = 1.0
    vf_canon = compute_vis_frac(H, C_canon, Hdot)
    vm_canon = compute_v_min(H, C_canon, Hdot)

    # Optimise: find C that maximises |vis_frac| (or vis_frac if we want alignment)
    # Parameterise C as the first m rows of an orthogonal matrix Q
    # Q = exp(skew-symmetric A) applied to identity
    # But for a Grassmannian optimisation, use the Stiefel manifold parameterisation

    # Simpler: optimise over m x n matrices C with orthonormal rows
    # Use scipy with a projection step

    best_vf = vf_canon if vf_canon is not None else 0
    best_C = C_canon.copy()
    best_vm = vm_canon

    rng = np.random.default_rng(hash(label) % (2**31))

    # Try 30 random starting points
    for trial in range(30):
        # Random orthogonal frame
        Q = np.linalg.qr(rng.standard_normal((n_vib, n_vib)))[0]
        C_trial = Q[:m_opt, :]

        vf = compute_vis_frac(H, C_trial, Hdot)
        if vf is not None and vf > best_vf:
            best_vf = vf
            best_C = C_trial.copy()
            best_vm = compute_v_min(H, C_trial, Hdot)

    # Try the eigenvector-based observer (align with largest eigenvalues of H^{-1} Hdot)
    try:
        P = Hinv @ Hdot
        eigvals_P, eigvecs_P = eigh(sym(P))
        # Top m eigenvectors (largest eigenvalues)
        C_eig_top = eigvecs_P[:, -m_opt:].T
        vf = compute_vis_frac(H, C_eig_top, Hdot)
        if vf is not None and vf > best_vf:
            best_vf = vf
            best_C = C_eig_top.copy()
            best_vm = compute_v_min(H, C_eig_top, Hdot)

        # Most positive V eigenvalues (for V > 0 objective)
        # V = L^T Hdot L for each candidate C
    except:
        pass

    # Try the nomogeo closure-adapted observer
    try:
        from nomogeo import closure_adapted_observer
        eigvals_H, eigvecs_H = eigh(H)
        # Use H as a single-element family
        family = [Hdot]  # the perturbation as the "family member"
        # Actually closure_adapted_observer expects a precision family
        # Let's use H eigenbasis directly
        # The adapted observer for H itself
        result = closure_adapted_observer([H], m_opt)
        C_adapted = result.observer
        vf = compute_vis_frac(H, C_adapted, Hdot)
        vm = compute_v_min(H, C_adapted, Hdot)
        vf_adapted = vf
        vm_adapted = vm
    except:
        vf_adapted = None
        vm_adapted = None

    opt_result = {
        'vf_canon': vf_canon,
        'vm_canon': vm_canon,
        'vf_optimised': best_vf,
        'vm_optimised': best_vm,
        'vf_adapted': vf_adapted,
        'vm_adapted': vm_adapted,
        'improvement': best_vf - (vf_canon or 0),
    }

    return {'label': label, 'formula': formula, 'n_vib': n_vib,
            'amb_rate': amb_rate, 'm_sweep': m_sweep, 'opt': opt_result}


def main():
    t0 = time.time()
    print("=" * 70)
    print("Targets 2+3: m-Sweep + Observer Optimisation on Real Molecules")
    print("=" * 70)

    # Load matched pairs (500 per shard for broader coverage)
    per_shard = 500
    vac_mols = {}
    sol_mols = {}

    print(f"\n  Loading vacuum (500/shard)...")
    for shard in range(5):
        for row in iter_arrow_records("vacuum", shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: vac_mols[mol['label']] = mol
            except: pass
    print(f"  Loaded {len(vac_mols)} vacuum molecules")

    print(f"  Loading water (500/shard)...")
    for shard in range(5):
        for row in iter_arrow_records("water", shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']: sol_mols[mol['label']] = mol
            except: pass
    print(f"  Loaded {len(sol_mols)} water molecules")

    matched = set(vac_mols.keys()) & set(sol_mols.keys())
    pairs = []
    for label in matched:
        v = vac_mols[label]; s = sol_mols[label]
        if v['n_vib'] == s['n_vib']:
            pairs.append((label, v['h_vib'], s['h_vib'], v['n_vib'], v['formula']))
    print(f"  Matched pairs: {len(pairs)}")

    # Process all pairs
    print(f"\n  Processing ({N_WORKERS} workers)...")
    with Pool(N_WORKERS) as pool:
        results = [r for r in pool.imap_unordered(process_molecule, pairs, chunksize=5) if r is not None]
    print(f"  Processed {len(results)} molecules")

    # --- TARGET 2: m-sweep aggregate ---
    print(f"\n  === TARGET 2: Information Capture Curve ===")
    ratio_bins = defaultdict(list)
    for r in results:
        n = r['n_vib']
        for pt in r['m_sweep']:
            vf = pt['vis_frac']
            if abs(vf) < 100:
                ratio = pt['m'] / n
                bin_idx = min(int(ratio * 20), 19)
                ratio_bins[bin_idx].append(vf)

    print(f"  {'m/n':>6s}  {'n':>6s}  {'median':>8s}  {'p25':>8s}  {'p75':>8s}")
    curve_data = []
    for i in range(20):
        if ratio_bins[i]:
            vals = np.array(ratio_bins[i])
            mid = (i+0.5)/20
            med = np.median(vals)
            p25 = np.percentile(vals, 25)
            p75 = np.percentile(vals, 75)
            print(f"  {mid:6.3f}  {len(vals):6d}  {med:8.4f}  {p25:8.4f}  {p75:8.4f}")
            curve_data.append({'ratio': mid, 'n': len(vals), 'median': med, 'p25': p25, 'p75': p75})

    # Find 50% capture point
    for cd in curve_data:
        if cd['median'] >= 0.5:
            print(f"\n  50% capture at m/n ~ {cd['ratio']:.3f}")
            break

    # --- TARGET 3: Observer optimisation ---
    print(f"\n  === TARGET 3: Observer Optimisation ===")
    opt_results = [r for r in results if r['opt'] is not None]

    if opt_results:
        canon_vfs = [r['opt']['vf_canon'] for r in opt_results if r['opt']['vf_canon'] is not None]
        opt_vfs = [r['opt']['vf_optimised'] for r in opt_results]
        improvements = [r['opt']['improvement'] for r in opt_results]

        print(f"  Molecules with optimisation: {len(opt_results)}")
        print(f"  Canonical observer (m=3, soft modes):")
        print(f"    Median vis_frac: {np.median(canon_vfs):.4f}")
        print(f"    Mean vis_frac:   {np.mean(canon_vfs):.4f}")

        print(f"  Optimised observer (best of 30 random + eigenvector):")
        print(f"    Median vis_frac: {np.median(opt_vfs):.4f}")
        print(f"    Mean vis_frac:   {np.mean(opt_vfs):.4f}")

        print(f"  Improvement (optimised - canonical):")
        print(f"    Median: {np.median(improvements):.4f}")
        print(f"    Mean:   {np.mean(improvements):.4f}")
        print(f"    Max:    {np.max(improvements):.4f}")

        # V > 0 improvement
        canon_vpos = sum(1 for r in opt_results if r['opt']['vm_canon'] is not None and r['opt']['vm_canon'] > 1e-6)
        opt_vpos = sum(1 for r in opt_results if r['opt']['vm_optimised'] is not None and r['opt']['vm_optimised'] > 1e-6)
        print(f"\n  V > 0 (source law active):")
        print(f"    Canonical: {canon_vpos}/{len(opt_results)} ({100*canon_vpos/len(opt_results):.1f}%)")
        print(f"    Optimised: {opt_vpos}/{len(opt_results)} ({100*opt_vpos/len(opt_results):.1f}%)")

        # Adapted observer comparison
        adapted_vfs = [r['opt']['vf_adapted'] for r in opt_results if r['opt']['vf_adapted'] is not None]
        if adapted_vfs:
            adapted_vpos = sum(1 for r in opt_results if r['opt']['vm_adapted'] is not None and r['opt']['vm_adapted'] > 1e-6)
            print(f"  Closure-adapted observer:")
            print(f"    Median vis_frac: {np.median(adapted_vfs):.4f}")
            print(f"    V > 0: {adapted_vpos}/{len(adapted_vfs)}")

        # Show individual molecules with largest improvement
        opt_results_sorted = sorted(opt_results, key=lambda r: -r['opt']['improvement'])
        print(f"\n  Top 10 improvements:")
        print(f"  {'Label':>18s}  {'n':>4s}  {'canon':>8s}  {'optimised':>10s}  {'improve':>8s}  {'V>0':>4s}")
        for r in opt_results_sorted[:10]:
            o = r['opt']
            vp = 'Y' if o['vm_optimised'] and o['vm_optimised'] > 1e-6 else 'N'
            print(f"  {r['label']:>18s}  {r['n_vib']:4d}  {o['vf_canon'] or 0:8.4f}  "
                  f"{o['vf_optimised']:10.4f}  {o['improvement']:8.4f}  {vp:>4s}")

    t1 = time.time()

    # Write outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_DIR / "capture_curve_full.json", 'w') as f:
        json.dump(curve_data, f, indent=2)
    with open(OUTPUT_DIR / "observer_optimisation.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['label','formula','n_vib','amb_rate','vf_canon','vm_canon',
                     'vf_optimised','vm_optimised','vf_adapted','vm_adapted','improvement'])
        for r in results:
            if r['opt']:
                o = r['opt']
                w.writerow([r['label'], r['formula'], r['n_vib'], r['amb_rate'],
                           o['vf_canon'], o['vm_canon'], o['vf_optimised'], o['vm_optimised'],
                           o['vf_adapted'], o['vm_adapted'], o['improvement']])

    print(f"\n  Elapsed: {t1-t0:.1f}s")


if __name__ == "__main__":
    main()
