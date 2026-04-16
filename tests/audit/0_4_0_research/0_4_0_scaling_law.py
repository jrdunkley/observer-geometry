"""
Target 2: Information capture curve — how does visible fraction scale with m?

For each molecule-solvent pair, sweep m from 1 to n_vib-1 and compute
the visible fraction. This gives the "information capture curve" that
characterises how much solvent response each additional mode captures.

Uses the existing sweep CSV data.
"""

import csv
import sys
import numpy as np
from numpy.linalg import inv, eigh, svd, slogdet
from scipy.linalg import null_space
from pathlib import Path
from multiprocessing import Pool

N_WORKERS = 20
DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
OUTPUT_DIR = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/papers/qm9_sweep_outputs")

ATOMIC_MASS = {1: 1.00782503223, 6: 12.0, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273}


def sym(M):
    return 0.5 * (M + M.T)


def iter_arrow_records(split, shard, limit=None):
    import pyarrow.ipc as ipc
    path = DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow"
    if not path.exists():
        return
    seen = 0
    with ipc.open_stream(path) as reader:
        for batch in reader:
            for row in batch.to_pylist():
                yield row
                seen += 1
                if limit and seen >= limit:
                    return


def vibrational_hessian(row):
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    n_atoms = int(atomic_numbers.shape[0])
    masses = np.array([ATOMIC_MASS[int(z)] for z in atomic_numbers])
    h_raw = np.asarray(row["hessian"], dtype=float)
    if h_raw.shape == (n_atoms, 3, n_atoms, 3):
        h_cart = sym(h_raw.reshape(3*n_atoms, 3*n_atoms))
    elif h_raw.shape == (n_atoms, n_atoms, 3, 3):
        h_cart = sym(h_raw.transpose(0,2,1,3).reshape(3*n_atoms, 3*n_atoms))
    else:
        raise ValueError(f"bad shape {h_raw.shape}")
    diag_m = np.repeat(masses, 3)
    ism = 1.0 / np.sqrt(diag_m)
    h_mw = sym((ism[:,None] * h_cart) * ism[None,:])
    weights = np.sqrt(masses)
    center = np.average(positions, axis=0, weights=masses)
    centered = positions - center
    cols = []
    for axis in range(3):
        vec = np.zeros((n_atoms, 3)); vec[:, axis] = weights; cols.append(vec.reshape(-1))
    for axis in np.eye(3):
        vec = np.cross(axis[None,:], centered) * weights[:,None]; cols.append(vec.reshape(-1))
    raw = np.column_stack(cols)
    u, sv, _ = np.linalg.svd(raw, full_matrices=False)
    rank = int(np.sum(sv > 1e-10 * max(1, float(np.max(sv)))))
    vib_basis = null_space(u[:,:rank].T)
    h_vib = sym(vib_basis.T @ h_mw @ vib_basis)
    eigvals = np.linalg.eigvalsh(h_vib)
    spd = float(np.min(eigvals)) > 1e-6
    return {'label': row.get('label', '?'), 'h_vib': h_vib, 'n_vib': h_vib.shape[0], 'spd': spd,
            'n_atoms': n_atoms}


def compute_capture_curve(args):
    """For one molecule pair, compute vis_frac at each m from 1 to n-1."""
    label, h_vac, h_sol, solvent, n_vib = args

    Hdot = h_sol - h_vac
    H = h_vac
    try:
        Hinv = inv(H)
    except:
        return None

    amb_rate = float(np.trace(inv(H) @ Hdot))
    if abs(amb_rate) < 0.01:
        return None  # degenerate denominator

    results = []
    for m in range(1, n_vib):
        C = np.zeros((m, n_vib))
        for i in range(m):
            C[i, i] = 1.0

        try:
            Phi = inv(C @ Hinv @ C.T)
            dHinv = -Hinv @ Hdot @ Hinv
            dPhi = sym(-Phi @ (C @ dHinv @ C.T) @ Phi)
            vis_rate = float(np.trace(inv(Phi) @ dPhi))
        except:
            vis_rate = float('nan')

        results.append({
            'm': m,
            'vis_rate': vis_rate,
            'vis_frac': vis_rate / amb_rate,
        })

    return {
        'label': label,
        'solvent': solvent,
        'n_vib': n_vib,
        'amb_rate': amb_rate,
        'curve': results,
    }


def load_pairs(solvent, max_per_shard=100):
    """Load matched vacuum-solvent pairs."""
    vac_mols = {}
    sol_mols = {}

    for shard in range(5):
        for row in iter_arrow_records("vacuum", shard, max_per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']:
                    vac_mols[mol['label']] = mol
            except:
                pass

    for shard in range(5):
        for row in iter_arrow_records(solvent, shard, max_per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']:
                    sol_mols[mol['label']] = mol
            except:
                pass

    matched = set(vac_mols.keys()) & set(sol_mols.keys())
    pairs = []
    for label in matched:
        v = vac_mols[label]
        s = sol_mols[label]
        if v['n_vib'] == s['n_vib']:
            pairs.append((label, v['h_vib'], s['h_vib'], solvent, v['n_vib']))

    return pairs


def main():
    import time
    t0 = time.time()
    print("=" * 70)
    print("Target 2: Information Capture Curve (m-sweep)")
    print("=" * 70)

    # Load a sample of pairs (100 per shard = ~500 per split, enough for statistics)
    pairs = load_pairs("water", max_per_shard=100)
    print(f"  Loaded {len(pairs)} matched water-vacuum pairs")

    if not pairs:
        print("  No pairs found")
        return

    # Compute capture curves in parallel
    print(f"  Computing capture curves ({N_WORKERS} workers)...")
    with Pool(N_WORKERS) as pool:
        results = [r for r in pool.imap_unordered(compute_capture_curve, pairs, chunksize=5) if r is not None]

    print(f"  Computed {len(results)} capture curves")

    # Aggregate: for each m/n ratio, what is the median visible fraction?
    # Group by n_vib to avoid mixing different-sized molecules
    from collections import defaultdict

    # Normalised capture curve: m/n_vib vs vis_frac
    ratio_bins = np.linspace(0, 1, 21)  # 5% bins
    bin_data = defaultdict(list)

    for r in results:
        n = r['n_vib']
        for point in r['curve']:
            m = point['m']
            vf = point['vis_frac']
            if not np.isnan(vf) and abs(vf) < 100:
                ratio = m / n
                bin_idx = int(ratio * 20)
                bin_idx = min(bin_idx, 19)
                bin_data[bin_idx].append(vf)

    print(f"\n  === Information Capture Curve ===")
    print(f"  {'m/n':>6s}  {'n_pts':>6s}  {'median_vf':>10s}  {'p25':>8s}  {'p75':>8s}  {'mean':>8s}")
    medians = []
    ratios = []
    for i in range(20):
        if bin_data[i]:
            vals = np.array(bin_data[i])
            med = np.median(vals)
            p25 = np.percentile(vals, 25)
            p75 = np.percentile(vals, 75)
            mn = np.mean(vals)
            midpoint = (i + 0.5) / 20
            print(f"  {midpoint:6.3f}  {len(vals):6d}  {med:10.4f}  {p25:8.4f}  {p75:8.4f}  {mn:8.4f}")
            medians.append(med)
            ratios.append(midpoint)

    # Key statistics
    if medians:
        # Find the m/n where median vis_frac crosses 0.5
        half_idx = None
        for i in range(1, len(medians)):
            if medians[i-1] < 0.5 and medians[i] >= 0.5:
                half_idx = i
                break
            elif medians[i-1] >= 0.5 and i == 1:
                half_idx = 0
                break

        if half_idx is not None:
            print(f"\n  50% capture at m/n ~ {ratios[half_idx]:.2f}")

        # Marginal return: how much does each additional 5% of modes add?
        print(f"\n  Marginal information per 5% of modes:")
        for i in range(1, len(medians)):
            marginal = medians[i] - medians[i-1]
            print(f"    {ratios[i-1]:.2f} -> {ratios[i]:.2f}: +{marginal:.4f}")

    # Per-dimension analysis
    print(f"\n  Per-dimension capture at m=3:")
    dim_data = defaultdict(list)
    for r in results:
        n = r['n_vib']
        for point in r['curve']:
            if point['m'] == 3 and abs(point['vis_frac']) < 100 and not np.isnan(point['vis_frac']):
                dim_data[n].append(point['vis_frac'])

    for n_vib in sorted(dim_data.keys()):
        vals = dim_data[n_vib]
        if len(vals) >= 3:
            print(f"    n_vib={n_vib:3d}: median vis_frac={np.median(vals):.4f} (n={len(vals)})")

    t1 = time.time()

    # Write output
    out_path = OUTPUT_DIR / "capture_curve.csv"
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['label', 'solvent', 'n_vib', 'amb_rate', 'm', 'vis_rate', 'vis_frac'])
        for r in results:
            for point in r['curve']:
                writer.writerow([r['label'], r['solvent'], r['n_vib'], r['amb_rate'],
                                point['m'], point['vis_rate'], point['vis_frac']])

    print(f"\n  Output: {out_path}")
    print(f"  Elapsed: {t1-t0:.1f}s")


if __name__ == "__main__":
    main()
