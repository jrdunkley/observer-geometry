"""
Full HessianQM9 solvent sweep — 20-worker parallel execution.

For every molecule appearing in vacuum AND at least one solvent,
compute the conservation law information budget.

Outputs:
  - CSV with per-molecule diagnostics
  - JSON summary statistics
  - Markdown report

Usage:
  python 0_4_0_qm9_full_sweep.py
"""

import json
import csv
import math
import time
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det, slogdet
from scipy.linalg import sqrtm, null_space

# ---- Configuration ----

N_WORKERS = 20
M_VISIBLE = 3  # observe the 3 softest modes
PER_SHARD = None  # None = load all records in shard

DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
OUTPUT_DIR = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/papers/qm9_sweep_outputs")

SOLVENTS = ["thf", "toluene", "water"]

ATOMIC_MASS = {1: 1.00782503223, 6: 12.0, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273}


# ---- Data loading ----

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
                if limit is not None and seen >= limit:
                    return


def sym(M):
    return 0.5 * (M + M.T)


def hessian_blocks_to_cartesian(blocks, n_atoms):
    array = np.asarray(blocks, dtype=float)
    if array.shape == (n_atoms, 3, n_atoms, 3):
        return sym(array.reshape(3 * n_atoms, 3 * n_atoms))
    if array.shape == (n_atoms, n_atoms, 3, 3):
        return sym(array.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms))
    raise ValueError(f"unexpected Hessian shape {array.shape}")


def rigid_motion_basis(positions, masses, tol=1e-10):
    weights = np.sqrt(masses)
    center = np.average(positions, axis=0, weights=masses)
    centered = positions - center
    columns = []
    for axis in range(3):
        vec = np.zeros((positions.shape[0], 3))
        vec[:, axis] = weights
        columns.append(vec.reshape(-1))
    for axis in np.eye(3):
        vec = np.cross(axis[None, :], centered) * weights[:, None]
        columns.append(vec.reshape(-1))
    raw = np.column_stack(columns)
    u, sv, _ = np.linalg.svd(raw, full_matrices=False)
    cutoff = max(tol, tol * max(1.0, float(np.max(sv))))
    rank = int(np.sum(sv > cutoff))
    return u[:, :rank]


def vibrational_hessian(row):
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    n_atoms = int(atomic_numbers.shape[0])
    masses = np.array([ATOMIC_MASS[int(z)] for z in atomic_numbers])

    h_cart = hessian_blocks_to_cartesian(row["hessian"], n_atoms)
    diag_m = np.repeat(masses, 3)
    inv_sqrt_m = 1.0 / np.sqrt(diag_m)
    h_mw = sym((inv_sqrt_m[:, None] * h_cart) * inv_sqrt_m[None, :])

    rigid = rigid_motion_basis(positions, masses)
    vib_basis = null_space(rigid.T)
    h_vib = sym(vib_basis.T @ h_mw @ vib_basis)

    eigenvalues = np.linalg.eigvalsh(h_vib)
    min_eig = float(np.min(eigenvalues)) if eigenvalues.size else -1.0

    # Formula string
    counts = {}
    for z in atomic_numbers:
        counts[int(z)] = counts.get(int(z), 0) + 1
    order = [6, 1, 7, 8, 9]
    sym_map = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
    formula = ""
    for z in order:
        c = counts.pop(z, 0)
        if c:
            formula += sym_map[z] + (str(c) if c > 1 else "")

    return {
        'label': row.get('label', 'unknown'),
        'formula': formula,
        'h_vib': h_vib,
        'n_vib': h_vib.shape[0],
        'n_atoms': n_atoms,
        'min_eig': min_eig,
        'spd': min_eig > 1e-6,
    }


# ---- Source law computation ----

def compute_budget(h_vac, h_sol, m):
    """Compute the information budget for a vacuum-solvent pair."""
    n = h_vac.shape[0]
    if n < m + 1:
        return None

    C = np.zeros((m, n))
    for i in range(m):
        C[i, i] = 1.0

    H = h_vac
    Hdot = h_sol - h_vac

    try:
        Hinv = inv(H)
        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C)
        Z = Vt[m:].T
        R = Z.T @ H @ Z
        Rinv = inv(R)
    except np.linalg.LinAlgError:
        return None

    V = L.T @ Hdot @ L
    U_h = Z.T @ Hdot @ Z
    B = L.T @ Hdot @ Z

    # Full conservation with connection
    dHinv = -Hinv @ Hdot @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    vis_rate = np.trace(inv(Phi) @ dPhi)
    hid_rate = np.trace(Rinv @ U_h)
    amb_rate = np.trace(Hinv @ Hdot)
    cons_err = abs(vis_rate + hid_rate - amb_rate)

    v_eigvals = sorted(eigh(V)[0])
    phi_eigvals = sorted(eigh(Phi)[0])

    # Hidden defect
    Qhat_trace = np.trace(B @ Rinv @ B.T)

    # A_cpl if V > 0
    a_cpl_min = None
    a_cpl_max = None
    a_cpl_trace = None
    if v_eigvals[0] > 1e-6:
        try:
            Vsqrt = np.real(sqrtm(V))
            Vsqrt_inv = inv(Vsqrt)
            theta = -Rinv @ B.T
            W = theta.T @ B.T + B @ theta
            A_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)
            a_eigvals = eigh(A_cpl)[0]
            a_cpl_min = float(np.min(a_eigvals))
            a_cpl_max = float(np.max(a_eigvals))
            a_cpl_trace = float(np.trace(A_cpl))
        except:
            pass

    return {
        'n_vib': n,
        'phi_min': float(phi_eigvals[0]),
        'phi_max': float(phi_eigvals[-1]),
        'v_min': float(v_eigvals[0]),
        'v_max': float(v_eigvals[-1]),
        'vis_rate': float(vis_rate),
        'hid_rate': float(hid_rate),
        'amb_rate': float(amb_rate),
        'cons_err': float(cons_err),
        'vis_frac': float(vis_rate / amb_rate) if abs(amb_rate) > 1e-15 else float('nan'),
        'qhat_trace': float(Qhat_trace),
        'v_positive': v_eigvals[0] > 1e-6,
        'a_cpl_min': a_cpl_min,
        'a_cpl_max': a_cpl_max,
        'a_cpl_trace': a_cpl_trace,
    }


# ---- Worker functions ----

def load_shard(args):
    """Load one shard of one split. Returns {label: vibrational_hessian_dict}."""
    split, shard = args
    results = {}
    for row in iter_arrow_records(split, shard, PER_SHARD):
        try:
            mol = vibrational_hessian(row)
            if mol['spd']:
                results[mol['label']] = mol
        except:
            pass
    return split, shard, results


def process_pair(args):
    """Process one molecule across vacuum and solvent."""
    label, h_vac, h_sol, solvent, formula, n_atoms = args
    budget = compute_budget(h_vac, h_sol, M_VISIBLE)
    if budget is None:
        return None
    budget['label'] = label
    budget['formula'] = formula
    budget['solvent'] = solvent
    budget['n_atoms'] = n_atoms
    return budget


def main():
    t0 = time.time()
    print("=" * 70)
    print(f"HessianQM9 Full Solvent Sweep ({N_WORKERS} workers)")
    print("=" * 70)

    if not DATA_ROOT.exists():
        print(f"  ERROR: Data root not found: {DATA_ROOT}")
        return 1

    # Phase 1: Load all shards in parallel
    print(f"\n  Phase 1: Loading Arrow data...")
    load_tasks = []
    for split in ["vacuum"] + SOLVENTS:
        for shard in range(5):
            if (DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow").exists():
                load_tasks.append((split, shard))

    print(f"  {len(load_tasks)} shard files to load across {1 + len(SOLVENTS)} splits")

    all_mols = {}  # {split: {label: mol_dict}}
    for split in ["vacuum"] + SOLVENTS:
        all_mols[split] = {}

    with Pool(min(N_WORKERS, len(load_tasks))) as pool:
        for split, shard, results in pool.imap_unordered(load_shard, load_tasks):
            all_mols[split].update(results)
            print(f"    Loaded {split}/shard-{shard}: {len(results)} SPD molecules")

    vac_labels = set(all_mols["vacuum"].keys())
    print(f"\n  Vacuum: {len(vac_labels)} SPD molecules")
    for solvent in SOLVENTS:
        sol_labels = set(all_mols[solvent].keys())
        overlap = vac_labels & sol_labels
        print(f"  {solvent}: {len(all_mols[solvent])} SPD, {len(overlap)} matched with vacuum")

    # Phase 2: Compute budgets in parallel
    print(f"\n  Phase 2: Computing information budgets...")
    pair_tasks = []
    for solvent in SOLVENTS:
        matched = vac_labels & set(all_mols[solvent].keys())
        for label in matched:
            vac = all_mols["vacuum"][label]
            sol = all_mols[solvent][label]
            if vac['n_vib'] == sol['n_vib']:
                pair_tasks.append((
                    label, vac['h_vib'], sol['h_vib'],
                    solvent, vac['formula'], vac['n_atoms'],
                ))

    print(f"  {len(pair_tasks)} valid pairs to process")

    results = []
    with Pool(N_WORKERS) as pool:
        for budget in pool.imap_unordered(process_pair, pair_tasks, chunksize=10):
            if budget is not None:
                results.append(budget)

    t1 = time.time()
    print(f"\n  Computed {len(results)} budgets in {t1 - t0:.1f}s")

    if not results:
        print("  No results. Check data.")
        return 1

    # Phase 3: Write outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # CSV
    csv_path = OUTPUT_DIR / "solvent_information_budgets.csv"
    fieldnames = [
        'label', 'formula', 'solvent', 'n_atoms', 'n_vib',
        'vis_rate', 'hid_rate', 'amb_rate', 'vis_frac', 'cons_err',
        'phi_min', 'phi_max', 'v_min', 'v_max', 'v_positive',
        'qhat_trace', 'a_cpl_min', 'a_cpl_max', 'a_cpl_trace',
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in sorted(results, key=lambda x: (x['solvent'], x['label'])):
            writer.writerow({k: r.get(k) for k in fieldnames})

    # Summary statistics
    cons_errs = [r['cons_err'] for r in results]
    vis_fracs = [r['vis_frac'] for r in results if not math.isnan(r['vis_frac'])]
    v_pos = sum(1 for r in results if r['v_positive'])

    per_solvent = {}
    for solvent in SOLVENTS:
        sol_results = [r for r in results if r['solvent'] == solvent]
        if sol_results:
            sol_vis_fracs = [r['vis_frac'] for r in sol_results if not math.isnan(r['vis_frac'])]
            per_solvent[solvent] = {
                'count': len(sol_results),
                'v_positive': sum(1 for r in sol_results if r['v_positive']),
                'mean_vis_frac': float(np.mean(sol_vis_fracs)) if sol_vis_fracs else None,
                'median_vis_frac': float(np.median(sol_vis_fracs)) if sol_vis_fracs else None,
                'max_cons_err': float(max(r['cons_err'] for r in sol_results)),
            }

    summary = {
        'total_pairs': len(results),
        'max_conservation_error': float(max(cons_errs)),
        'v_positive_count': v_pos,
        'v_positive_fraction': v_pos / len(results) if results else 0,
        'mean_vis_fraction': float(np.mean(vis_fracs)) if vis_fracs else None,
        'median_vis_fraction': float(np.median(vis_fracs)) if vis_fracs else None,
        'per_solvent': per_solvent,
        'elapsed_seconds': t1 - t0,
        'm_visible': M_VISIBLE,
    }

    json_path = OUTPUT_DIR / "solvent_sweep_summary.json"
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)

    # Markdown report
    md_path = OUTPUT_DIR / "solvent_sweep_report.md"
    with open(md_path, 'w') as f:
        f.write("# HessianQM9 Full Solvent Sweep\n\n")
        f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Workers:** {N_WORKERS}\n")
        f.write(f"**Visible modes (m):** {M_VISIBLE}\n")
        f.write(f"**Elapsed:** {t1-t0:.1f}s\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"| Metric | Value |\n|--------|-------|\n")
        f.write(f"| Total pairs | {len(results)} |\n")
        f.write(f"| Max conservation error | {max(cons_errs):.2e} |\n")
        f.write(f"| V > 0 (source law active) | {v_pos}/{len(results)} ({100*v_pos/len(results):.1f}%) |\n")
        f.write(f"| Mean visible fraction | {np.mean(vis_fracs):.4f} |\n")
        f.write(f"| Median visible fraction | {np.median(vis_fracs):.4f} |\n\n")
        f.write(f"## Per-solvent\n\n")
        f.write(f"| Solvent | Pairs | V>0 | Mean vis frac | Max cons err |\n")
        f.write(f"|---------|-------|-----|---------------|-------------|\n")
        for sol, stats in per_solvent.items():
            f.write(f"| {sol} | {stats['count']} | {stats['v_positive']} | "
                    f"{stats['mean_vis_frac']:.4f} | {stats['max_cons_err']:.2e} |\n")

    print(f"\n  Outputs written to {OUTPUT_DIR}/")
    print(f"    {csv_path.name}")
    print(f"    {json_path.name}")
    print(f"    {md_path.name}")

    # Print summary
    print(f"\n  === RESULTS ===")
    print(f"  Total pairs: {len(results)}")
    print(f"  Max conservation error: {max(cons_errs):.2e}")
    print(f"  V > 0: {v_pos}/{len(results)} ({100*v_pos/len(results):.1f}%)")
    print(f"  Mean visible fraction: {np.mean(vis_fracs):.4f}")
    print(f"  Median visible fraction: {np.median(vis_fracs):.4f}")
    for sol, stats in per_solvent.items():
        print(f"  {sol}: {stats['count']} pairs, mean vis_frac = {stats['mean_vis_frac']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
