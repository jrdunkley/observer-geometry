"""
Source law on HessianQM9 full Hessian data.

For each molecule appearing in multiple solvent environments, the solvent
perturbation defines a path through SPD(n_vib). The source law tracks how
visible precision changes along this path.

Strategy:
1. Load molecules from vacuum and one solvent (e.g., THF)
2. Match by label (same molecule in two environments)
3. For each pair: H_vacuum, H_thf are both vibrational Hessians (SPD)
4. Path: H(t) = (1-t)*H_vac + t*H_thf, Hdot = H_thf - H_vac
5. Observer C selects a subset of vibrational modes
6. Compute source law, conservation, curvature
"""

import sys
import math
import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det, slogdet
from scipy.linalg import sqrtm, null_space
from pathlib import Path

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok

def sym(M):
    return 0.5 * (M + M.T)

# --- Data loading (adapted from molecular_vibrational_atlas.py) ---

BASE = Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src/utils/operationalise/attachment_layer")
DATA_ROOT = BASE / "incoming_data" / "hessian_qm9"

# Try multiple possible directory structures
POSSIBLE_ROOTS = [
    Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict"),
    DATA_ROOT / "hessian_qm9_DatasetDict",
    DATA_ROOT / "DatasetQM9_DatasetDict",
    DATA_ROOT,
]

ATOMIC_MASS = {1: 1.00782503223, 6: 12.0, 7: 14.00307400443, 8: 15.99491461957, 9: 18.99840316273}


def find_data_root():
    for root in POSSIBLE_ROOTS:
        if root.exists():
            # Check for split directories
            for split in ["vacuum", "thf", "toluene", "water"]:
                if (root / split).exists():
                    return root
    return None


def iter_arrow_records(data_root, split, shard, limit):
    try:
        import pyarrow.ipc as ipc
    except ImportError:
        print("  pyarrow not available, skipping Arrow data")
        return

    path = data_root / split / f"data-{shard:05d}-of-00005.arrow"
    if not path.exists():
        return
    seen = 0
    with ipc.open_stream(path) as reader:
        for batch in reader:
            for row in batch.to_pylist():
                yield row
                seen += 1
                if seen >= limit:
                    return


def hessian_blocks_to_cartesian(blocks, n_atoms):
    array = np.asarray(blocks, dtype=float)
    if array.shape == (n_atoms, 3, n_atoms, 3):
        return sym(array.reshape(3 * n_atoms, 3 * n_atoms))
    if array.shape == (n_atoms, n_atoms, 3, 3):
        return sym(array.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms))
    raise ValueError(f"unexpected Hessian shape {array.shape}")


def mass_weight_hessian(h_cart, masses):
    diag_m = np.repeat(masses, 3)
    inv_sqrt_m = 1.0 / np.sqrt(diag_m)
    return sym((inv_sqrt_m[:, None] * h_cart) * inv_sqrt_m[None, :]), inv_sqrt_m


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
    """Extract vibrational SPD Hessian from an Arrow record."""
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    n_atoms = int(atomic_numbers.shape[0])
    masses = np.array([ATOMIC_MASS[int(z)] for z in atomic_numbers])

    h_cart = hessian_blocks_to_cartesian(row["hessian"], n_atoms)
    h_mw, _ = mass_weight_hessian(h_cart, masses)

    rigid = rigid_motion_basis(positions, masses)
    vib_basis = null_space(rigid.T)
    h_vib = sym(vib_basis.T @ h_mw @ vib_basis)

    eigenvalues = np.linalg.eigvalsh(h_vib)
    min_eig = float(np.min(eigenvalues))

    formula = ''.join(
        f"{chr(z)}{c}" if c > 1 else chr(z)
        for z, c in sorted(
            {int(z): int(np.sum(atomic_numbers == z)) for z in np.unique(atomic_numbers)}.items()
        )
    )

    return {
        'label': row.get('label', row.get('mol_id', 'unknown')),
        'formula': formula,
        'h_vib': h_vib,
        'n_vib': h_vib.shape[0],
        'min_eig': min_eig,
        'max_eig': float(np.max(eigenvalues)),
        'spd': min_eig > 1e-6,
        'positions': positions,
        'masses': masses,
    }


# --- Source law computation ---

def source_law_solvent_pair(mol_vac, mol_sol, m=3):
    """
    Compute source law between vacuum and solvent Hessians.
    """
    H_vac = mol_vac['h_vib']
    H_sol = mol_sol['h_vib']
    n = H_vac.shape[0]

    if n != H_sol.shape[0]:
        return None
    if not mol_vac['spd'] or not mol_sol['spd']:
        return None

    # Observer: lowest m vibrational modes
    m = min(m, n // 2)
    if m < 1:
        return None

    C = np.zeros((m, n))
    for i in range(m):
        C[i, i] = 1.0

    H = H_vac
    Hdot = H_sol - H_vac

    Hinv = inv(H)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C)
    Z = Vt[m:].T
    R = Z.T @ H @ Z
    Rinv = inv(R)

    V = L.T @ Hdot @ L
    B = L.T @ Hdot @ Z
    U_h = Z.T @ Hdot @ Z

    phi_eigvals = sorted(eigh(Phi)[0])
    v_eigvals = sorted(eigh(V)[0])

    # Conservation
    vis_rate = np.trace(inv(Phi) @ V)  # + 2*Tr(alpha), but alpha contribution is small for diagonal C
    hid_rate = np.trace(Rinv @ U_h)
    amb_rate = np.trace(Hinv @ Hdot)

    # Full conservation with connection
    dHinv = -Hinv @ Hdot @ Hinv
    dPhi = -Phi @ (C @ dHinv @ C.T) @ Phi
    alpha = inv(Phi) @ (L.T @ H @ (dHinv @ C.T @ Phi + Hinv @ C.T @ dPhi))
    vis_rate_full = np.trace(inv(Phi) @ dPhi)  # = Tr(Phi^{-1} V) + 2 Tr(alpha)
    cons_err = abs(vis_rate_full + hid_rate - amb_rate)

    # A_cpl
    a_cpl_eigvals = None
    if np.min(v_eigvals) > 1e-6:
        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        theta = -Rinv @ B.T
        W = theta.T @ B.T + B @ theta  # Hddot = 0 for linear path
        A_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)
        a_cpl_eigvals = sorted(eigh(A_cpl)[0])

    return {
        'label': mol_vac['label'],
        'formula': mol_vac['formula'],
        'n_vib': n,
        'm': m,
        'phi_min': phi_eigvals[0],
        'v_min': v_eigvals[0],
        'v_max': v_eigvals[-1],
        'vis_rate': vis_rate_full,
        'hid_rate': hid_rate,
        'amb_rate': amb_rate,
        'cons_err': cons_err,
        'a_cpl_eigvals': a_cpl_eigvals,
        'vis_frac': vis_rate_full / amb_rate if abs(amb_rate) > 1e-15 else float('nan'),
    }


def main():
    print("=" * 70)
    print("Source Law on HessianQM9 (Full Hessian, Solvent Perturbation)")
    print("=" * 70)

    data_root = find_data_root()
    if data_root is None:
        # List what's actually in the data directory
        print(f"\n  Data root not found. Checking directories...")
        for p in POSSIBLE_ROOTS:
            print(f"    {p}: exists={p.exists()}")
            if p.exists():
                for item in sorted(p.iterdir())[:10]:
                    print(f"      {item.name}")
        print("\n  Waiting for extraction to complete.")
        return

    print(f"\n  Data root: {data_root}")

    # Check available splits
    available = []
    for split in ["vacuum", "thf", "toluene", "water"]:
        split_dir = data_root / split
        if split_dir.exists():
            n_files = len(list(split_dir.glob("*.arrow")))
            available.append((split, n_files))
            print(f"  Split '{split}': {n_files} Arrow files")

    if len(available) < 2:
        print("  Need at least 2 splits (vacuum + solvent). Waiting for extraction.")
        return

    # Load molecules from all shards (first 200 per shard to limit memory)
    per_shard = 200
    solvent = "thf"

    print(f"\n  Loading vacuum (all shards, {per_shard} per shard)...")
    vac_mols = {}
    for shard in range(5):
        for row in iter_arrow_records(data_root, "vacuum", shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']:
                    vac_mols[mol['label']] = mol
            except:
                pass
    print(f"  Loaded {len(vac_mols)} SPD vacuum molecules")

    print(f"  Loading {solvent} (all shards, {per_shard} per shard)...")
    sol_mols = {}
    for shard in range(5):
        for row in iter_arrow_records(data_root, solvent, shard, per_shard):
            try:
                mol = vibrational_hessian(row)
                if mol['spd']:
                    sol_mols[mol['label']] = mol
            except:
                pass
    print(f"  Loaded {len(sol_mols)} SPD {solvent} molecules")

    # Match molecules
    matched = set(vac_mols.keys()) & set(sol_mols.keys())
    print(f"  Matched molecules: {len(matched)}")

    if not matched:
        print("  No matched molecules found.")
        return

    # Compute source law for matched pairs
    results = []
    for label in sorted(matched):
        result = source_law_solvent_pair(vac_mols[label], sol_mols[label], m=3)
        if result is not None:
            results.append(result)

    if not results:
        print("  No valid source law computations.")
        return

    # Conservation check
    cons_errs = [r['cons_err'] for r in results]
    max_cons = max(cons_errs)
    report(f"Conservation across {len(results)} solvent pairs", max_cons)

    # Print results
    print(f"\n  {'Label':>10s}  {'n':>4s}  {'m':>3s}  {'vis_rate':>10s}  {'hid_rate':>10s}  "
          f"{'vis_frac':>8s}  {'V>0':>4s}")
    for r in results[:30]:
        v_pos = 'Y' if r['v_min'] > 1e-6 else 'N'
        print(f"  {str(r['label'])[:10]:>10s}  {r['n_vib']:4d}  {r['m']:3d}  "
              f"{r['vis_rate']:10.4f}  {r['hid_rate']:10.4f}  "
              f"{r['vis_frac']:8.3f}  {v_pos:>4s}")

    # Statistics
    vis_rates = [r['vis_rate'] for r in results]
    hid_rates = [r['hid_rate'] for r in results]
    vis_fracs = [r['vis_frac'] for r in results if not np.isnan(r['vis_frac'])]
    v_positive = sum(1 for r in results if r['v_min'] > 1e-6)

    print(f"\n  Summary ({solvent} vs vacuum):")
    print(f"    Molecules analysed: {len(results)}")
    print(f"    V > 0 (source law active): {v_positive}/{len(results)}")
    print(f"    Mean visible rate: {np.mean(vis_rates):.6f}")
    print(f"    Mean hidden rate: {np.mean(hid_rates):.6f}")
    print(f"    Mean visible fraction: {np.mean(vis_fracs):.4f}")
    print(f"    Max conservation error: {max_cons:.2e}")

    # A_cpl analysis for V > 0 cases
    acpl_results = [r for r in results if r['a_cpl_eigvals'] is not None]
    if acpl_results:
        print(f"\n  A_cpl analysis ({len(acpl_results)} molecules with V > 0):")
        for r in acpl_results[:10]:
            eigstr = ', '.join(f'{e:.4f}' for e in r['a_cpl_eigvals'])
            print(f"    {str(r['label'])[:10]}: A_cpl eigvals = [{eigstr}]")

    print(f"\n  RESULTS: {passes} passed, {fails} failed")


if __name__ == "__main__":
    main()
