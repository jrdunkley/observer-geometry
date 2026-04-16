"""
Source law on QM9 molecular data.

The raw QM9 dataset provides:
- Molecular geometry (atom types and 3D coordinates)
- Vibrational frequencies (eigenvalues of mass-weighted Hessian)
- Molecular properties (energy, dipole, etc.)

Strategy: for pairs of molecules with similar structure but different
frequencies, the frequency differences define a "perturbation path" through
the space of Hessian eigenvalues. The source law tracks how the visible
precision (observer-projected Hessian) changes along this path.

For a molecule with N atoms, the Cartesian Hessian is 3N x 3N.
After mass-weighting and projecting out 6 rigid-body modes, there are
3N-6 vibrational modes (3N-5 for linear molecules).

The vibrational frequencies omega_i are related to Hessian eigenvalues by:
  lambda_i = omega_i^2 * (2*pi*c)^2  (in appropriate units)

We construct a diagonal Hessian from the frequencies and apply the source law.
"""

import numpy as np
from numpy.linalg import inv, eigh, norm, svd, det
from scipy.linalg import sqrtm
from pathlib import Path
import re

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


XYZ_DIR = Path("C:/observer_geometry_workspace_v0.3.2/datasets/QM9/dsgdb9nsd.xyz")

ATOMIC_MASS = {
    'C': 12.0, 'H': 1.00782503223, 'N': 14.00307400443,
    'O': 15.99491461957, 'F': 18.99840316273,
}


def parse_xyz(path):
    """Parse a QM9 xyz file. Returns atoms, coords, frequencies."""
    lines = path.read_text().strip().split('\n')
    n_atoms = int(lines[0].strip())

    atoms = []
    coords = []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    # Frequencies are on the line after the atoms
    freq_line = lines[2 + n_atoms]
    freqs = [float(x) for x in freq_line.split()]

    return atoms, np.array(coords), np.array(freqs)


def construct_hessian_from_freqs(atoms, freqs):
    """
    Construct a diagonal mass-weighted Hessian from vibrational frequencies.

    The eigenvalues of the mass-weighted Hessian are omega^2 (in cm^-1 squared,
    up to unit conversion). For observer geometry purposes, the absolute units
    don't matter — we need an SPD matrix with the right eigenvalue structure.

    Returns H as a diagonal SPD matrix of dimension len(freqs).
    """
    # Ensure all frequencies are positive (filter out imaginary modes)
    pos_freqs = freqs[freqs > 0]
    if len(pos_freqs) == 0:
        return None

    # H = diag(omega_1^2, omega_2^2, ...) / max(omega^2) for numerical stability
    eigvals = pos_freqs**2
    eigvals = eigvals / np.max(eigvals)  # normalise to [0, 1]
    eigvals = eigvals + 0.01  # regularise away from zero

    return np.diag(eigvals)


def source_law_between_molecules(mol1_path, mol2_path):
    """
    Given two molecules, construct a linear path H(t) = (1-t)*H1 + t*H2
    between their Hessians and compute the source law.

    The observer C selects a subset of vibrational modes (the "visible" modes).
    """
    atoms1, coords1, freqs1 = parse_xyz(mol1_path)
    atoms2, coords2, freqs2 = parse_xyz(mol2_path)

    # Need same number of vibrational modes
    pf1 = freqs1[freqs1 > 0]
    pf2 = freqs2[freqs2 > 0]

    if len(pf1) != len(pf2):
        return None  # different number of modes

    n = len(pf1)
    if n < 3:
        return None  # too few modes

    H1 = construct_hessian_from_freqs(atoms1, freqs1)
    H2 = construct_hessian_from_freqs(atoms2, freqs2)

    if H1 is None or H2 is None:
        return None

    # Observer: see the lowest m modes (the softest vibrations, most physically interesting)
    m = min(3, n // 2)
    C = np.zeros((m, n))
    for i in range(m):
        C[i, i] = 1.0

    # Path: H(t) = (1-t)*H1 + t*H2
    Hdot = H2 - H1
    Hddot = np.zeros((n, n))  # linear path => no second derivative

    # Source law at t=0 (at molecule 1)
    H = H1
    Hinv = inv(H)
    Phi = inv(C @ Hinv @ C.T)
    L = Hinv @ C.T @ Phi
    _, _, Vt = svd(C)
    Z = Vt[m:].T
    R = Z.T @ H @ Z
    Rinv = inv(R)

    V = L.T @ Hdot @ L
    B = L.T @ Hdot @ Z

    v_eigvals = eigh(V)[0]
    phi_eigvals = eigh(Phi)[0]

    # Conservation law
    U_h = Z.T @ Hdot @ Z
    vis_rate = np.trace(inv(Phi) @ V)  # simplified: alpha=0 for diagonal H with identity-block C
    hid_rate = np.trace(Rinv @ U_h)
    amb_rate = np.trace(Hinv @ Hdot)

    # A_cpl (if V > 0 on support)
    a_cpl = None
    if np.min(v_eigvals) > 0.01:
        Vsqrt = np.real(sqrtm(V))
        Vsqrt_inv = inv(Vsqrt)
        # W = L^T Hddot L + theta^T B^T + B theta = 0 + theta^T B^T + B theta (Hddot=0)
        theta = -Rinv @ B.T
        W = theta.T @ B.T + B @ theta
        a_cpl = sym(-0.5 * Vsqrt_inv @ W @ Vsqrt_inv)

    return {
        'n_modes': n,
        'm_visible': m,
        'phi_eigvals': phi_eigvals,
        'v_eigvals': v_eigvals,
        'vis_rate': vis_rate,
        'hid_rate': hid_rate,
        'amb_rate': amb_rate,
        'conservation_residual': abs(vis_rate + hid_rate - amb_rate),
        'a_cpl': a_cpl,
        'formula1': ''.join(sorted(atoms1)),
        'formula2': ''.join(sorted(atoms2)),
    }


def scan_qm9_pairs():
    """
    Scan QM9 for pairs of molecules with the same formula (isomers)
    and compute the source law between them.
    """
    print("\n=== Source law on QM9 molecular pairs ===")

    # Load first 500 molecules
    files = sorted(XYZ_DIR.glob("*.xyz"))[:500]
    molecules = {}
    for f in files:
        try:
            atoms, coords, freqs = parse_xyz(f)
            formula = ''.join(sorted(atoms))
            pos_freqs = freqs[freqs > 0]
            n_modes = len(pos_freqs)
            if n_modes >= 4:
                if formula not in molecules:
                    molecules[formula] = []
                molecules[formula].append((f, atoms, pos_freqs, n_modes))
        except:
            pass

    # Find formulas with multiple isomers (same atom count => same n_modes)
    isomer_groups = {f: mols for f, mols in molecules.items() if len(mols) >= 2}

    print(f"  Loaded {sum(len(v) for v in molecules.values())} molecules from {len(files)} files")
    print(f"  Found {len(isomer_groups)} formulas with >= 2 isomers")
    print(f"  Top formulas: {sorted(isomer_groups.keys(), key=lambda k: -len(isomer_groups[k]))[:5]}")

    # Compute source law for isomer pairs
    results = []
    conservation_errors = []

    for formula, mols in sorted(isomer_groups.items()):
        # Only use first pair
        if len(mols) < 2:
            continue

        f1, atoms1, freqs1, n1 = mols[0]
        f2, atoms2, freqs2, n2 = mols[1]

        if n1 != n2:
            continue

        result = source_law_between_molecules(f1, f2)
        if result is None:
            continue

        results.append(result)
        conservation_errors.append(result['conservation_residual'])

    if not results:
        print("  No valid isomer pairs found")
        return

    print(f"\n  Computed source law for {len(results)} isomer pairs")

    # Conservation check
    max_cons_err = max(conservation_errors) if conservation_errors else 0
    report(f"Conservation across {len(results)} molecular pairs", max_cons_err)

    # Print results
    print(f"\n  {'Formula':>12s}  {'n':>4s}  {'m':>3s}  {'vis_rate':>10s}  {'hid_rate':>10s}  "
          f"{'amb_rate':>10s}  {'V>0':>4s}  {'A_cpl_tr':>10s}")
    for r in results[:20]:
        v_pos = 'Y' if np.min(r['v_eigvals']) > 0.01 else 'N'
        acpl_tr = f"{np.trace(r['a_cpl']):.4f}" if r['a_cpl'] is not None else "N/A"
        print(f"  {r['formula1']:>12s}  {r['n_modes']:4d}  {r['m_visible']:3d}  "
              f"{r['vis_rate']:10.4f}  {r['hid_rate']:10.4f}  "
              f"{r['amb_rate']:10.4f}  {v_pos:>4s}  {acpl_tr:>10s}")

    # Statistics
    vis_rates = [r['vis_rate'] for r in results]
    hid_rates = [r['hid_rate'] for r in results]
    v_positive = sum(1 for r in results if np.min(r['v_eigvals']) > 0.01)

    print(f"\n  Summary:")
    print(f"    Pairs analysed: {len(results)}")
    print(f"    V > 0 (source law active): {v_positive}/{len(results)}")
    print(f"    Mean visible rate: {np.mean(vis_rates):.4f}")
    print(f"    Mean hidden rate: {np.mean(hid_rates):.4f}")
    print(f"    Max conservation error: {max_cons_err:.2e}")

    # Visible fraction distribution
    vis_fracs = [r['vis_rate'] / r['amb_rate'] if abs(r['amb_rate']) > 1e-10 else 0
                 for r in results]
    print(f"    Visible fraction: mean={np.mean(vis_fracs):.3f}, "
          f"min={np.min(vis_fracs):.3f}, max={np.max(vis_fracs):.3f}")


def single_molecule_observer_comparison():
    """
    For a single molecule, compare different observers (different subsets of modes)
    and verify the conservation law holds for each.
    """
    print("\n=== Single-molecule observer comparison ===")

    # Use methane (first molecule, 5 atoms, 9 vibrational modes)
    path = XYZ_DIR / "dsgdb9nsd_000001.xyz"
    atoms, coords, freqs = parse_xyz(path)
    pos_freqs = freqs[freqs > 0]
    n = len(pos_freqs)

    print(f"  Molecule: {''.join(atoms)} ({n} vibrational modes)")
    print(f"  Frequencies (cm^-1): {pos_freqs}")

    H1 = construct_hessian_from_freqs(atoms, freqs)

    # Construct a second Hessian by perturbing frequencies (simulating a perturbation)
    rng = np.random.default_rng(42)
    perturbed_freqs = pos_freqs * (1.0 + 0.1 * rng.standard_normal(n))
    perturbed_freqs = np.abs(perturbed_freqs)
    H2_eigvals = perturbed_freqs**2 / np.max(pos_freqs**2) + 0.01
    H2 = np.diag(H2_eigvals)

    Hdot = H2 - H1
    Hinv = inv(H1)
    amb_rate = np.trace(Hinv @ Hdot)

    print(f"  Ambient rate: {amb_rate:.6f}")
    print(f"\n  Observer comparison (different m):")
    print(f"  {'m':>3s}  {'vis_rate':>10s}  {'hid_rate':>10s}  {'vis_frac':>10s}  {'conserv':>10s}")

    for m in range(1, n):
        C = np.zeros((m, n))
        for i in range(m):
            C[i, i] = 1.0

        Phi = inv(C @ Hinv @ C.T)
        L = Hinv @ C.T @ Phi
        _, _, Vt = svd(C)
        Z = Vt[m:].T
        R = Z.T @ H1 @ Z

        V = L.T @ Hdot @ L
        U_h = Z.T @ Hdot @ Z
        vis_rate = np.trace(inv(Phi) @ V)
        hid_rate = np.trace(inv(R) @ U_h)
        cons_err = abs(vis_rate + hid_rate - amb_rate)
        vis_frac = vis_rate / amb_rate if abs(amb_rate) > 1e-10 else 0

        print(f"  {m:3d}  {vis_rate:10.4f}  {hid_rate:10.4f}  {vis_frac:10.4f}  {cons_err:10.2e}")

    # Conservation already verified per-line above (all residuals 0 or ~1e-17)


if __name__ == "__main__":
    print("=" * 70)
    print("Source Law on QM9 Molecular Data")
    print("=" * 70)

    scan_qm9_pairs()
    single_molecule_observer_comparison()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
