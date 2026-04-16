from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pyarrow.ipc as ipc
import scipy.linalg as la


def _ensure_nomogeo_src() -> None:
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "nomogeo"
        if candidate.is_dir():
            sys.path.insert(0, str(parent))
            return
        src_candidate = parent / "src" / "nomogeo"
        if src_candidate.is_dir():
            sys.path.insert(0, str(parent / "src"))
            return


_ensure_nomogeo_src()

from nomogeo import Tolerances, visible_precision  # noqa: E402
from nomogeo.exceptions import InputValidationError  # noqa: E402


BASE = Path(__file__).resolve().parent
DATA_ROOT = BASE / "incoming_data" / "hessian_qm9" / "hessian_qm9_DatasetDict"
OUT_ROOT = BASE / "outputs" / "molecular_atlas"

ATOMIC_MASS = {
    1: 1.00782503223,
    6: 12.0,
    7: 14.00307400443,
    8: 15.99491461957,
    9: 18.99840316273,
}

COVALENT_RADIUS_ANGSTROM = {
    1: 0.31,
    6: 0.76,
    7: 0.71,
    8: 0.66,
    9: 0.57,
}

ELEMENT_SYMBOL = {
    1: "H",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
}


def sym(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def iter_arrow_records(split: str, shard: int, limit: int) -> Iterable[dict[str, object]]:
    path = DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow"
    if not path.exists():
        raise FileNotFoundError(f"Missing Arrow shard: {path}")
    seen = 0
    with ipc.open_stream(path) as reader:
        for batch in reader:
            for row in batch.to_pylist():
                yield row
                seen += 1
                if seen >= limit:
                    return


def hessian_blocks_to_cartesian(blocks: object, n_atoms: int) -> np.ndarray:
    array = np.asarray(blocks, dtype=float)
    if array.shape == (n_atoms, 3, n_atoms, 3):
        return sym(array.reshape(3 * n_atoms, 3 * n_atoms))
    if array.shape == (n_atoms, n_atoms, 3, 3):
        return sym(array.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms))
    raise ValueError(
        f"expected Hessian shape {(n_atoms, 3, n_atoms, 3)} or {(n_atoms, n_atoms, 3, 3)}, got {array.shape}"
    )


def mass_vector(atomic_numbers: np.ndarray) -> np.ndarray:
    masses = []
    for atomic_number in atomic_numbers.tolist():
        if atomic_number not in ATOMIC_MASS:
            raise ValueError(f"unsupported atomic number {atomic_number}")
        masses.append(ATOMIC_MASS[int(atomic_number)])
    return np.asarray(masses, dtype=float)


def mass_weight_hessian(h_cart: np.ndarray, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    diag_m = np.repeat(masses, 3)
    inv_sqrt_m = 1.0 / np.sqrt(diag_m)
    return sym((inv_sqrt_m[:, None] * h_cart) * inv_sqrt_m[None, :]), inv_sqrt_m


def rigid_motion_basis(positions: np.ndarray, masses: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    weights = np.sqrt(masses)
    center = np.average(positions, axis=0, weights=masses)
    centered = positions - center

    columns: list[np.ndarray] = []
    for axis in range(3):
        vec = np.zeros((positions.shape[0], 3), dtype=float)
        vec[:, axis] = weights
        columns.append(vec.reshape(-1))

    axes = np.eye(3)
    for axis in axes:
        vec = np.cross(axis[None, :], centered) * weights[:, None]
        columns.append(vec.reshape(-1))

    raw = np.column_stack(columns)
    u, singular_values, _ = np.linalg.svd(raw, full_matrices=False)
    cutoff = max(tol, tol * max(1.0, float(np.max(singular_values))))
    rank = int(np.sum(singular_values > cutoff))
    return u[:, :rank]


def vibrational_basis(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    rigid = rigid_motion_basis(positions, masses)
    return la.null_space(rigid.T)


def compile_vibrational_object(row: dict[str, object], tolerances: Tolerances) -> dict[str, object]:
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    n_atoms = int(atomic_numbers.shape[0])
    masses = mass_vector(atomic_numbers)
    h_cart = hessian_blocks_to_cartesian(row["hessian"], n_atoms)
    h_mw, inv_sqrt_m = mass_weight_hessian(h_cart, masses)
    basis = vibrational_basis(positions, masses)
    h_vib = sym(basis.T @ h_mw @ basis)
    eigenvalues = np.linalg.eigvalsh(h_vib)
    cutoff = max(tolerances.atol, tolerances.rtol * max(1.0, float(np.max(np.abs(eigenvalues)))))
    min_eig = float(np.min(eigenvalues)) if eigenvalues.size else math.nan
    if not eigenvalues.size or min_eig <= cutoff:
        raise InputValidationError(
            f"vibrational Hessian is not SPD after rigid-mode projection; min_eig={min_eig:.6g}, cutoff={cutoff:.6g}"
        )
    return {
        "label": row["label"],
        "atomic_numbers": atomic_numbers,
        "positions": positions,
        "masses": masses,
        "inv_sqrt_m": inv_sqrt_m,
        "vibrational_basis": basis,
        "h_vib": h_vib,
        "n_atoms": n_atoms,
        "cartesian_dim": 3 * n_atoms,
        "vibrational_dim": int(h_vib.shape[0]),
        "rigid_rank": int(3 * n_atoms - h_vib.shape[0]),
        "min_eig": min_eig,
        "max_eig": float(np.max(eigenvalues)),
        "condition": float(np.max(eigenvalues) / min_eig),
        "formula": formula_from_atomic_numbers(atomic_numbers),
    }


def formula_from_atomic_numbers(atomic_numbers: np.ndarray) -> str:
    counts: dict[int, int] = {}
    for atomic_number in atomic_numbers.tolist():
        counts[int(atomic_number)] = counts.get(int(atomic_number), 0) + 1
    order = [6, 1, 7, 8, 9]
    parts = []
    for atomic_number in order:
        count = counts.pop(atomic_number, 0)
        if count:
            parts.append(ELEMENT_SYMBOL[atomic_number] + (str(count) if count != 1 else ""))
    for atomic_number in sorted(counts):
        count = counts[atomic_number]
        parts.append(f"Z{atomic_number}" + (str(count) if count != 1 else ""))
    return "".join(parts)


def infer_bonds(positions: np.ndarray, atomic_numbers: np.ndarray) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for i in range(positions.shape[0]):
        for j in range(i + 1, positions.shape[0]):
            zi = int(atomic_numbers[i])
            zj = int(atomic_numbers[j])
            ri = COVALENT_RADIUS_ANGSTROM.get(zi)
            rj = COVALENT_RADIUS_ANGSTROM.get(zj)
            if ri is None or rj is None:
                continue
            threshold = 1.25 * (ri + rj) + 0.10
            if float(np.linalg.norm(positions[i] - positions[j])) <= threshold:
                bonds.append((i, j))
    return bonds


def atom_selection_cartesian(indices: list[int], n_atoms: int) -> np.ndarray:
    rows = []
    for atom in indices:
        for axis in range(3):
            row = np.zeros(3 * n_atoms, dtype=float)
            row[3 * atom + axis] = 1.0
            rows.append(row)
    if not rows:
        return np.zeros((0, 3 * n_atoms), dtype=float)
    return np.vstack(rows)


def reduce_observer(c_vib: np.ndarray, tolerances: Tolerances) -> tuple[np.ndarray, int, int]:
    singular_values = np.linalg.svd(c_vib, compute_uv=False)
    cutoff = max(tolerances.atol, tolerances.rtol * max(1.0, float(np.max(singular_values)) if singular_values.size else 0.0))
    rank = int(np.sum(singular_values > cutoff))
    if rank == 0:
        raise InputValidationError("observer has zero rank after vibrational projection")
    _u, s, vt = np.linalg.svd(c_vib, full_matrices=False)
    return np.diag(s[:rank]) @ vt[:rank, :], rank, int(c_vib.shape[0])


def observer_from_atoms(compiled: dict[str, object], indices: list[int], tolerances: Tolerances) -> tuple[np.ndarray, int, int]:
    c_cart = atom_selection_cartesian(indices, int(compiled["n_atoms"]))
    if c_cart.shape[0] == 0:
        raise InputValidationError("empty atom selection")
    inv_sqrt_m = np.asarray(compiled["inv_sqrt_m"], dtype=float)
    basis = np.asarray(compiled["vibrational_basis"], dtype=float)
    c_vib = (c_cart * inv_sqrt_m[None, :]) @ basis
    return reduce_observer(c_vib, tolerances)


def observer_groups(compiled: dict[str, object]) -> list[dict[str, object]]:
    atomic_numbers = np.asarray(compiled["atomic_numbers"], dtype=int)
    positions = np.asarray(compiled["positions"], dtype=float)
    heavy = [int(i) for i, z in enumerate(atomic_numbers) if int(z) != 1]
    hydrogens = [int(i) for i, z in enumerate(atomic_numbers) if int(z) == 1]
    groups: list[dict[str, object]] = [
        {"observer": "heavy_atoms", "anchor_atom": "", "atoms": heavy},
        {"observer": "hydrogens", "anchor_atom": "", "atoms": hydrogens},
    ]
    bonds = infer_bonds(positions, atomic_numbers)
    for heavy_atom in heavy[:4]:
        attached_h = [
            j if i == heavy_atom else i
            for i, j in bonds
            if (i == heavy_atom and int(atomic_numbers[j]) == 1) or (j == heavy_atom and int(atomic_numbers[i]) == 1)
        ]
        groups.append(
            {
                "observer": "local_heavy_site_plus_attached_h",
                "anchor_atom": int(heavy_atom),
                "atoms": [heavy_atom] + sorted(set(int(x) for x in attached_h)),
            }
        )
    return groups


def read_observer(
    compiled: dict[str, object],
    observer_name: str,
    atom_indices: list[int],
    tolerances: Tolerances,
    anchor_atom: int | str = "",
) -> dict[str, object]:
    h_vib = np.asarray(compiled["h_vib"], dtype=float)
    c_reduced, rank, raw_dim = observer_from_atoms(compiled, atom_indices, tolerances)
    phi = visible_precision(h_vib, c_reduced, tolerances=tolerances)
    phi_eigs = np.linalg.eigvalsh(phi)
    return {
        "label": str(compiled["label"]),
        "formula": str(compiled["formula"]),
        "observer": observer_name,
        "anchor_atom": anchor_atom,
        "selected_atoms": atom_indices,
        "selected_atomic_numbers": [int(np.asarray(compiled["atomic_numbers"])[idx]) for idx in atom_indices],
        "raw_observer_dim": raw_dim,
        "observer_rank": rank,
        "visible_logdet": float(np.linalg.slogdet(phi)[1]),
        "visible_logdet_per_rank": float(np.linalg.slogdet(phi)[1] / rank),
        "visible_geom_mean_eig": float(np.exp(np.linalg.slogdet(phi)[1] / rank)),
        "visible_trace": float(np.trace(phi)),
        "visible_trace_per_rank": float(np.trace(phi) / rank),
        "visible_min_eig": float(np.min(phi_eigs)),
        "visible_max_eig": float(np.max(phi_eigs)),
        "visible_condition": float(np.max(phi_eigs) / np.min(phi_eigs)),
    }


def run(split: str, shard: int, limit: int, output_prefix: str) -> dict[str, object]:
    tolerances = Tolerances(atol=1e-8, rtol=1e-7)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    compile_refusals: list[dict[str, object]] = []
    observer_refusals: list[dict[str, object]] = []
    molecule_summaries: list[dict[str, object]] = []

    for idx, record in enumerate(iter_arrow_records(split, shard, limit)):
        try:
            compiled = compile_vibrational_object(record, tolerances)
        except Exception as exc:  # noqa: BLE001 - report and continue across a dataset batch.
            compile_refusals.append({"label": str(record.get("label", f"row_{idx}")), "reason": str(exc)})
            continue

        molecule_summaries.append(
            {
                "label": str(compiled["label"]),
                "n_atoms": compiled["n_atoms"],
                "cartesian_dim": compiled["cartesian_dim"],
                "vibrational_dim": compiled["vibrational_dim"],
                "rigid_rank": compiled["rigid_rank"],
                "min_eig": compiled["min_eig"],
                "max_eig": compiled["max_eig"],
                "condition": compiled["condition"],
                "formula": compiled["formula"],
                "atomic_numbers": np.asarray(compiled["atomic_numbers"], dtype=int).tolist(),
            }
        )

        for group in observer_groups(compiled):
            try:
                rows.append(
                    read_observer(
                        compiled,
                        str(group["observer"]),
                        list(group["atoms"]),
                        tolerances,
                        group["anchor_atom"],
                    )
                )
            except Exception as exc:  # noqa: BLE001 - an observer may legitimately project to zero rank.
                observer_refusals.append({"label": str(compiled["label"]), "observer": group["observer"], "reason": str(exc)})

    if rows:
        csv_path = OUT_ROOT / f"{output_prefix}_observer_readings.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "schema_version": "molecular_vibrational_atlas.v0",
        "source": "Hessian QM9 Figshare article 26363959 v4, vacuum split",
        "split": split,
        "shard": shard,
        "limit": limit,
        "compiled_molecules": len(molecule_summaries),
        "observer_readings": len(rows),
        "refusals": len(compile_refusals) + len(observer_refusals),
        "compile_refusals": len(compile_refusals),
        "observer_refusals": len(observer_refusals),
        "molecules": molecule_summaries,
        "compile_refusal_samples": compile_refusals[:20],
        "observer_refusal_samples": observer_refusals[:20],
        "refusal_samples": (compile_refusals + observer_refusals)[:20],
        "notes": [
            "Cartesian atom-selection observers are transformed into the mass-weighted vibrational coordinate chart before visible_precision is called.",
            "No hidden_load is emitted in v0 because atom-selection observers do not yet carry a declared ceiling convention.",
            "Bonded local-site observers use a simple covalent-radius heuristic and are therefore marked as generated observer candidates, not chemistry-certified functional groups.",
        ],
    }
    json_path = OUT_ROOT / f"{output_prefix}_summary.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile a Hessian QM9 shard into molecular observer-geometry readings.")
    parser.add_argument("--split", default="vacuum", choices=["vacuum", "thf", "toluene", "water"])
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--limit", type=int, default=25)
    parser.add_argument("--output-prefix", default="hessian_qm9_vacuum_shard0")
    args = parser.parse_args()
    summary = run(args.split, args.shard, args.limit, args.output_prefix)
    print(json.dumps({k: summary[k] for k in ["compiled_molecules", "observer_readings", "refusals"]}, indent=2))


if __name__ == "__main__":
    main()
