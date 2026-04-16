from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import scipy.linalg as la

from molecular_vibrational_atlas import (
    DATA_ROOT,
    OUT_ROOT,
    Tolerances,
    compile_vibrational_object,
    hessian_blocks_to_cartesian,
    mass_vector,
    mass_weight_hessian,
    observer_from_atoms,
    observer_groups,
    sym,
)


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

from nomogeo import closure_scores  # noqa: E402


def shard_path(split: str, shard: int) -> Path:
    return DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow"


def label_set(split: str, shard: int) -> set[str]:
    labels: set[str] = set()
    with ipc.open_stream(shard_path(split, shard)) as reader:
        for batch in reader:
            labels.update(str(label) for label in batch.column("label").to_pylist())
    return labels


def records_by_label(split: str, shard: int, labels: set[str]) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    with ipc.open_stream(shard_path(split, shard)) as reader:
        for batch in reader:
            batch_labels = [str(label) for label in batch.column("label").to_pylist()]
            indices = [idx for idx, label in enumerate(batch_labels) if label in labels]
            if indices:
                selected = batch.take(pa.array(indices))
                for row in selected.to_pylist():
                    records[str(row["label"])] = row
    return records


def kabsch_rotation(source: np.ndarray, target: np.ndarray, masses: np.ndarray) -> tuple[np.ndarray, float]:
    source_center = np.average(source, axis=0, weights=masses)
    target_center = np.average(target, axis=0, weights=masses)
    x = source - source_center
    y = target - target_center
    w = np.sqrt(masses)[:, None]
    covariance = (x * w).T @ (y * w)
    u, _s, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1.0
        rotation = u @ vt
    residual = (x @ rotation) - y
    rmsd = float(np.sqrt(np.sum(masses[:, None] * residual * residual) / np.sum(masses)))
    return rotation, rmsd


def solvent_hessian_in_reference_vibrational_chart(
    solvent_row: dict[str, object],
    reference: dict[str, object],
) -> tuple[np.ndarray, float]:
    atomic_numbers = np.asarray(solvent_row["atomic_numbers"], dtype=int)
    reference_atomic_numbers = np.asarray(reference["atomic_numbers"], dtype=int)
    if atomic_numbers.tolist() != reference_atomic_numbers.tolist():
        raise ValueError("atomic numbers do not match across solvent records")

    masses = mass_vector(atomic_numbers)
    positions = np.asarray(solvent_row["positions"], dtype=float)
    reference_positions = np.asarray(reference["positions"], dtype=float)
    rotation, rmsd = kabsch_rotation(positions, reference_positions, masses)

    h_cart = hessian_blocks_to_cartesian(solvent_row["hessian"], int(atomic_numbers.shape[0]))
    h_mw, _inv_sqrt_m = mass_weight_hessian(h_cart, masses)
    transform = np.kron(np.eye(int(atomic_numbers.shape[0])), rotation.T)
    h_aligned = sym(transform @ h_mw @ transform.T)
    vib_basis = np.asarray(reference["vibrational_basis"], dtype=float)
    return sym(vib_basis.T @ h_aligned @ vib_basis), rmsd


def observer_plane_from_c(H: np.ndarray, C: np.ndarray, tolerances: Tolerances) -> np.ndarray:
    eigvals, eigvecs = la.eigh(H)
    cutoff = max(tolerances.atol, tolerances.rtol * max(1.0, float(np.max(np.abs(eigvals)))))
    if float(np.min(eigvals)) <= cutoff:
        raise ValueError("H must be SPD for observer-plane conversion")
    h_inv_half = sym((eigvecs * (1.0 / np.sqrt(eigvals))) @ eigvecs.T)
    raw = h_inv_half @ C.T
    u, singular_values, _vt = np.linalg.svd(raw, full_matrices=False)
    rank = int(np.sum(singular_values > cutoff))
    if rank == 0:
        raise ValueError("observer has zero H-whitened rank")
    return u[:, :rank]


def run(limit: int, shard: int, output_prefix: str) -> dict[str, object]:
    tolerances = Tolerances(atol=1e-8, rtol=1e-7)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    common_labels = sorted(set.intersection(*(label_set(split, shard) for split in ["vacuum", "thf", "toluene", "water"])))[:limit]
    target_labels = set(common_labels)
    split_maps = {split: records_by_label(split, shard, target_labels) for split in ["vacuum", "thf", "toluene", "water"]}

    rows: list[dict[str, object]] = []
    refusals: list[dict[str, object]] = []
    alignment_rmsd: list[float] = []

    for label in common_labels:
        try:
            reference = compile_vibrational_object(split_maps["vacuum"][label], tolerances)
            h_ref = np.asarray(reference["h_vib"], dtype=float)
            perturbations = []
            rmsds = {}
            for split in ["thf", "toluene", "water"]:
                h_solvent, rmsd = solvent_hessian_in_reference_vibrational_chart(split_maps[split][label], reference)
                perturbations.append(sym(h_solvent - h_ref))
                rmsds[split] = rmsd
                alignment_rmsd.append(rmsd)

            for group in observer_groups(reference):
                if group["observer"] == "local_heavy_site_plus_attached_h" and int(group["anchor_atom"]) > 1:
                    continue
                c_reduced, rank, raw_dim = observer_from_atoms(reference, list(group["atoms"]), tolerances)
                plane = observer_plane_from_c(h_ref, c_reduced, tolerances)
                scores = closure_scores(h_ref, perturbations, plane, tolerances=tolerances)
                rows.append(
                    {
                        "label": label,
                        "formula": str(reference["formula"]),
                        "observer": group["observer"],
                        "anchor_atom": group["anchor_atom"],
                        "raw_observer_dim": raw_dim,
                        "observer_rank": rank,
                        "thf_alignment_rmsd": rmsds["thf"],
                        "toluene_alignment_rmsd": rmsds["toluene"],
                        "water_alignment_rmsd": rmsds["water"],
                        "solvent_leakage": scores.leakage,
                        "solvent_visible_score": scores.visible_score,
                        "solvent_eta": scores.eta,
                        "solvent_total_curvature": scores.total_curvature,
                    }
                )
        except Exception as exc:  # noqa: BLE001 - this is an empirical admission/refusal pass.
            refusals.append({"label": label, "reason": str(exc)})

    if rows:
        csv_path = OUT_ROOT / f"{output_prefix}_solvent_closure.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "schema_version": "molecular_solvent_perturbation_probe.v0",
        "source": "Hessian QM9 Figshare article 26363959 v4; first shard across vacuum, thf, toluene, water",
        "declaration": "Each solvent Hessian is mass-weighted, mass-weighted-Kabsch aligned to the vacuum geometry, and pulled back to the vacuum vibrational chart before forming solvent perturbations.",
        "limit": limit,
        "shard": shard,
        "common_labels": len(common_labels),
        "solvent_closure_rows": len(rows),
        "refusals": len(refusals),
        "alignment_rmsd_median": float(np.median(alignment_rmsd)) if alignment_rmsd else None,
        "alignment_rmsd_max": float(np.max(alignment_rmsd)) if alignment_rmsd else None,
        "refusal_samples": refusals[:20],
        "notes": [
            "This is a declared same-label perturbation probe, not a claim that solvent-optimized geometries share a canonical global chart without alignment.",
            "The pullback convention should be stress-tested before using solvent leakage as chemistry-facing evidence.",
            "No hidden_load is emitted; this probe only licenses closure/leakage of a declared perturbation family.",
        ],
    }
    (OUT_ROOT / f"{output_prefix}_solvent_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Probe solvent perturbation leakage on Hessian QM9.")
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--output-prefix", default="hessian_qm9_solvent_shard0")
    args = parser.parse_args()
    summary = run(args.limit, args.shard, args.output_prefix)
    print(json.dumps({k: summary[k] for k in ["common_labels", "solvent_closure_rows", "refusals", "alignment_rmsd_median", "alignment_rmsd_max"]}, indent=2))


if __name__ == "__main__":
    main()
