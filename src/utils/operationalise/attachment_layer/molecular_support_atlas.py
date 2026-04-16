from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as stats
from collections import Counter
from pathlib import Path

import numpy as np

from molecular_vibrational_atlas import (
    DATA_ROOT,
    ELEMENT_SYMBOL,
    OUT_ROOT,
    Tolerances,
    formula_from_atomic_numbers,
    hessian_blocks_to_cartesian,
    iter_arrow_records,
    mass_vector,
    mass_weight_hessian,
    rigid_motion_basis,
    sym,
    vibrational_basis,
)


def inertia(values: np.ndarray, cutoff: float) -> tuple[int, int, int]:
    negative = int(np.sum(values < -cutoff))
    near_zero = int(np.sum(np.abs(values) <= cutoff))
    positive = int(np.sum(values > cutoff))
    return negative, near_zero, positive


def frequency_magnitudes(row: dict[str, object]) -> np.ndarray:
    freq = np.asarray(row["frequencies"], dtype=float)
    if freq.ndim == 1:
        return np.abs(freq)
    return np.sqrt(np.sum(freq * freq, axis=1))


def log_spectrum_correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    left = np.sort(np.asarray(left, dtype=float))
    right = np.sort(np.asarray(right, dtype=float))[-left.shape[0] :]
    keep = (left > 0.0) & (right > 0.0)
    if int(np.sum(keep)) < 3:
        return None
    return float(np.corrcoef(np.log(left[keep]), np.log(right[keep]))[0, 1])


def element_participation(mode: np.ndarray, atomic_numbers: np.ndarray) -> dict[str, float | str]:
    atom_power = np.sum(mode.reshape(-1, 3) ** 2, axis=1)
    total = float(np.sum(atom_power))
    if total <= 0.0:
        return {
            "soft_mode_heavy_share": math.nan,
            "soft_mode_hydrogen_share": math.nan,
            "soft_mode_hetero_share": math.nan,
            "soft_mode_top_element": "",
            "soft_mode_top_atom_share": math.nan,
        }
    heavy = float(np.sum(atom_power[atomic_numbers != 1]) / total)
    hydrogen = float(np.sum(atom_power[atomic_numbers == 1]) / total)
    hetero = float(np.sum(atom_power[(atomic_numbers != 1) & (atomic_numbers != 6)]) / total)
    top_atom = int(np.argmax(atom_power))
    top_z = int(atomic_numbers[top_atom])
    return {
        "soft_mode_heavy_share": heavy,
        "soft_mode_hydrogen_share": hydrogen,
        "soft_mode_hetero_share": hetero,
        "soft_mode_top_element": ELEMENT_SYMBOL.get(top_z, f"Z{top_z}"),
        "soft_mode_top_atom_share": float(atom_power[top_atom] / total),
    }


def compile_support_record(row: dict[str, object], split: str, shard: int, row_index: int, tolerances: Tolerances) -> dict[str, object]:
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    n_atoms = int(atomic_numbers.shape[0])
    masses = mass_vector(atomic_numbers)
    h_cart = hessian_blocks_to_cartesian(row["hessian"], n_atoms)
    h_mw, _inv_sqrt_m = mass_weight_hessian(h_cart, masses)
    rigid = rigid_motion_basis(positions, masses)
    vib = vibrational_basis(positions, masses)
    h_vib = sym(vib.T @ h_mw @ vib)
    raw_eigs = np.linalg.eigvalsh(h_mw)
    vib_eigs, vib_vecs = np.linalg.eigh(h_vib)
    raw_cutoff = max(tolerances.atol, tolerances.rtol * max(1.0, float(np.max(np.abs(raw_eigs)))))
    vib_cutoff = max(tolerances.atol, tolerances.rtol * max(1.0, float(np.max(np.abs(vib_eigs)))))
    raw_neg, raw_zero, raw_pos = inertia(raw_eigs, raw_cutoff)
    vib_neg, vib_zero, vib_pos = inertia(vib_eigs, vib_cutoff)

    h_norm = max(float(np.linalg.norm(h_mw, ord="fro")), tolerances.atol)
    rigid_block_norm = float(np.linalg.norm(rigid.T @ h_mw @ rigid, ord="fro"))
    rigid_cross_norm = float(np.linalg.norm(rigid.T @ h_mw @ vib, ord="fro"))
    vib_block_norm = float(np.linalg.norm(h_vib, ord="fro"))

    positive_vib = vib_eigs[vib_eigs > vib_cutoff]
    if positive_vib.size:
        softest_eig = float(np.min(positive_vib))
        stiffest_eig = float(np.max(positive_vib))
        condition = float(stiffest_eig / softest_eig)
        softness_logmean = float(-np.mean(np.log(positive_vib)))
    else:
        softest_eig = math.nan
        stiffest_eig = math.nan
        condition = math.nan
        softness_logmean = math.nan

    soft_idx = int(np.argmax(vib_eigs > vib_cutoff)) if positive_vib.size else -1
    if soft_idx >= 0:
        soft_mode = vib @ vib_vecs[:, soft_idx]
        mode_stats = element_participation(soft_mode, atomic_numbers)
    else:
        mode_stats = element_participation(np.zeros(3 * n_atoms), atomic_numbers)

    freq_mag = frequency_magnitudes(row)
    sqrt_vib = np.sqrt(np.maximum(positive_vib, 0.0))
    corr = log_spectrum_correlation(sqrt_vib, freq_mag)
    counts = Counter(int(z) for z in atomic_numbers.tolist())
    formula = formula_from_atomic_numbers(atomic_numbers)
    status = "vibrational_spd" if vib_neg == 0 and vib_zero == 0 and positive_vib.size == vib_eigs.size else "vibrational_indefinite_or_singular"
    if rigid.shape[1] == 5:
        rigid_geometry = "linear_or_rotationally_degenerate"
    elif rigid.shape[1] == 6:
        rigid_geometry = "nonlinear"
    else:
        rigid_geometry = "unexpected_rigid_rank"

    out = {
        "split": split,
        "shard": shard,
        "row_index": row_index,
        "label": str(row["label"]),
        "formula": formula,
        "n_atoms": n_atoms,
        "cartesian_dim": 3 * n_atoms,
        "rigid_rank": int(rigid.shape[1]),
        "rigid_geometry": rigid_geometry,
        "vibrational_dim": int(vib.shape[1]),
        "heavy_count": int(sum(v for z, v in counts.items() if z != 1)),
        "hydrogen_count": int(counts[1]),
        "hetero_count": int(sum(v for z, v in counts.items() if z not in (1, 6))),
        "raw_neg_count": raw_neg,
        "raw_zero_count": raw_zero,
        "raw_pos_count": raw_pos,
        "vib_neg_count": vib_neg,
        "vib_zero_count": vib_zero,
        "vib_pos_count": vib_pos,
        "raw_min_eig": float(np.min(raw_eigs)),
        "raw_max_eig": float(np.max(raw_eigs)),
        "vib_min_eig": float(np.min(vib_eigs)),
        "vib_max_eig": float(np.max(vib_eigs)),
        "softest_positive_eig": softest_eig,
        "stiffest_positive_eig": stiffest_eig,
        "vib_condition": condition,
        "softness_logmean": softness_logmean,
        "rigid_block_norm_ratio": rigid_block_norm / h_norm,
        "rigid_cross_norm_ratio": rigid_cross_norm / h_norm,
        "vib_block_norm_ratio": vib_block_norm / h_norm,
        "frequency_log_spectrum_corr": corr if corr is not None else "",
        "support_status": status,
    }
    out.update(mode_stats)
    return out


def summarize_rows(rows: list[dict[str, object]], split: str, shard: int, limit: int, output_prefix: str) -> dict[str, object]:
    def nums(key: str) -> list[float]:
        values = []
        for row in rows:
            value = row[key]
            if value == "" or value is None:
                continue
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                values.append(number)
        return values

    def quantiles(key: str) -> dict[str, float] | None:
        values = sorted(nums(key))
        if not values:
            return None
        return {
            "min": values[0],
            "q01": values[int(0.01 * (len(values) - 1))],
            "q10": values[int(0.10 * (len(values) - 1))],
            "median": stats.median(values),
            "q90": values[int(0.90 * (len(values) - 1))],
            "q99": values[int(0.99 * (len(values) - 1))],
            "max": values[-1],
        }

    status_counts = Counter(str(row["support_status"]) for row in rows)
    rigid_counts = Counter(str(row["rigid_geometry"]) for row in rows)
    top_elements = Counter(str(row["soft_mode_top_element"]) for row in rows)
    summary = {
        "schema_version": "molecular_support_atlas.v0",
        "source": "Hessian QM9 Figshare article 26363959 v4",
        "split": split,
        "shard": shard,
        "limit": limit,
        "rows": len(rows),
        "output_prefix": output_prefix,
        "status_counts": dict(status_counts),
        "rigid_geometry_counts": dict(rigid_counts),
        "soft_mode_top_element_counts": dict(top_elements),
        "vib_condition": quantiles("vib_condition"),
        "softest_positive_eig": quantiles("softest_positive_eig"),
        "rigid_block_norm_ratio": quantiles("rigid_block_norm_ratio"),
        "rigid_cross_norm_ratio": quantiles("rigid_cross_norm_ratio"),
        "frequency_log_spectrum_corr": quantiles("frequency_log_spectrum_corr"),
        "soft_mode_hydrogen_share": quantiles("soft_mode_hydrogen_share"),
        "soft_mode_hetero_share": quantiles("soft_mode_hetero_share"),
        "notes": [
            "This is a support/provenance atlas, not a hidden-load computation.",
            "The support_status field records whether the projected vibrational block is SPD under the declared tolerance.",
            "Soft-mode participation is computed in the mass-weighted Cartesian vector obtained by pulling the softest positive vibrational eigenvector back through the vibrational basis.",
        ],
    }
    return summary


def run(split: str, shard: int, limit: int, output_prefix: str) -> dict[str, object]:
    tolerances = Tolerances(atol=1e-8, rtol=1e-7)
    if not (DATA_ROOT / split / f"data-{shard:05d}-of-00005.arrow").exists():
        raise FileNotFoundError(f"missing extracted shard for split={split}, shard={shard}")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for idx, row in enumerate(iter_arrow_records(split, shard, limit)):
        try:
            rows.append(compile_support_record(row, split, shard, idx, tolerances))
        except Exception as exc:  # noqa: BLE001 - keep dataset scan moving, but record exact failure.
            failures.append({"row_index": idx, "label": str(row.get("label", f"row_{idx}")), "reason": str(exc)})

    if rows:
        csv_path = OUT_ROOT / f"{output_prefix}_support_rows.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    summary = summarize_rows(rows, split, shard, limit, output_prefix)
    summary["failures"] = len(failures)
    summary["failure_samples"] = failures[:20]
    (OUT_ROOT / f"{output_prefix}_support_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build molecular support/nullspace provenance summaries from Hessian QM9.")
    parser.add_argument("--split", default="vacuum", choices=["vacuum", "thf", "toluene", "water"])
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--output-prefix", default="hessian_qm9_vacuum_support")
    args = parser.parse_args()
    summary = run(args.split, args.shard, args.limit, args.output_prefix)
    print(json.dumps({k: summary[k] for k in ["rows", "failures", "status_counts", "rigid_geometry_counts"]}, indent=2))


if __name__ == "__main__":
    main()
