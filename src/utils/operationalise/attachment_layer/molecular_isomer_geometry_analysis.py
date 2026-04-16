from __future__ import annotations

import csv
import json
import math
import statistics as stats
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(__file__).resolve().parent
OUT_ROOT = BASE / "outputs" / "molecular_atlas"
SPLITS = ("vacuum", "thf", "toluene", "water")


def safe_float(value: str | int | float | None) -> float | None:
    if value in ("", None):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def log_value(value: str | int | float | None, floor: float = 1.0e-12) -> float | None:
    number = safe_float(value)
    if number is None:
        return None
    return math.log(max(number, floor))


def quantiles(values: list[float]) -> dict[str, float]:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return {}
    return {
        "min": clean[0],
        "q10": clean[int(0.10 * (len(clean) - 1))],
        "median": stats.median(clean),
        "q90": clean[int(0.90 * (len(clean) - 1))],
        "max": clean[-1],
    }


def load_observer_rows() -> dict[str, dict[str, object]]:
    molecules: dict[str, dict[str, object]] = defaultdict(lambda: {"local_site_gmeans": [], "local_site_conditions": []})
    for shard in range(5):
        path = OUT_ROOT / f"hessian_qm9_vacuum_shard{shard}_full_observer_readings.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                label = row["label"]
                mol = molecules[label]
                mol["label"] = label
                mol["formula"] = row["formula"]
                observer = row["observer"]
                if observer == "heavy_atoms":
                    mol["heavy_gmean"] = safe_float(row["visible_geom_mean_eig"])
                    mol["heavy_trace_per_rank"] = safe_float(row["visible_trace_per_rank"])
                    mol["heavy_condition"] = safe_float(row["visible_condition"])
                    mol["heavy_rank"] = safe_float(row["observer_rank"])
                elif observer == "hydrogens":
                    mol["hydrogen_gmean"] = safe_float(row["visible_geom_mean_eig"])
                    mol["hydrogen_trace_per_rank"] = safe_float(row["visible_trace_per_rank"])
                    mol["hydrogen_condition"] = safe_float(row["visible_condition"])
                    mol["hydrogen_rank"] = safe_float(row["observer_rank"])
                elif observer == "local_heavy_site_plus_attached_h":
                    value = safe_float(row["visible_geom_mean_eig"])
                    condition = safe_float(row["visible_condition"])
                    if value is not None:
                        mol["local_site_gmeans"].append(value)
                    if condition is not None:
                        mol["local_site_conditions"].append(condition)
    return molecules


def add_support_patterns(molecules: dict[str, dict[str, object]]) -> None:
    path = OUT_ROOT / "molecular_support_label_patterns_v0.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            label = row["label"]
            if label not in molecules:
                continue
            mol = molecules[label]
            mol["support_pattern"] = row["pattern"]
            for split in SPLITS:
                mol[f"{split}_softest_positive_eig"] = safe_float(row[f"{split}_softest_positive_eig"])
                mol[f"{split}_vib_min_eig"] = safe_float(row[f"{split}_vib_min_eig"])


def finalize_molecule_features(molecules: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mol in molecules.values():
        if "heavy_gmean" not in mol:
            continue
        local = sorted(float(x) for x in mol["local_site_gmeans"])
        if local:
            mol["local_site_count"] = len(local)
            mol["local_gmean_min"] = min(local)
            mol["local_gmean_median"] = stats.median(local)
            mol["local_gmean_max"] = max(local)
            mol["local_gmean_log_range"] = math.log(max(max(local), 1.0e-12)) - math.log(max(min(local), 1.0e-12))
        else:
            mol["local_site_count"] = 0
            mol["local_gmean_min"] = None
            mol["local_gmean_median"] = None
            mol["local_gmean_max"] = None
            mol["local_gmean_log_range"] = None

        hv = mol.get("heavy_gmean")
        hy = mol.get("hydrogen_gmean")
        mol["heavy_hydrogen_log_ratio"] = None if hv is None or hy in (None, 0.0) else math.log(max(float(hv), 1.0e-12) / max(float(hy), 1.0e-12))
        vac = mol.get("vacuum_softest_positive_eig")
        for split in ("thf", "toluene", "water"):
            other = mol.get(f"{split}_softest_positive_eig")
            mol[f"{split}_soft_vacuum_log_ratio"] = (
                None if vac in (None, 0.0) or other in (None, 0.0) else math.log(max(float(other), 1.0e-12) / max(float(vac), 1.0e-12))
            )
        rows.append(mol)
    return rows


def write_molecule_fingerprints(rows: list[dict[str, object]]) -> None:
    fields = [
        "label",
        "formula",
        "support_pattern",
        "heavy_gmean",
        "hydrogen_gmean",
        "heavy_hydrogen_log_ratio",
        "local_site_count",
        "local_gmean_min",
        "local_gmean_median",
        "local_gmean_max",
        "local_gmean_log_range",
        "vacuum_softest_positive_eig",
        "thf_soft_vacuum_log_ratio",
        "toluene_soft_vacuum_log_ratio",
        "water_soft_vacuum_log_ratio",
    ]
    path = OUT_ROOT / "molecular_isomer_fingerprints_v0.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def formula_summary(rows: list[dict[str, object]], min_count: int) -> list[dict[str, object]]:
    by_formula: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_formula[str(row["formula"])].append(row)
    summaries: list[dict[str, object]] = []
    for formula, group in by_formula.items():
        if len(group) < min_count:
            continue
        ratios = [float(row["heavy_hydrogen_log_ratio"]) for row in group if row.get("heavy_hydrogen_log_ratio") is not None]
        local_ranges = [float(row["local_gmean_log_range"]) for row in group if row.get("local_gmean_log_range") is not None]
        water_ratios = [float(row["water_soft_vacuum_log_ratio"]) for row in group if row.get("water_soft_vacuum_log_ratio") is not None]
        support = Counter(str(row.get("support_pattern", "")) for row in group)
        support_entropy = 0.0
        for count in support.values():
            p = count / len(group)
            support_entropy -= p * math.log(p)
        ratio_q = quantiles(ratios)
        local_q = quantiles(local_ranges)
        water_q = quantiles(water_ratios)
        summaries.append(
            {
                "formula": formula,
                "n": len(group),
                "support_pattern_count": len(support),
                "support_entropy": support_entropy,
                "dominant_support_pattern": support.most_common(1)[0][0],
                "heavy_hydrogen_log_ratio_spread": ratio_q.get("q90", 0.0) - ratio_q.get("q10", 0.0),
                "heavy_hydrogen_log_ratio_median": ratio_q.get("median", ""),
                "local_log_range_median": local_q.get("median", ""),
                "water_soft_log_ratio_median": water_q.get("median", ""),
            }
        )
    return sorted(
        summaries,
        key=lambda item: (
            float(item["support_entropy"]),
            float(item["heavy_hydrogen_log_ratio_spread"]),
            float(item["n"]),
        ),
        reverse=True,
    )


def feature_vector(row: dict[str, object], keys: tuple[str, ...]) -> tuple[float, ...] | None:
    out = []
    for key in keys:
        if key.endswith("_log"):
            source = key[:-4]
            value = log_value(row.get(source))
        else:
            value = safe_float(row.get(key))
        if value is None:
            return None
        out.append(float(value))
    return tuple(out)


def distance(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    return math.sqrt(sum((a - b) * (a - b) for a, b in zip(left, right, strict=True)))


def coarse_collapse_candidates(rows: list[dict[str, object]], min_count: int) -> list[dict[str, object]]:
    by_formula: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_formula[str(row["formula"])].append(row)

    coarse_keys = ("heavy_gmean_log", "heavy_condition_log")
    rich_keys = (
        "heavy_gmean_log",
        "hydrogen_gmean_log",
        "local_gmean_median_log",
        "local_gmean_log_range",
        "vacuum_softest_positive_eig_log",
        "water_soft_vacuum_log_ratio",
    )
    candidates: list[dict[str, object]] = []
    for formula, group in by_formula.items():
        if len(group) < min_count:
            continue
        vectors = []
        for row in group:
            coarse = feature_vector(row, coarse_keys)
            rich = feature_vector(row, rich_keys)
            if coarse is None or rich is None:
                continue
            vectors.append((row, coarse, rich))
        if len(vectors) < 2:
            continue
        vectors.sort(key=lambda item: item[1][0])
        best = None
        for idx, (row, coarse, rich) in enumerate(vectors):
            for other_row, other_coarse, other_rich in vectors[max(0, idx - 12) : min(len(vectors), idx + 13)]:
                if row is other_row:
                    continue
                coarse_dist = distance(coarse, other_coarse)
                rich_dist = distance(rich, other_rich)
                score = rich_dist / max(coarse_dist, 1.0e-6)
                if best is None or score > best["score"]:
                    best = {
                        "formula": formula,
                        "n": len(group),
                        "left_label": row["label"],
                        "right_label": other_row["label"],
                        "left_pattern": row.get("support_pattern", ""),
                        "right_pattern": other_row.get("support_pattern", ""),
                        "coarse_distance": coarse_dist,
                        "rich_distance": rich_dist,
                        "score": score,
                        "left_heavy_gmean": row.get("heavy_gmean", ""),
                        "right_heavy_gmean": other_row.get("heavy_gmean", ""),
                        "left_local_median": row.get("local_gmean_median", ""),
                        "right_local_median": other_row.get("local_gmean_median", ""),
                        "left_water_soft_ratio": row.get("water_soft_vacuum_log_ratio", ""),
                        "right_water_soft_ratio": other_row.get("water_soft_vacuum_log_ratio", ""),
                    }
        if best is not None:
            candidates.append(best)
    return sorted(candidates, key=lambda item: float(item["score"]), reverse=True)


def write_formula_summary(items: list[dict[str, object]]) -> None:
    path = OUT_ROOT / "molecular_isomer_formula_summary_v0.csv"
    fields = [
        "formula",
        "n",
        "support_pattern_count",
        "support_entropy",
        "dominant_support_pattern",
        "heavy_hydrogen_log_ratio_spread",
        "heavy_hydrogen_log_ratio_median",
        "local_log_range_median",
        "water_soft_log_ratio_median",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in items:
            writer.writerow({field: item.get(field, "") for field in fields})


def write_candidates(items: list[dict[str, object]]) -> None:
    path = OUT_ROOT / "molecular_isomer_coarse_collapse_candidates_v0.csv"
    fields = [
        "formula",
        "n",
        "left_label",
        "right_label",
        "left_pattern",
        "right_pattern",
        "coarse_distance",
        "rich_distance",
        "score",
        "left_heavy_gmean",
        "right_heavy_gmean",
        "left_local_median",
        "right_local_median",
        "left_water_soft_ratio",
        "right_water_soft_ratio",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in items:
            writer.writerow({field: item.get(field, "") for field in fields})


def write_markdown(formulas: list[dict[str, object]], candidates: list[dict[str, object]], rows: list[dict[str, object]]) -> None:
    pattern_counts = Counter(str(row.get("support_pattern", "")) for row in rows)
    lines = [
        "# Molecular Isomer Geometry Analysis v0",
        "",
        "Status: generated from the full vacuum observer atlas and the four-environment molecular support atlas.",
        "",
        "This is not a theorem-level false-collapse result. Formula is a coarse descriptor, not a module observer map. The analysis asks a narrower empirical question: when formula is held fixed, does observer/support geometry still vary enough to create useful lead sets?",
        "",
        "Because the observer fingerprints come from the vacuum `visible_precision` atlas, this report only includes molecules admitted by the vacuum vibrational-SPD gate. The separate support analysis covers all records and all 16 environment support patterns.",
        "",
        "## Scope",
        "",
        f"- molecule fingerprints with vacuum observer geometry: `{len(rows)}`",
        f"- formulas with at least 20 records: `{len(formulas)}`",
        f"- support patterns represented inside this vacuum-admitted subset: `{len(pattern_counts)}`",
        "",
        "## Formula Groups With Rich Internal Structure",
        "",
        "Ranked by support-pattern entropy, then heavy/hydrogen observer-ratio spread.",
        "",
        "| Formula | n | support patterns | entropy | dominant pattern | heavy/H log-ratio spread | local range median | water soft log-ratio median |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for item in formulas[:20]:
        lines.append(
            "| {formula} | {n} | {support_pattern_count} | {support_entropy:.3f} | {dominant_support_pattern} | {heavy_hydrogen_log_ratio_spread:.3f} | {local_log_range_median:.3f} | {water_soft_log_ratio_median:.3f} |".format(
                **item
            )
        )

    lines.extend(
        [
            "",
            "## Coarse-Collapse Candidates",
            "",
            "These are same-formula pairs that are close under a coarse heavy-atom observer fingerprint but separated by a richer fingerprint using hydrogen, local-site, and water-softening coordinates. They are lead pairs for later structural inspection, not final claims.",
            "",
            "| Formula | Pair | Patterns | coarse distance | rich distance | score | local medians | water soft ratios |",
            "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for item in candidates[:20]:
        lines.append(
            "| {formula} | {left_label} / {right_label} | {left_pattern} / {right_pattern} | {coarse_distance:.4g} | {rich_distance:.4g} | {score:.3g} | {left_local_median:.4g} / {right_local_median:.4g} | {left_water_soft_ratio:.4g} / {right_water_soft_ratio:.4g} |".format(
                **item
            )
        )

    lines.extend(
        [
            "",
            "## Scientific Read",
            "",
            "- Same formula does not fix observer/support geometry. Several large formula classes have multiple support patterns and substantial observer-ratio spread.",
            "- The best next concrete lead is not arbitrary clustering; it is pair inspection inside high-entropy formula groups, where formula-level collapse is strongest and observer/support separation is largest.",
            "- The candidate pairs should be joined to molecular graph/structure data next. Without graph inspection, we cannot say whether the separation corresponds to known isomer motifs, conformer differences, or numerical quirks.",
            "- This strengthens the molecular atlas programme because it shows the reader is not merely reproducing elemental composition. It is finding within-formula geometry that needs explanation.",
            "",
            "## Boundaries",
            "",
            "- No chemical identity or functional group classifier was used.",
            "- Formula is not a linear observer, so `false collapse` remains an analogy here.",
            "- Pair distances are exploratory fingerprints, not invariant metrics from the module theory.",
            "- The support side is static across environments, not a continuous support-event path.",
        ]
    )
    (OUT_ROOT / "molecular_isomer_geometry_analysis_v0.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    molecules = load_observer_rows()
    add_support_patterns(molecules)
    rows = finalize_molecule_features(molecules)
    write_molecule_fingerprints(rows)
    formulas = formula_summary(rows, min_count=20)
    candidates = coarse_collapse_candidates(rows, min_count=20)
    write_formula_summary(formulas)
    write_candidates(candidates)
    report = {
        "schema_version": "molecular_isomer_geometry_analysis.v0",
        "fingerprints": len(rows),
        "formula_summaries": len(formulas),
        "coarse_collapse_candidates": len(candidates),
        "top_formula_summaries": formulas[:25],
        "top_candidates": candidates[:25],
    }
    (OUT_ROOT / "molecular_isomer_geometry_analysis_v0.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(formulas, candidates, rows)
    print(json.dumps({k: report[k] for k in ["fingerprints", "formula_summaries", "coarse_collapse_candidates"]}, indent=2))


if __name__ == "__main__":
    main()
