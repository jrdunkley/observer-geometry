from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow.ipc as ipc

from molecular_vibrational_atlas import DATA_ROOT, ELEMENT_SYMBOL, OUT_ROOT, infer_bonds


BASE = Path(__file__).resolve().parent
TOP_PAIR_LIMIT = 25
WATER_REFUSAL_LIMIT = 50


def safe_float(value: str | int | float | None) -> float | None:
    if value in ("", None):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def symbol(atomic_number: int) -> str:
    return ELEMENT_SYMBOL.get(int(atomic_number), f"Z{int(atomic_number)}")


def compact_counter(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter))


def graph_components(nodes: list[int], edges: list[tuple[int, int]]) -> int:
    if not nodes:
        return 0
    node_set = set(nodes)
    adjacency = {node: set() for node in nodes}
    for i, j in edges:
        if i in node_set and j in node_set:
            adjacency[i].add(j)
            adjacency[j].add(i)
    seen: set[int] = set()
    components = 0
    for node in nodes:
        if node in seen:
            continue
        components += 1
        stack = [node]
        seen.add(node)
        while stack:
            current = stack.pop()
            for nxt in adjacency[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
    return components


def graph_signature(row: dict[str, object]) -> dict[str, object]:
    atomic_numbers = np.asarray(row["atomic_numbers"], dtype=int)
    positions = np.asarray(row["positions"], dtype=float)
    bonds = infer_bonds(positions, atomic_numbers)
    n_atoms = int(atomic_numbers.shape[0])

    atom_counts = Counter(symbol(int(z)) for z in atomic_numbers.tolist())
    bond_counts: Counter[str] = Counter()
    degree = Counter()
    heavy_degree = Counter()
    attached_h = Counter()
    hetero_neighbours = Counter()

    for i, j in bonds:
        si = symbol(int(atomic_numbers[i]))
        sj = symbol(int(atomic_numbers[j]))
        bond_counts["-".join(sorted((si, sj)))] += 1
        degree[i] += 1
        degree[j] += 1

    heavy_nodes = [int(i) for i, z in enumerate(atomic_numbers) if int(z) != 1]
    heavy_bonds = [(i, j) for i, j in bonds if int(atomic_numbers[i]) != 1 and int(atomic_numbers[j]) != 1]
    for i, j in heavy_bonds:
        heavy_degree[i] += 1
        heavy_degree[j] += 1
        if int(atomic_numbers[j]) not in (1, 6):
            hetero_neighbours[i] += 1
        if int(atomic_numbers[i]) not in (1, 6):
            hetero_neighbours[j] += 1

    for i, j in bonds:
        if int(atomic_numbers[i]) != 1 and int(atomic_numbers[j]) == 1:
            attached_h[i] += 1
        if int(atomic_numbers[j]) != 1 and int(atomic_numbers[i]) == 1:
            attached_h[j] += 1

    local_env: Counter[str] = Counter()
    heavy_degree_hist: Counter[str] = Counter()
    for node in heavy_nodes:
        z = int(atomic_numbers[node])
        hd = int(heavy_degree[node])
        ah = int(attached_h[node])
        hn = int(hetero_neighbours[node])
        local_env[f"{symbol(z)}:hd{hd}:H{ah}:het{hn}"] += 1
        heavy_degree_hist[str(hd)] += 1

    total_components = graph_components(list(range(n_atoms)), bonds)
    heavy_components = graph_components(heavy_nodes, heavy_bonds)
    total_cyclomatic = len(bonds) - n_atoms + total_components
    heavy_cyclomatic = len(heavy_bonds) - len(heavy_nodes) + heavy_components if heavy_nodes else 0

    return {
        "label": str(row["label"]),
        "atom_count_signature": compact_counter(atom_counts),
        "bond_signature": compact_counter(bond_counts),
        "local_env_signature": compact_counter(local_env),
        "heavy_degree_hist": compact_counter(heavy_degree_hist),
        "n_atoms": n_atoms,
        "n_bonds": len(bonds),
        "heavy_atoms": len(heavy_nodes),
        "heavy_bonds": len(heavy_bonds),
        "total_components": total_components,
        "heavy_components": heavy_components,
        "total_cyclomatic": total_cyclomatic,
        "heavy_cyclomatic": heavy_cyclomatic,
    }


def load_candidate_pairs(limit: int) -> list[dict[str, str]]:
    path = OUT_ROOT / "molecular_isomer_coarse_collapse_candidates_v0.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))[:limit]


def load_water_refusals(limit: int) -> list[dict[str, str]]:
    path = OUT_ROOT / "molecular_support_label_patterns_v0.csv"
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["pattern"] == "SSSI":
                rows.append(row)
    rows.sort(key=lambda row: safe_float(row.get("water_vib_min_eig")) or 0.0)
    return rows[:limit]


def collect_vacuum_rows(labels: set[str]) -> dict[str, dict[str, object]]:
    found: dict[str, dict[str, object]] = {}
    for shard in range(5):
        path = DATA_ROOT / "vacuum" / f"data-{shard:05d}-of-00005.arrow"
        with ipc.open_stream(path) as reader:
            for batch in reader:
                slim = batch.select(["label", "atomic_numbers", "positions"])
                for row in slim.to_pylist():
                    label = str(row["label"])
                    if label in labels:
                        found[label] = row
        if len(found) == len(labels):
            break
    return found


def pair_rows(candidates: list[dict[str, str]], signatures: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in candidates:
        left = signatures.get(item["left_label"])
        right = signatures.get(item["right_label"])
        if left is None or right is None:
            continue
        graph_same = (
            left["bond_signature"] == right["bond_signature"]
            and left["local_env_signature"] == right["local_env_signature"]
            and left["heavy_cyclomatic"] == right["heavy_cyclomatic"]
        )
        rows.append(
            {
                **item,
                "graph_same_by_signature": graph_same,
                "left_bond_signature": left["bond_signature"],
                "right_bond_signature": right["bond_signature"],
                "left_local_env_signature": left["local_env_signature"],
                "right_local_env_signature": right["local_env_signature"],
                "left_heavy_degree_hist": left["heavy_degree_hist"],
                "right_heavy_degree_hist": right["heavy_degree_hist"],
                "left_heavy_cyclomatic": left["heavy_cyclomatic"],
                "right_heavy_cyclomatic": right["heavy_cyclomatic"],
            }
        )
    return rows


def water_rows(refusals: list[dict[str, str]], signatures: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in refusals:
        sig = signatures.get(item["label"])
        if sig is None:
            continue
        rows.append(
            {
                "label": item["label"],
                "formula": item["formula"],
                "pattern": item["pattern"],
                "water_vib_min_eig": item["water_vib_min_eig"],
                "vacuum_softest_positive_eig": item["vacuum_softest_positive_eig"],
                "thf_softest_positive_eig": item["thf_softest_positive_eig"],
                "toluene_softest_positive_eig": item["toluene_softest_positive_eig"],
                "bond_signature": sig["bond_signature"],
                "local_env_signature": sig["local_env_signature"],
                "heavy_degree_hist": sig["heavy_degree_hist"],
                "heavy_cyclomatic": sig["heavy_cyclomatic"],
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_markdown(pair_items: list[dict[str, object]], water_items: list[dict[str, object]], missing: list[str]) -> None:
    same_graph = sum(1 for row in pair_items if row["graph_same_by_signature"])
    water_formulas = Counter(str(row["formula"]) for row in water_items)
    water_local = Counter(str(row["local_env_signature"]) for row in water_items)

    lines = [
        "# Molecular Graph Lead Inspection v0",
        "",
        "Status: generated from top observer-geometry lead pairs and water-only support refusals, joined back to the vacuum molecular graph using the existing covalent-radius heuristic.",
        "",
        "This is a graph-inspection layer, not a chemical classifier. The bond graph is heuristic and is used to triage leads before any chemistry-facing claim.",
        "",
        "## Scope",
        "",
        f"- top coarse-collapse candidate pairs inspected: `{len(pair_items)}`",
        f"- top water-only support refusals inspected: `{len(water_items)}`",
        f"- labels missing from vacuum Arrow scan: `{len(missing)}`",
        "",
        "## Pair Inspection",
        "",
        f"Among the inspected same-formula candidate pairs, `{same_graph}` of `{len(pair_items)}` have the same simple graph signature under bond/local-environment/cyclomatic comparison. These are the strongest leads for geometry separating molecules that are not already trivially separated by this graph heuristic.",
        "",
        "| Formula | Pair | graph same | score | bond signatures | local environment signatures |",
        "| --- | --- | ---: | ---: | --- | --- |",
    ]
    for row in pair_items[:10]:
        score = safe_float(row.get("score"))
        display = dict(row)
        display["score_display"] = "" if score is None else f"{score:.4g}"
        lines.append(
            "| {formula} | {left_label} / {right_label} | {graph_same_by_signature} | {score_display} | {left_bond_signature} // {right_bond_signature} | {left_local_env_signature} // {right_local_env_signature} |".format(
                **display,
            )
        )

    lines.extend(
        [
            "",
            "## Water-Only Refusal Inspection",
            "",
            "The `SSSI` pattern remains the cleanest environment-specific refusal lead: admitted in vacuum, THF, and toluene, refused in water under the declared support gate.",
            "",
            "Top formulas among the inspected water-only refusals:",
            "",
            "| Formula | Count |",
            "| --- | ---: |",
        ]
    )
    for formula, count in water_formulas.most_common(12):
        lines.append(f"| {formula} | {count} |")

    lines.extend(
        [
            "",
            "Most common local graph signatures among the inspected water-only refusals:",
            "",
            "| Local signature | Count |",
            "| --- | ---: |",
        ]
    )
    for local, count in water_local.most_common(10):
        lines.append(f"| {local} | {count} |")

    lines.extend(
        [
            "",
            "## Read",
            "",
            "- The graph join is now in place, so the molecular atlas can distinguish leads that are only formula-level from leads that survive a first graph-level sanity check.",
            "- If a high-scoring pair has different graph signatures, the next interpretation is simple: the reader is probably detecting ordinary isomer graph structure.",
            "- If a high-scoring pair has the same graph signature, the next interpretation is more interesting: the separation may be conformer geometry, local curvature, solvent sensitivity, or a weakness of the heuristic graph.",
            "- The water-only refusal set should next be inspected by graph family and by the magnitude of the water projected minimum eigenvalue, before any solvent-chemistry claim is made.",
            "",
            "## Boundaries",
            "",
            "- This pass uses vacuum positions for the graph even when inspecting water support refusals.",
            "- The graph is inferred from covalent radii, not supplied as a certified chemistry graph.",
            "- No hidden-load call is licensed by this graph layer.",
        ]
    )
    if missing:
        lines.extend(["", "Missing labels:", "", ", ".join(missing)])
    (OUT_ROOT / "molecular_graph_lead_inspection_v0.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    candidates = load_candidate_pairs(TOP_PAIR_LIMIT)
    refusals = load_water_refusals(WATER_REFUSAL_LIMIT)
    labels: set[str] = set()
    for item in candidates:
        labels.add(item["left_label"])
        labels.add(item["right_label"])
    for item in refusals:
        labels.add(item["label"])

    records = collect_vacuum_rows(labels)
    signatures = {label: graph_signature(row) for label, row in records.items()}
    missing = sorted(labels - set(signatures))

    pairs = pair_rows(candidates, signatures)
    waters = water_rows(refusals, signatures)

    write_csv(OUT_ROOT / "molecular_graph_pair_inspection_v0.csv", pairs)
    write_csv(OUT_ROOT / "molecular_graph_water_refusal_examples_v0.csv", waters)
    report = {
        "schema_version": "molecular_graph_lead_inspection.v0",
        "pair_candidates_inspected": len(pairs),
        "water_refusals_inspected": len(waters),
        "missing_labels": missing,
        "same_graph_pair_count": sum(1 for row in pairs if row["graph_same_by_signature"]),
        "top_pair_rows": pairs[:10],
        "water_refusal_formula_counts": Counter(str(row["formula"]) for row in waters).most_common(25),
        "water_refusal_local_signature_counts": Counter(str(row["local_env_signature"]) for row in waters).most_common(25),
    }
    (OUT_ROOT / "molecular_graph_lead_inspection_v0.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(pairs, waters, missing)
    print(json.dumps({k: report[k] for k in ["pair_candidates_inspected", "water_refusals_inspected", "same_graph_pair_count"]}, indent=2))


if __name__ == "__main__":
    main()
