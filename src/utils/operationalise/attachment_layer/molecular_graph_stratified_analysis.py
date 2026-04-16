from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.ipc as ipc

from molecular_graph_lead_inspection import graph_signature
from molecular_isomer_geometry_analysis import distance, feature_vector
from molecular_vibrational_atlas import DATA_ROOT, OUT_ROOT


FINGERPRINT_PATH = OUT_ROOT / "molecular_isomer_fingerprints_v0.csv"
MIN_GROUP_SIZE = 2


def load_fingerprints() -> dict[str, dict[str, str]]:
    with FINGERPRINT_PATH.open(newline="", encoding="utf-8") as handle:
        return {row["label"]: row for row in csv.DictReader(handle)}


def collect_graph_signatures(labels: set[str]) -> dict[str, dict[str, object]]:
    found: dict[str, dict[str, object]] = {}
    for shard in range(5):
        path = DATA_ROOT / "vacuum" / f"data-{shard:05d}-of-00005.arrow"
        with ipc.open_stream(path) as reader:
            for batch in reader:
                slim = batch.select(["label", "atomic_numbers", "positions"])
                for row in slim.to_pylist():
                    label = str(row["label"])
                    if label in labels:
                        found[label] = graph_signature(row)
        if len(found) == len(labels):
            break
    return found


def graph_key(signature: dict[str, object]) -> str:
    return "|".join(
        [
            str(signature["atom_count_signature"]),
            str(signature["bond_signature"]),
            str(signature["local_env_signature"]),
            str(signature["heavy_degree_hist"]),
            str(signature["heavy_cyclomatic"]),
        ]
    )


def best_pairs_by_graph(
    fingerprints: dict[str, dict[str, str]], signatures: dict[str, dict[str, object]]
) -> tuple[list[dict[str, object]], dict[str, object]]:
    rich_keys = (
        "heavy_gmean_log",
        "hydrogen_gmean_log",
        "local_gmean_median_log",
        "local_gmean_log_range",
        "vacuum_softest_positive_eig_log",
        "water_soft_vacuum_log_ratio",
    )
    groups: dict[str, list[tuple[dict[str, str], tuple[float, ...]]]] = defaultdict(list)
    graph_meta: dict[str, dict[str, object]] = {}
    skipped_no_vector = 0
    skipped_no_graph = 0
    for label, row in fingerprints.items():
        sig = signatures.get(label)
        if sig is None:
            skipped_no_graph += 1
            continue
        vector = feature_vector(row, rich_keys)
        if vector is None:
            skipped_no_vector += 1
            continue
        key = graph_key(sig)
        groups[key].append((row, vector))
        graph_meta[key] = sig

    pairs: list[dict[str, object]] = []
    group_sizes = Counter()
    for key, items in groups.items():
        group_sizes[len(items)] += 1
        if len(items) < MIN_GROUP_SIZE:
            continue
        best: tuple[float, dict[str, str], dict[str, str]] | None = None
        for i in range(len(items)):
            left, left_vec = items[i]
            for j in range(i + 1, len(items)):
                right, right_vec = items[j]
                d = distance(left_vec, right_vec)
                if best is None or d > best[0]:
                    best = (d, left, right)
        if best is None:
            continue
        d, left, right = best
        meta = graph_meta[key]
        pairs.append(
            {
                "formula": left["formula"],
                "graph_group_size": len(items),
                "left_label": left["label"],
                "right_label": right["label"],
                "rich_distance": d,
                "left_pattern": left.get("support_pattern", ""),
                "right_pattern": right.get("support_pattern", ""),
                "left_heavy_gmean": left.get("heavy_gmean", ""),
                "right_heavy_gmean": right.get("heavy_gmean", ""),
                "left_local_median": left.get("local_gmean_median", ""),
                "right_local_median": right.get("local_gmean_median", ""),
                "left_water_soft_ratio": left.get("water_soft_vacuum_log_ratio", ""),
                "right_water_soft_ratio": right.get("water_soft_vacuum_log_ratio", ""),
                "bond_signature": meta["bond_signature"],
                "local_env_signature": meta["local_env_signature"],
                "heavy_degree_hist": meta["heavy_degree_hist"],
                "heavy_cyclomatic": meta["heavy_cyclomatic"],
            }
        )

    pairs.sort(key=lambda row: float(row["rich_distance"]), reverse=True)
    summary = {
        "fingerprints": len(fingerprints),
        "graph_signatures": len(signatures),
        "usable_graph_vector_rows": sum(len(items) for items in groups.values()),
        "graph_groups": len(groups),
        "multi_record_graph_groups": sum(1 for items in groups.values() if len(items) >= 2),
        "skipped_no_vector": skipped_no_vector,
        "skipped_no_graph": skipped_no_graph,
        "largest_graph_group": max((len(items) for items in groups.values()), default=0),
        "group_size_histogram": {str(k): v for k, v in sorted(group_sizes.items())},
    }
    return pairs, summary


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


def write_markdown(pairs: list[dict[str, object]], summary: dict[str, object]) -> None:
    lines = [
        "# Molecular Graph-Stratified Analysis v0",
        "",
        "Status: same-heuristic-graph search over the vacuum-admitted molecular observer fingerprints.",
        "",
        "This pass asks a sharper question than same-formula comparison: after holding fixed a simple inferred graph signature, do observer/support fingerprints still vary?",
        "",
        "The graph signature uses atom counts, inferred bond counts, heavy-atom degree histogram, local heavy-atom environments, and heavy cyclomatic count. It is a triage signature, not a certified molecular identity.",
        "",
        "## Scope",
        "",
        f"- molecule fingerprints loaded: `{summary['fingerprints']}`",
        f"- graph signatures recovered: `{summary['graph_signatures']}`",
        f"- usable graph/vector rows: `{summary['usable_graph_vector_rows']}`",
        f"- graph groups: `{summary['graph_groups']}`",
        f"- multi-record graph groups: `{summary['multi_record_graph_groups']}`",
        f"- largest graph group: `{summary['largest_graph_group']}`",
        f"- skipped for missing rich vector: `{summary['skipped_no_vector']}`",
        f"- skipped for missing graph: `{summary['skipped_no_graph']}`",
        "",
        "## Strongest Same-Graph Separations",
        "",
        "| Formula | Pair | group size | rich distance | patterns | local medians | water soft ratios | graph summary |",
        "| --- | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in pairs[:25]:
        lines.append(
            "| {formula} | {left_label} / {right_label} | {graph_group_size} | {rich_distance:.4g} | {left_pattern} / {right_pattern} | {left_local_median} / {right_local_median} | {left_water_soft_ratio} / {right_water_soft_ratio} | {bond_signature}; {heavy_degree_hist}; cyc={heavy_cyclomatic} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Read",
            "",
            "- Same-formula top pairs were mostly ordinary graph-isomer separations. This graph-stratified pass removes that easy explanation for the listed pairs under the current heuristic.",
            "- These pairs are therefore better leads for conformer geometry, local curvature differences, graph-heuristic weakness, or dataset convention effects.",
            "- The result still does not license theorem-level false-collapse language. It creates a sharper queue for manual chemistry inspection.",
            "",
            "## Boundaries",
            "",
            "- The graph signature is inferred from vacuum coordinates with a covalent-radius heuristic.",
            "- Same graph signature does not prove identical chemical graph or conformer identity.",
            "- Rich distance is an exploratory fingerprint distance, not a module invariant.",
        ]
    )
    (OUT_ROOT / "molecular_graph_stratified_analysis_v0.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    fingerprints = load_fingerprints()
    signatures = collect_graph_signatures(set(fingerprints))
    pairs, summary = best_pairs_by_graph(fingerprints, signatures)
    write_csv(OUT_ROOT / "molecular_graph_stratified_pairs_v0.csv", pairs)
    report = {
        "schema_version": "molecular_graph_stratified_analysis.v0",
        "summary": summary,
        "top_pairs": pairs[:100],
    }
    (OUT_ROOT / "molecular_graph_stratified_analysis_v0.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(pairs, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
