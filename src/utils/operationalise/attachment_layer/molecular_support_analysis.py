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


def load_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for split in SPLITS:
        for shard in range(5):
            path = OUT_ROOT / f"hessian_qm9_{split}_shard{shard}_support_full_support_rows.csv"
            with path.open(newline="", encoding="utf-8") as handle:
                rows.extend(csv.DictReader(handle))
    return rows


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def quantiles(values: list[float]) -> dict[str, float]:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return {}
    return {
        "min": clean[0],
        "q01": clean[int(0.01 * (len(clean) - 1))],
        "q10": clean[int(0.10 * (len(clean) - 1))],
        "median": stats.median(clean),
        "q90": clean[int(0.90 * (len(clean) - 1))],
        "q99": clean[int(0.99 * (len(clean) - 1))],
        "max": clean[-1],
    }


def per_split(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for split in SPLITS:
        split_rows = [row for row in rows if row["split"] == split]
        spd = [row for row in split_rows if row["support_status"] == "vibrational_spd"]
        indefinite = [row for row in split_rows if row["support_status"] != "vibrational_spd"]
        out[split] = {
            "rows": len(split_rows),
            "spd": len(spd),
            "indefinite_or_singular": len(indefinite),
            "spd_rate": len(spd) / len(split_rows),
            "vib_condition": quantiles([f(row, "vib_condition") for row in split_rows]),
            "softest_positive_eig": quantiles([f(row, "softest_positive_eig") for row in split_rows]),
            "soft_mode_hydrogen_share": quantiles([f(row, "soft_mode_hydrogen_share") for row in split_rows]),
            "soft_mode_hetero_share": quantiles([f(row, "soft_mode_hetero_share") for row in split_rows]),
            "frequency_log_spectrum_corr": quantiles(
                [f(row, "frequency_log_spectrum_corr") for row in split_rows if row["frequency_log_spectrum_corr"]]
            ),
            "soft_mode_top_element": dict(Counter(row["soft_mode_top_element"] for row in split_rows)),
            "rigid_geometry": dict(Counter(row["rigid_geometry"] for row in split_rows)),
        }
    return out


def status_patterns(rows: list[dict[str, str]]) -> dict[str, int]:
    by_label: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by_label[row["label"]][row["split"]] = row
    counts: Counter[str] = Counter()
    for split_rows in by_label.values():
        if all(split in split_rows for split in SPLITS):
            pattern = "".join("S" if split_rows[split]["support_status"] == "vibrational_spd" else "I" for split in SPLITS)
            counts[pattern] += 1
    return dict(counts.most_common())


def formula_spreads(rows: list[dict[str, str]], min_count: int = 20) -> list[dict[str, object]]:
    formula_split: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for row in rows:
        bucket = formula_split[row["formula"]][row["split"]]
        bucket[1] += 1
        if row["support_status"] == "vibrational_spd":
            bucket[0] += 1
    out = []
    for formula, by_split in formula_split.items():
        if not all(by_split[split][1] >= min_count for split in SPLITS):
            continue
        rates = {split: by_split[split][0] / by_split[split][1] for split in SPLITS}
        counts = {split: by_split[split][1] for split in SPLITS}
        out.append({"formula": formula, "spread": max(rates.values()) - min(rates.values()), "rates": rates, "counts": counts})
    return sorted(out, key=lambda item: float(item["spread"]), reverse=True)


def environment_soft_deltas(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    by_label: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by_label[row["label"]][row["split"]] = row
    out: dict[str, dict[str, object]] = {}
    for target in ("thf", "toluene", "water"):
        values: list[tuple[float, str, str, str, str]] = []
        for label, split_rows in by_label.items():
            if not all(split in split_rows for split in SPLITS):
                continue
            vac = f(split_rows["vacuum"], "softest_positive_eig")
            other = f(split_rows[target], "softest_positive_eig")
            if vac <= 0.0 or other <= 0.0:
                continue
            values.append(
                (
                    math.log(other / vac),
                    label,
                    split_rows["vacuum"]["formula"],
                    split_rows["vacuum"]["support_status"],
                    split_rows[target]["support_status"],
                )
            )
        sorted_values = sorted(v[0] for v in values)
        out[target] = {
            "n": len(values),
            "log_ratio": quantiles(sorted_values),
            "lowest": values and [dict(log_ratio=v[0], label=v[1], formula=v[2], vacuum_status=v[3], target_status=v[4]) for v in sorted(values)[:10]],
            "highest": values and [dict(log_ratio=v[0], label=v[1], formula=v[2], vacuum_status=v[3], target_status=v[4]) for v in sorted(values)[-10:]],
        }
    return out


def mode_participation(rows: list[dict[str, str]]) -> dict[str, object]:
    by_split: dict[str, dict[str, object]] = {}
    for split in SPLITS:
        split_rows = [row for row in rows if row["split"] == split]
        by_split[split] = {
            "top_element_counts": dict(Counter(row["soft_mode_top_element"] for row in split_rows)),
            "hydrogen_share": quantiles([f(row, "soft_mode_hydrogen_share") for row in split_rows]),
            "hetero_share": quantiles([f(row, "soft_mode_hetero_share") for row in split_rows]),
        }
    return by_split


def write_pattern_csv(patterns: dict[str, int]) -> None:
    path = OUT_ROOT / "molecular_support_status_patterns_v0.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pattern", "count", "meaning"])
        writer.writeheader()
        for pattern, count in patterns.items():
            meaning = ", ".join(f"{split}:{state}" for split, state in zip(SPLITS, pattern, strict=True))
            writer.writerow({"pattern": pattern, "count": count, "meaning": meaning})


def write_formula_csv(spreads: list[dict[str, object]]) -> None:
    path = OUT_ROOT / "molecular_support_formula_spreads_v0.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["formula", "spread"] + [f"{split}_rate" for split in SPLITS] + [f"{split}_count" for split in SPLITS]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in spreads:
            row = {"formula": item["formula"], "spread": item["spread"]}
            row.update({f"{split}_rate": item["rates"][split] for split in SPLITS})
            row.update({f"{split}_count": item["counts"][split] for split in SPLITS})
            writer.writerow(row)


def write_label_pattern_csv(rows: list[dict[str, str]]) -> None:
    by_label: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by_label[row["label"]][row["split"]] = row
    path = OUT_ROOT / "molecular_support_label_patterns_v0.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["label", "formula", "pattern"] + [f"{split}_softest_positive_eig" for split in SPLITS] + [
            f"{split}_vib_min_eig" for split in SPLITS
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for label in sorted(by_label):
            split_rows = by_label[label]
            if not all(split in split_rows for split in SPLITS):
                continue
            pattern = "".join("S" if split_rows[split]["support_status"] == "vibrational_spd" else "I" for split in SPLITS)
            row = {"label": label, "formula": split_rows["vacuum"]["formula"], "pattern": pattern}
            row.update({f"{split}_softest_positive_eig": split_rows[split]["softest_positive_eig"] for split in SPLITS})
            row.update({f"{split}_vib_min_eig": split_rows[split]["vib_min_eig"] for split in SPLITS})
            writer.writerow(row)


def write_environment_extremes_csv(report: dict[str, object]) -> None:
    path = OUT_ROOT / "molecular_support_environment_soft_extremes_v0.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["environment", "tail", "rank", "log_ratio", "label", "formula", "vacuum_status", "target_status"])
        writer.writeheader()
        for environment, payload in report["environment_soft_deltas"].items():
            for tail in ("lowest", "highest"):
                for rank, item in enumerate(payload[tail], start=1):
                    writer.writerow(
                        {
                            "environment": environment,
                            "tail": tail,
                            "rank": rank,
                            "log_ratio": item["log_ratio"],
                            "label": item["label"],
                            "formula": item["formula"],
                            "vacuum_status": item["vacuum_status"],
                            "target_status": item["target_status"],
                        }
                    )


def md_table_split(split_stats: dict[str, dict[str, object]]) -> str:
    lines = [
        "| Split | Rows | SPD | Indefinite/singular | SPD rate | soft eig median | condition median | soft H share median | top soft element |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for split in SPLITS:
        item = split_stats[split]
        top_element = max(item["soft_mode_top_element"].items(), key=lambda kv: kv[1])[0]
        lines.append(
            "| {split} | {rows} | {spd} | {indef} | {rate:.4f} | {soft:.6g} | {cond:.3f} | {hshare:.3f} | {top} |".format(
                split=split,
                rows=item["rows"],
                spd=item["spd"],
                indef=item["indefinite_or_singular"],
                rate=item["spd_rate"],
                soft=item["softest_positive_eig"]["median"],
                cond=item["vib_condition"]["median"],
                hshare=item["soft_mode_hydrogen_share"]["median"],
                top=top_element,
            )
        )
    return "\n".join(lines)


def write_markdown(report: dict[str, object]) -> None:
    lines = [
        "# Molecular Support Atlas v0 Analysis",
        "",
        "Status: generated aggregate analysis over all four Hessian QM9 environments.",
        "",
        "This is a support/provenance layer, not a hidden-load or visible-precision layer. It uses the information removed from the first SPD reader: rigid nullspace, vibrational support, signed projected inertia, soft modes, and spectrum validation.",
        "",
        "## Split Summary",
        "",
        md_table_split(report["per_split"]),
        "",
        "## Cross-Environment Support Patterns",
        "",
        "Pattern order is `vacuum`, `thf`, `toluene`, `water`; `S` means projected vibrational SPD and `I` means indefinite or singular under the declared tolerance.",
        "",
        "| Pattern | Count |",
        "| --- | ---: |",
    ]
    for pattern, count in list(report["status_patterns"].items())[:16]:
        lines.append(f"| {pattern} | {count} |")

    lines.extend(
        [
            "",
            "## Formula Families With Large Environment Spread",
            "",
            "These are formula-level support admission rates, not isomer-resolved chemistry. They identify where solvent/environment changes most strongly alter the projected support gate.",
            "",
            "| Formula | Spread | Vacuum | THF | Toluene | Water |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for item in report["formula_spreads"][:20]:
        rates = item["rates"]
        lines.append(
            f"| {item['formula']} | {item['spread']:.3f} | {rates['vacuum']:.3f} | {rates['thf']:.3f} | {rates['toluene']:.3f} | {rates['water']:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Soft-Mode Environment Drift",
            "",
            "Log ratio is `log(softest_positive_eig_environment / softest_positive_eig_vacuum)`. Negative medians mean the environment tends to soften the lowest positive mode under this compiler convention.",
            "",
            "| Environment | n | q10 | median | q90 |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for target in ("thf", "toluene", "water"):
        q = report["environment_soft_deltas"][target]["log_ratio"]
        lines.append(f"| {target} | {report['environment_soft_deltas'][target]['n']} | {q['q10']:.3f} | {q['median']:.3f} | {q['q90']:.3f} |")

    lines.extend(
        [
            "",
            "## Main Signals",
            "",
            "- Toluene is the cleanest projected-support environment in this pass: SPD rate `0.9909` versus vacuum `0.9745`, THF `0.9699`, and water `0.9560`.",
            "- Water has the strongest softening and the highest indefinite/singular rate. This is a compelling lead, but it is not yet a solvent-physics claim; it may include optimization, finite-difference, tolerance, or chart effects.",
            "- The softest-mode top element shifts by environment. Vacuum is mostly carbon-top (`20,127`) with oxygen second (`16,011`); THF, toluene, and water become oxygen-top more often.",
            "- Softest modes become more hydrogen-involved in solvents: median hydrogen participation rises from `0.179` in vacuum to `0.215` in THF, `0.207` in toluene, and `0.222` in water.",
            "- Most molecules are stable across all four environments (`SSSS`: `37,799`), but the neighboring refusal patterns are informative. `SSSI` (`1,446`) is especially important: clean in vacuum/THF/toluene, refused in water.",
            "",
            "## Boundaries",
            "",
            "- This analysis does not call `hidden_load`; no molecular ceiling has been declared.",
            "- The support event language here is static across separate optimized Hessians, not a continuous reaction or solvent path.",
            "- Formula-level spreads mix many isomers and must not be read as a functional-group theorem.",
            "- Rigid geometry is reported from the mass-weighted rigid-motion basis; all records in this pass appear nonlinear under this numerical criterion.",
            "",
            "## Next Scientific Moves",
            "",
            "1. Inspect `SSSI` and `SISS` labels: water-only and THF-only support failures are the cleanest environment-specific leads.",
            "2. Add isomer grouping: same formula, different observer/support fingerprints, then test false collapse under coarse formula-level views.",
            "3. Validate solvent softening against dataset frequencies directly, not only projected Hessian eigenvalues.",
            "4. Add a careful molecular graph layer so local-site observers can be grouped by atom environment rather than only atom type and attached-H count.",
            "5. Decide a molecular ceiling convention only after observer charts and graph groups are stable.",
        ]
    )
    (OUT_ROOT / "molecular_support_analysis_v0.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows = load_rows()
    report = {
        "schema_version": "molecular_support_analysis.v0",
        "rows": len(rows),
        "per_split": per_split(rows),
        "status_patterns": status_patterns(rows),
        "formula_spreads": formula_spreads(rows),
        "environment_soft_deltas": environment_soft_deltas(rows),
        "mode_participation": mode_participation(rows),
    }
    write_pattern_csv(report["status_patterns"])
    write_formula_csv(report["formula_spreads"])
    write_label_pattern_csv(rows)
    write_environment_extremes_csv(report)
    (OUT_ROOT / "molecular_support_analysis_v0.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(report)
    print(json.dumps({"rows": report["rows"], "patterns": len(report["status_patterns"]), "formula_spreads": len(report["formula_spreads"])}, indent=2))


if __name__ == "__main__":
    main()
