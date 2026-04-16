from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

import run_card


def _load_card(card_id: str) -> dict[str, Any]:
    path = Path(__file__).resolve().parent / "cards" / f"{card_id}.json"
    return run_card.load_card(path)


def _scalar(readout: dict[str, Any], field: str) -> float:
    value = readout[field]
    if isinstance(value, list):
        return float(np.asarray(value, dtype=float).reshape(-1)[0])
    return float(value)


def _spring_card(base: dict[str, Any], kc: float) -> dict[str, Any]:
    card = copy.deepcopy(base)
    params = dict(card["default_parameters"])
    params["kc"] = float(kc)
    k1 = float(params["k1"])
    k2 = float(params["k2"])
    card["default_parameters"] = params
    card["declared_objects"]["H"] = [[k1 + kc, -kc], [-kc, k2 + kc]]
    card["declared_objects"]["T"] = [[k1 + kc]]
    return card


def _lc_card(base: dict[str, Any], coupling_reciprocal: float) -> dict[str, Any]:
    card = copy.deepcopy(base)
    params = dict(card["default_parameters"])
    params["Lc"] = float(1.0 / coupling_reciprocal)
    L1 = float(params["L1"])
    L2 = float(params["L2"])
    g1 = 1.0 / L1
    g2 = 1.0 / L2
    gc = float(coupling_reciprocal)
    card["default_parameters"] = params
    card["declared_objects"]["H"] = [[g1 + gc, -gc], [-gc, g2 + gc]]
    card["declared_objects"]["T"] = [[g1 + gc]]
    return card


def _row(label: str, coupling: float, card: dict[str, Any]) -> dict[str, Any]:
    result = run_card.run_card(card)
    if result["refused_calls"]:
        raise RuntimeError(f"{label} coupling {coupling} refused calls: {result['refused_calls']}")
    readout = result["readout"]
    phi = _scalar(readout, "visible_precision")
    ceiling = _scalar(readout, "ceiling")
    gap = _scalar(readout, "ceiling_minus_visible_precision")
    hidden = _scalar(readout, "hidden_load")
    clock = float(readout["clock"])
    return {
        "family": label,
        "coupling": coupling,
        "visible_precision": phi,
        "ceiling": ceiling,
        "ceiling_minus_visible_precision": gap,
        "hidden_load": hidden,
        "clock": clock,
        "ceiling_dominance_min_eigenvalue": float(readout["ceiling_dominance_min_eigenvalue"]),
        "executed_calls": result["executed_calls"],
    }


def _assert_monotone(rows: list[dict[str, Any]], field: str) -> None:
    values = [float(row[field]) for row in rows]
    if any(b <= a for a, b in zip(values, values[1:])):
        raise AssertionError(f"{field} is not strictly increasing: {values}")


def _render_markdown(rows: list[dict[str, Any]]) -> str:
    lines = [
        "# Coupled Hidden-Renormalisation Analogue Sweep",
        "",
        "This generated check varies coupling in the coupled spring and coupled LC cards while keeping the same one-visible, one-hidden geometry.",
        "",
        "Expected invariant pattern: stronger hidden coupling increases the visible-block ceiling, the relaxation gap `T - Phi`, the hidden load, and the determinant clock.",
        "",
        "| family | coupling | Phi | T | T - Phi | hidden load | clock |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {family} | {coupling:.6g} | {visible_precision:.12g} | {ceiling:.12g} | "
            "{ceiling_minus_visible_precision:.12g} | {hidden_load:.12g} | {clock:.12g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "Interpretation: both substrates realise the same qualitative observer-geometry story. The measured coordinate sees a relaxed visible precision below the visible block ceiling, and the gap grows as hidden coupling is strengthened.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a small coupled spring / coupled LC analogue sweep.")
    parser.add_argument("--nomogeo-src", help="path to observer_geometry/src")
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "outputs" / "sweeps"))
    args = parser.parse_args()

    run_card.configure_nomogeo_import(args.nomogeo_src)
    spring = _load_card("mass_spring_coupled_observe_one")
    lc = _load_card("coupled_lc_resonators_observe_one_node")

    spring_couplings = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    lc_couplings = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0]

    spring_rows = [_row("coupled_spring", kc, _spring_card(spring, kc)) for kc in spring_couplings]
    lc_rows = [_row("coupled_lc", gc, _lc_card(lc, gc)) for gc in lc_couplings]

    for rows in (spring_rows, lc_rows):
        _assert_monotone(rows, "ceiling")
        _assert_monotone(rows, "ceiling_minus_visible_precision")
        _assert_monotone(rows, "hidden_load")
        _assert_monotone(rows, "clock")

    all_rows = spring_rows + lc_rows
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "coupled_hidden_renormalisation_sweep.json"
    csv_path = output_dir / "coupled_hidden_renormalisation_sweep.csv"
    md_path = output_dir / "coupled_hidden_renormalisation_sweep.md"

    payload = {
        "schema_version": "attachment_sweep.v0",
        "description": "Coupled spring and coupled LC hidden-renormalisation sweep.",
        "monotone_checks": ["ceiling", "ceiling_minus_visible_precision", "hidden_load", "clock"],
        "rows": all_rows,
    }
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "family",
                "coupling",
                "visible_precision",
                "ceiling",
                "ceiling_minus_visible_precision",
                "hidden_load",
                "clock",
                "ceiling_dominance_min_eigenvalue",
            ],
        )
        writer.writeheader()
        for row in all_rows:
            writer.writerow({key: row[key] for key in writer.fieldnames})
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write(_render_markdown(all_rows))

    print(json.dumps({"rows": len(all_rows), "outputs": [str(json_path), str(csv_path), str(md_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
