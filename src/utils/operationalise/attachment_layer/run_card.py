from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np


REQUIRED_FIELDS = {
    "card_id",
    "schema_version",
    "system",
    "attachment_type",
    "state_variables",
    "sensor_surface",
    "control_variables",
    "regime_assumptions",
    "equation_hook",
    "units_and_dimensions",
    "coordinate_note",
    "module_object",
    "derivation_status",
    "observer_map",
    "ceiling_or_reference",
    "allowed_nomogeo_calls",
    "nomogeo_readout",
    "positive_scope",
    "boundary",
    "failure_mode",
    "next_measurement",
    "default_parameters",
    "declared_objects",
}

ATTACHMENT_TYPES = {
    "exact_quadratic",
    "exact_quadratic_with_declared_ceiling",
    "exact_special_sector",
    "local_quadratic_approximation",
    "declared_fluctuation_model",
    "weighted_family",
    "support_path",
    "not_yet_attachable",
}

IMPLEMENTED_CALLS = {"visible_precision", "hidden_load"}

visible_precision = None
hidden_load = None


def _array(value: Any, name: str) -> np.ndarray:
    if value is None:
        raise ValueError(f"declared object {name!r} is missing")
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"declared object {name!r} must be a 2D array")
    return arr


def _symmetrize(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def configure_nomogeo_import(nomogeo_src: str | None) -> None:
    """Configure and import nomogeo without assuming one workspace layout."""
    candidates: list[tuple[str, Path]] = []
    if nomogeo_src:
        candidates.append(("--nomogeo-src", Path(nomogeo_src)))
    if os.environ.get("NOMOGEO_SRC"):
        candidates.append(("NOMOGEO_SRC", Path(os.environ["NOMOGEO_SRC"])))

    for parent in Path(__file__).resolve().parents:
        if (parent / "nomogeo").exists():
            candidates.append(("package source fallback", parent))
            break
        workspace_fallback = parent / "observer_geometry" / "src"
        if workspace_fallback.exists():
            candidates.append(("workspace fallback", workspace_fallback))
            break

    for source, candidate in candidates:
        if not candidate.exists():
            raise ValueError(f"{source} path does not exist: {candidate}")
        sys.path.insert(0, str(candidate))
        try:
            _import_nomogeo()
            return
        except ModuleNotFoundError:
            try:
                sys.path.remove(str(candidate))
            except ValueError:
                pass

    try:
        _import_nomogeo()
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Could not import nomogeo. Install nomogeo, pass --nomogeo-src <path-to-observer_geometry/src>, "
            "or set NOMOGEO_SRC."
        ) from exc


def _import_nomogeo() -> None:
    global hidden_load, visible_precision
    module = importlib.import_module("nomogeo")
    visible_precision = module.visible_precision
    hidden_load = module.hidden_load


def ceiling_dominates(T: np.ndarray, phi: np.ndarray) -> tuple[bool, float, float]:
    if T.shape != phi.shape:
        return False, float("nan"), 0.0
    gap = _symmetrize(T - phi)
    eigs = np.linalg.eigvalsh(gap)
    scale = max(1.0, float(np.linalg.norm(T, ord="fro")), float(np.linalg.norm(phi, ord="fro")))
    tol = 1.0e-10 * scale
    return bool(float(np.min(eigs)) >= -tol), float(np.min(eigs)), tol


def load_card(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        card = json.load(handle)
    missing = sorted(REQUIRED_FIELDS - set(card))
    if missing:
        raise ValueError(f"{path.name} missing required fields: {missing}")
    if card["schema_version"] != "attachment_card.v0":
        raise ValueError(f"{path.name} has unsupported schema_version {card['schema_version']!r}")
    if card["attachment_type"] not in ATTACHMENT_TYPES:
        raise ValueError(f"{path.name} has unsupported attachment_type {card['attachment_type']!r}")
    if "ideal" not in card["sensor_surface"] or "practical" not in card["sensor_surface"]:
        raise ValueError(f"{path.name} sensor_surface must include ideal and practical")
    if not isinstance(card["control_variables"], list) or not card["control_variables"]:
        raise ValueError(f"{path.name} control_variables must be a non-empty list")
    allowed = set(card["allowed_nomogeo_calls"])
    unknown_allowed = sorted(allowed - IMPLEMENTED_CALLS)
    if unknown_allowed:
        raise ValueError(f"{path.name} allows calls not implemented by v0 runner: {unknown_allowed}")
    if card["attachment_type"] == "not_yet_attachable" and allowed:
        raise ValueError(f"{path.name} is not_yet_attachable but allows module calls")
    return card


def run_card(card: dict[str, Any]) -> dict[str, Any]:
    requested = list(card["nomogeo_readout"].get("requested", []))
    allowed = set(card["allowed_nomogeo_calls"])
    declared = card["declared_objects"]
    result: dict[str, Any] = {
        "card_id": card["card_id"],
        "system": card["system"],
        "attachment_type": card["attachment_type"],
        "executed_calls": [],
        "refused_calls": [],
        "readout": {},
        "coordinate_note": card["coordinate_note"],
        "positive_scope": card["positive_scope"],
        "boundary": card["boundary"],
    }

    if card["attachment_type"] == "not_yet_attachable":
        result["refused_calls"].append(
            {
                "call": "*",
                "reason": "card is not_yet_attachable; no theorem-local kernel calls are licensed",
            }
        )
        return result

    phi: np.ndarray | None = None
    if "visible_precision" in requested or "hidden_load" in requested or "clock" in requested:
        if "visible_precision" not in allowed:
            result["refused_calls"].append({"call": "visible_precision", "reason": "not in allowed_nomogeo_calls"})
        else:
            H = _array(declared.get("H"), "H")
            C = _array(declared.get("C"), "C")
            if visible_precision is None:
                raise RuntimeError("nomogeo visible_precision is not configured")
            phi = visible_precision(H, C)
            result["readout"]["visible_precision"] = phi.tolist()
            result["executed_calls"].append("visible_precision")

    if "hidden_load" in requested or "clock" in requested:
        if "hidden_load" not in allowed:
            result["refused_calls"].append({"call": "hidden_load", "reason": "not in allowed_nomogeo_calls"})
        elif phi is None:
            result["refused_calls"].append({"call": "hidden_load", "reason": "visible precision was not computed"})
        elif "T" not in declared:
            result["refused_calls"].append({"call": "hidden_load", "reason": "declared ceiling/reference T is missing"})
        else:
            T = _array(declared.get("T"), "T")
            dominates, min_gap_eig, gap_tol = ceiling_dominates(T, phi)
            result["readout"]["ceiling"] = T.tolist()
            result["readout"]["ceiling_minus_visible_precision"] = _symmetrize(T - phi).tolist()
            result["readout"]["ceiling_dominance_min_eigenvalue"] = min_gap_eig
            result["readout"]["ceiling_dominance_tolerance"] = gap_tol
            if not dominates:
                result["refused_calls"].append(
                    {
                        "call": "hidden_load",
                        "reason": "declared ceiling does not dominate computed visible precision",
                    }
                )
                return result
            if hidden_load is None:
                raise RuntimeError("nomogeo hidden_load is not configured")
            load = hidden_load(T, phi)
            result["readout"]["hidden_load"] = load.lambda_.tolist()
            result["readout"]["hidden_rank"] = load.rank
            result["readout"]["clock"] = load.clock
            result["executed_calls"].append("hidden_load")

    return result


def write_outputs(card: dict[str, Any], result: dict[str, Any], output_dir: Path, markdown: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{card['card_id']}.json"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    if markdown:
        md_path = output_dir / f"{card['card_id']}.md"
        with md_path.open("w", encoding="utf-8") as handle:
            handle.write(render_markdown(card, result))


def render_markdown(card: dict[str, Any], result: dict[str, Any]) -> str:
    lines = [
        f"# {card['system']}",
        "",
        f"- card: `{card['card_id']}`",
        f"- attachment type: `{card['attachment_type']}`",
        f"- derivation status: {card['derivation_status']}",
        f"- equation hook: `{card['equation_hook']}`",
        f"- coordinate note: {card['coordinate_note']}",
        f"- allowed calls: {', '.join(card['allowed_nomogeo_calls']) or 'none'}",
        "",
        "## Scope",
        "",
        card["positive_scope"],
        "",
        card["boundary"],
        "",
        "## Sensor Surface",
        "",
        f"- ideal: {', '.join(card['sensor_surface']['ideal'])}",
        f"- practical: {', '.join(card['sensor_surface']['practical'])}",
        "",
        "## Readout",
        "",
        "```json",
        json.dumps(result["readout"], indent=2),
        "```",
        "",
        "## Refusals",
        "",
        "```json",
        json.dumps(result["refused_calls"], indent=2),
        "```",
        "",
    ]
    return "\n".join(lines)


def card_paths(args: argparse.Namespace) -> list[Path]:
    cards_dir = Path(__file__).resolve().parent / "cards"
    if args.all:
        return sorted(cards_dir.glob("*.json"))
    if args.card:
        candidate = Path(args.card)
        if candidate.exists():
            return [candidate]
        by_id = cards_dir / f"{args.card}.json"
        if by_id.exists():
            return [by_id]
    raise SystemExit("Provide --all or --card <card_id|path>")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run operational attachment cards.")
    parser.add_argument("--all", action="store_true", help="run every card in cards/")
    parser.add_argument("--card", help="card id or JSON path")
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "outputs"))
    parser.add_argument("--markdown", action="store_true", help="also emit generated Markdown summaries")
    parser.add_argument("--nomogeo-src", help="path to observer_geometry/src; overrides installed-package lookup")
    args = parser.parse_args()
    configure_nomogeo_import(args.nomogeo_src)

    overall: list[dict[str, Any]] = []
    output_dir = Path(args.output_dir)
    for path in card_paths(args):
        card = load_card(path)
        result = run_card(card)
        write_outputs(card, result, output_dir, args.markdown)
        overall.append(
            {
                "card_id": card["card_id"],
                "attachment_type": card["attachment_type"],
                "executed_calls": result["executed_calls"],
                "refused_calls": result["refused_calls"],
            }
        )

    print(json.dumps({"ran": len(overall), "cards": overall}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
