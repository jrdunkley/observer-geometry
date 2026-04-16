from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import run_card


def _fail(message: str) -> None:
    raise AssertionError(message)


def _read_cards(cards_dir: Path) -> list[dict[str, Any]]:
    cards: list[dict[str, Any]] = []
    for path in sorted(cards_dir.glob("*.json")):
        card = run_card.load_card(path)
        card["_path"] = str(path)
        cards.append(card)
    if not cards:
        _fail(f"no cards found in {cards_dir}")
    return cards


def _assert_nonempty_string(card: dict[str, Any], field: str) -> None:
    value = card.get(field)
    if not isinstance(value, str) or not value.strip():
        _fail(f"{card['card_id']}: {field} must be a non-empty string")


def validate_card_shape(card: dict[str, Any]) -> None:
    _assert_nonempty_string(card, "coordinate_note")
    _assert_nonempty_string(card, "positive_scope")
    _assert_nonempty_string(card, "boundary")
    _assert_nonempty_string(card, "failure_mode")
    _assert_nonempty_string(card, "next_measurement")

    sensor = card["sensor_surface"]
    if not isinstance(sensor["ideal"], list) or not sensor["ideal"]:
        _fail(f"{card['card_id']}: sensor_surface.ideal must be a non-empty list")
    if not isinstance(sensor["practical"], list) or not sensor["practical"]:
        _fail(f"{card['card_id']}: sensor_surface.practical must be a non-empty list")

    allowed = set(card["allowed_nomogeo_calls"])
    declared = card["declared_objects"]
    if "visible_precision" in allowed and not {"H", "C"} <= set(declared):
        _fail(f"{card['card_id']}: visible_precision allowed but H/C are not declared")
    if "hidden_load" in allowed and not {"H", "C", "T"} <= set(declared):
        _fail(f"{card['card_id']}: hidden_load allowed but H/C/T are not declared")
    if card["attachment_type"] == "not_yet_attachable":
        if allowed:
            _fail(f"{card['card_id']}: not_yet_attachable card must allow no calls")
        if declared:
            _fail(f"{card['card_id']}: not_yet_attachable card should not declare runnable objects in v0")


def validate_run_result(card: dict[str, Any], result: dict[str, Any]) -> None:
    allowed = set(card["allowed_nomogeo_calls"])
    executed = set(result["executed_calls"])
    if not executed <= allowed:
        _fail(f"{card['card_id']}: executed calls exceed allowed calls: {executed - allowed}")

    if card["attachment_type"] == "not_yet_attachable":
        if result["executed_calls"]:
            _fail(f"{card['card_id']}: boundary card executed calls")
        if not result["refused_calls"]:
            _fail(f"{card['card_id']}: boundary card did not record refusal")
        return

    requested = set(card["nomogeo_readout"].get("requested", []))
    if "visible_precision" in requested and "visible_precision" in allowed:
        if "visible_precision" not in result["readout"]:
            _fail(f"{card['card_id']}: missing visible_precision readout")

    if "hidden_load" in allowed:
        for field in (
            "ceiling",
            "ceiling_minus_visible_precision",
            "ceiling_dominance_min_eigenvalue",
            "ceiling_dominance_tolerance",
            "hidden_load",
            "hidden_rank",
            "clock",
        ):
            if field not in result["readout"]:
                _fail(f"{card['card_id']}: missing hidden-load readout field {field}")
        min_eig = float(result["readout"]["ceiling_dominance_min_eigenvalue"])
        tol = float(result["readout"]["ceiling_dominance_tolerance"])
        if min_eig < -tol:
            _fail(f"{card['card_id']}: ceiling dominance failed but hidden_load executed")


def validate_negative_ceiling(cards: list[dict[str, Any]]) -> None:
    hidden_cards = [card for card in cards if "hidden_load" in card["allowed_nomogeo_calls"]]
    if not hidden_cards:
        _fail("no hidden-load cards available for negative ceiling test")

    card = json.loads(json.dumps(hidden_cards[0]))
    phi_result = run_card.run_card(card)
    phi = np.asarray(phi_result["readout"]["visible_precision"], dtype=float)
    card["declared_objects"]["T"] = (phi - np.eye(phi.shape[0]) * 1.0).tolist()
    result = run_card.run_card(card)
    refused_reasons = {item["reason"] for item in result["refused_calls"]}
    if "declared ceiling does not dominate computed visible precision" not in refused_reasons:
        _fail("negative ceiling test did not refuse hidden_load")
    if "hidden_load" in result["executed_calls"]:
        _fail("negative ceiling test executed hidden_load")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate attachment-layer cards and runner invariants.")
    parser.add_argument("--nomogeo-src", help="path to observer_geometry/src")
    parser.add_argument("--cards-dir", default=str(Path(__file__).resolve().parent / "cards"))
    args = parser.parse_args()

    run_card.configure_nomogeo_import(args.nomogeo_src)
    cards = _read_cards(Path(args.cards_dir))
    for card in cards:
        validate_card_shape(card)
        result = run_card.run_card(card)
        validate_run_result(card, result)
    validate_negative_ceiling(cards)

    print(json.dumps({"validated_cards": [card["card_id"] for card in cards], "negative_ceiling_test": "passed"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
