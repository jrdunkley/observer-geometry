from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .engine import fetch_coinbase_hourly, run_observer_gate, synthetic_hourly_candles


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep fixed BTC observer-gate settings")
    parser.add_argument("--live", action="store_true", help="fetch public Coinbase BTC-USD hourly candles once")
    parser.add_argument("--limit", type=int, default=300, help="candle count")
    args = parser.parse_args()

    candles = fetch_coinbase_hourly(limit=args.limit) if args.live else synthetic_hourly_candles(args.limit)
    rows = []
    for train_window in [48, 72, 96, 120, 144]:
        for abstain_z in [0.0, 0.15, 0.3, 0.5, 0.8]:
            try:
                result = run_observer_gate(
                    candles,
                    train_window=train_window,
                    abstain_z=abstain_z,
                    source="coinbase" if args.live else "synthetic",
                )
                metrics = result["metrics"]
                rows.append(
                    {
                        "mode": "live_coinbase" if args.live else "synthetic_offline",
                        "first_candle": result["first_candle"],
                        "last_candle": result["last_candle"],
                        "last_close": result["last_close"],
                        "train_window": train_window,
                        "abstain_z": abstain_z,
                        "test_points": metrics["test_points"],
                        "selected_rmse": metrics["selected_rmse"],
                        "zero_return_rmse": metrics["zero_return_rmse"],
                        "all_feature_rmse": metrics["all_feature_rmse"],
                        "selected_direction_accuracy_when_speaking": metrics[
                            "selected_direction_accuracy_when_speaking"
                        ],
                        "last_return_direction_accuracy": metrics["last_return_direction_accuracy"],
                        "momentum_6_direction_accuracy": metrics["momentum_6_direction_accuracy"],
                        "selected_coverage": metrics["selected_coverage"],
                        "gate_passes": metrics["gate_passes"],
                        "gate_reason": metrics["gate_reason"],
                        "observer_wins": json.dumps(metrics["observer_wins"], sort_keys=True),
                    }
                )
            except Exception as exc:
                rows.append(
                    {
                        "mode": "live_coinbase" if args.live else "synthetic_offline",
                        "first_candle": "",
                        "last_candle": "",
                        "last_close": "",
                        "train_window": train_window,
                        "abstain_z": abstain_z,
                        "test_points": "",
                        "selected_rmse": "",
                        "zero_return_rmse": "",
                        "all_feature_rmse": "",
                        "selected_direction_accuracy_when_speaking": "",
                        "last_return_direction_accuracy": "",
                        "momentum_6_direction_accuracy": "",
                        "selected_coverage": "",
                        "gate_passes": False,
                        "gate_reason": f"{type(exc).__name__}: {exc}",
                        "observer_wins": "{}",
                    }
                )

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / ("live_sweep.csv" if args.live else "synthetic_sweep.csv")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    passed = [row for row in rows if row["gate_passes"] is True]
    best_rmse = min(
        (row for row in rows if row["selected_rmse"] != ""),
        key=lambda row: float(row["selected_rmse"]) / float(row["zero_return_rmse"]),
    )
    best_hit = max(
        (row for row in rows if row["selected_direction_accuracy_when_speaking"] not in {"", None}),
        key=lambda row: float(row["selected_direction_accuracy_when_speaking"] or 0.0),
    )

    print(f"sweep rows: {len(rows)}")
    print(f"gate passes: {len(passed)}")
    print(
        "best RMSE ratio selected/zero: {:.4f} at train_window={}, abstain_z={}".format(
            float(best_rmse["selected_rmse"]) / float(best_rmse["zero_return_rmse"]),
            best_rmse["train_window"],
            best_rmse["abstain_z"],
        )
    )
    print(
        "best hit rate when speaking: {} at train_window={}, abstain_z={}, coverage={}".format(
            best_hit["selected_direction_accuracy_when_speaking"],
            best_hit["train_window"],
            best_hit["abstain_z"],
            best_hit["selected_coverage"],
        )
    )
    print(f"output: {path}")


if __name__ == "__main__":
    main()
