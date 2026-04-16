from __future__ import annotations

import argparse
from pathlib import Path

from .engine import fetch_coinbase_hourly, run_observer_gate, synthetic_hourly_candles, write_outputs


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    parser = argparse.ArgumentParser(description="BTC observer-gate research harness")
    parser.add_argument("--live", action="store_true", help="fetch public Coinbase BTC-USD hourly candles")
    parser.add_argument("--limit", type=int, default=300, help="Coinbase candle limit, 50..300")
    parser.add_argument("--train-window", type=int, default=96, help="rolling train window in hours")
    parser.add_argument("--abstain-z", type=float, default=0.15, help="fixed speak/abstain z threshold")
    args = parser.parse_args()

    candles = fetch_coinbase_hourly(limit=args.limit) if args.live else synthetic_hourly_candles(args.limit)
    result = run_observer_gate(
        candles,
        train_window=args.train_window,
        abstain_z=args.abstain_z,
        source="coinbase" if args.live else "synthetic",
    )
    result["mode"] = "live_coinbase" if args.live else "synthetic_offline"
    write_outputs(result, OUT)
    print("BTC observer gate complete.")
    print(f"mode: {result['mode']}")
    print(f"candles: {result['candle_count']} ({result['first_candle']} to {result['last_candle']})")
    print(f"last close: {result['last_close']:.2f}")
    print(f"test points: {result['metrics']['test_points']}")
    print(f"selected RMSE: {result['metrics']['selected_rmse']:.6g}")
    print(f"zero-return RMSE: {result['metrics']['zero_return_rmse']:.6g}")
    print(f"coverage: {100.0 * result['metrics']['selected_coverage']:.1f}%")
    print(f"gate passes: {result['metrics']['gate_passes']} ({result['metrics']['gate_reason']})")
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
