from __future__ import annotations

import json
import math
import statistics
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from nomogeo import Tolerances, hidden_load, visible_precision
from nomogeo.exceptions import NomogeoError


FEATURES = (
    "ret_1",
    "mom_3",
    "mom_6",
    "mom_12",
    "mom_24",
    "vol_6",
    "vol_24",
    "range_1",
    "body_1",
    "volume_z_24",
)

OBSERVERS: dict[str, tuple[str, ...]] = {
    "impulse": ("ret_1", "mom_3", "mom_6"),
    "trend_context": ("mom_3", "mom_12", "mom_24"),
    "volatility_state": ("vol_6", "vol_24", "range_1"),
    "flow_pressure": ("volume_z_24", "body_1", "ret_1"),
    "mixed_local": ("ret_1", "mom_6", "vol_6", "volume_z_24", "range_1"),
}


@dataclass(frozen=True)
class Candle:
    time: int
    low: float
    high: float
    open: float
    close: float
    volume: float


def utc_iso(ts: int) -> str:
    return datetime.fromtimestamp(int(ts), tz=UTC).isoformat().replace("+00:00", "Z")


def fetch_coinbase_hourly(product: str = "BTC-USD", limit: int = 300, timeout: float = 12.0) -> list[Candle]:
    """Fetch recent hourly candles from the public Coinbase Exchange API."""
    if limit < 50 or limit > 300:
        raise ValueError("Coinbase candle limit must be in [50, 300] for this harness")
    query = urllib.parse.urlencode({"granularity": 3600, "limit": int(limit)})
    url = f"https://api.exchange.coinbase.com/products/{urllib.parse.quote(product)}/candles?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "nomogeo-btc-observer-gate/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode("utf-8"))
    candles = [Candle(*(float(x) if j else int(x) for j, x in enumerate(row))) for row in payload]
    return sorted(candles, key=lambda c: c.time)[-limit:]


def synthetic_hourly_candles(n: int = 300, seed: int = 20260411) -> list[Candle]:
    """Deterministic regime-switching candles for offline validation."""
    rng = np.random.default_rng(seed)
    price = 70000.0
    candles: list[Candle] = []
    ts0 = 1773600000
    vol_state = 0.006
    for i in range(n):
        if i in {90, 180, 240}:
            vol_state *= 1.7 if i != 180 else 0.45
        drift = 0.00035 * math.sin(i / 17.0) - 0.00025 * math.cos(i / 41.0)
        shock = rng.normal(0.0, vol_state)
        ret = drift + shock
        open_ = price
        close = max(1000.0, open_ * math.exp(ret))
        wick = abs(rng.normal(0.0, vol_state * 0.75))
        high = max(open_, close) * math.exp(wick)
        low = min(open_, close) * math.exp(-wick)
        volume = 150.0 + 80.0 * abs(shock) / max(vol_state, 1e-12) + rng.gamma(2.0, 25.0)
        candles.append(Candle(ts0 + 3600 * i, low, high, open_, close, volume))
        price = close
    return candles


def build_rows(candles: list[Candle]) -> list[dict[str, Any]]:
    if len(candles) < 80:
        raise ValueError("need at least 80 hourly candles")
    closes = np.array([c.close for c in candles], dtype=float)
    opens = np.array([c.open for c in candles], dtype=float)
    highs = np.array([c.high for c in candles], dtype=float)
    lows = np.array([c.low for c in candles], dtype=float)
    volumes = np.array([c.volume for c in candles], dtype=float)
    log_close = np.log(closes)
    returns = np.zeros_like(log_close)
    returns[1:] = np.diff(log_close)
    log_volume = np.log1p(volumes)

    rows: list[dict[str, Any]] = []
    for t in range(25, len(candles) - 1):
        vol_slice = log_volume[t - 23 : t + 1]
        vol_std = float(np.std(vol_slice))
        if vol_std <= 1e-12:
            vol_std = 1.0
        features = {
            "ret_1": float(returns[t]),
            "mom_3": float(np.sum(returns[t - 2 : t + 1])),
            "mom_6": float(np.sum(returns[t - 5 : t + 1])),
            "mom_12": float(np.sum(returns[t - 11 : t + 1])),
            "mom_24": float(np.sum(returns[t - 23 : t + 1])),
            "vol_6": float(np.std(returns[t - 5 : t + 1])),
            "vol_24": float(np.std(returns[t - 23 : t + 1])),
            "range_1": float(math.log(highs[t] / lows[t])),
            "body_1": float(math.log(closes[t] / opens[t])),
            "volume_z_24": float((log_volume[t] - float(np.mean(vol_slice))) / vol_std),
        }
        rows.append(
            {
                "time": candles[t].time,
                "time_iso": utc_iso(candles[t].time),
                "close": float(closes[t]),
                "target_time": candles[t + 1].time,
                "target_time_iso": utc_iso(candles[t + 1].time),
                "target_return": float(returns[t + 1]),
                "features": features,
            }
        )
    return rows


def _standardize(train_x: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(train_x, axis=0)
    scale = np.std(train_x, axis=0)
    scale = np.where(scale <= 1e-12, 1.0, scale)
    return (train_x - mean) / scale, (x - mean) / scale, scale


def _ridge_predict(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, ridge: float = 1.0e-3) -> float:
    y_mean = float(np.mean(train_y))
    centered_y = train_y - y_mean
    xtx = train_x.T @ train_x
    beta = np.linalg.solve(xtx + ridge * np.eye(train_x.shape[1]), train_x.T @ centered_y)
    return float(y_mean + test_x @ beta)


def _geometry(train_x_all: np.ndarray, train_y: np.ndarray, observer: tuple[str, ...]) -> dict[str, float | str]:
    xz, _unused, _scale = _standardize(train_x_all, train_x_all[0])
    y_mean = float(np.mean(train_y))
    y_std = float(np.std(train_y))
    if y_std <= 1e-12:
        y_std = 1.0
    yz = ((train_y - y_mean) / y_std).reshape(-1, 1)
    z = np.hstack([xz, yz])
    cov = np.cov(z, rowvar=False, bias=False)
    cov = 0.5 * (cov + cov.T) + 1.0e-6 * np.eye(cov.shape[0])
    h = np.linalg.inv(cov)
    h = 0.5 * (h + h.T)
    target_idx = len(FEATURES)
    indices = [FEATURES.index(name) for name in observer] + [target_idx]
    c = np.zeros((len(indices), h.shape[0]))
    c[np.arange(len(indices)), indices] = 1.0
    phi = visible_precision(h, c, tolerances=Tolerances(atol=1e-9, rtol=1e-8))
    cond_var = float(1.0 / max(phi[-1, -1], 1e-12))
    r2 = float(max(0.0, min(0.999999, 1.0 - cond_var)))
    ceiling = h[np.ix_(indices, indices)]
    try:
        load = hidden_load(ceiling, phi, tolerances=Tolerances(atol=1e-8, rtol=1e-7))
        clock = float(load.clock)
        load_rank = float(load.rank)
        status = "ok"
    except NomogeoError as exc:
        clock = float("nan")
        load_rank = float("nan")
        status = f"load_unavailable:{type(exc).__name__}"
    return {
        "visible_r2": r2,
        "conditional_std": math.sqrt(max(cond_var, 0.0)),
        "hidden_clock": clock,
        "hidden_rank": load_rank,
        "geometry_status": status,
    }


def run_observer_gate(
    candles: list[Candle],
    train_window: int = 96,
    abstain_z: float = 0.15,
    source: str = "provided",
) -> dict[str, Any]:
    rows = build_rows(candles)
    if len(rows) <= train_window + 5:
        raise ValueError("not enough feature rows for requested train window")
    x_all = np.array([[row["features"][name] for name in FEATURES] for row in rows], dtype=float)
    y = np.array([row["target_return"] for row in rows], dtype=float)

    decisions: list[dict[str, Any]] = []
    observer_wins = {name: 0 for name in OBSERVERS}
    all_feature_mse: list[float] = []
    selected_mse: list[float] = []
    zero_mse: list[float] = []
    last_dir_hits: list[int] = []
    mom_dir_hits: list[int] = []
    selected_dir_hits: list[int] = []
    selected_spoke = 0

    for i in range(train_window, len(rows)):
        train_slice = slice(i - train_window, i)
        train_x_all = x_all[train_slice]
        train_y = y[train_slice]
        test_x_all = x_all[i]
        actual = float(y[i])
        train_std = float(np.std(train_y))
        if train_std <= 1e-12:
            train_std = 1.0

        observer_records = []
        for name, features in OBSERVERS.items():
            geom = _geometry(train_x_all, train_y, features)
            cols = [FEATURES.index(f) for f in features]
            train_x, test_x, _scale = _standardize(train_x_all[:, cols], test_x_all[cols])
            pred = _ridge_predict(train_x, train_y, test_x)
            observer_records.append(
                {
                    "name": name,
                    "features": features,
                    "prediction": pred,
                    "prediction_z": pred / train_std,
                    **geom,
                }
            )

        # Fixed, predeclared rule: choose the observer with the strongest visible
        # local prediction geometry. The forecast model itself is then tested
        # one step ahead, outside the train window.
        chosen = max(observer_records, key=lambda r: (float(r["visible_r2"]), -len(r["features"])))
        observer_wins[str(chosen["name"])] += 1
        selected_pred = float(chosen["prediction"])
        selected_mse.append((selected_pred - actual) ** 2)
        zero_mse.append(actual**2)

        all_train_x, all_test_x, _scale = _standardize(train_x_all, test_x_all)
        all_pred = _ridge_predict(all_train_x, train_y, all_test_x)
        all_feature_mse.append((all_pred - actual) ** 2)

        last_pred = float(test_x_all[FEATURES.index("ret_1")])
        mom_pred = float(test_x_all[FEATURES.index("mom_6")])
        if actual != 0.0:
            last_dir_hits.append(int(math.copysign(1.0, last_pred) == math.copysign(1.0, actual)))
            mom_dir_hits.append(int(math.copysign(1.0, mom_pred) == math.copysign(1.0, actual)))
            if abs(selected_pred / train_std) >= abstain_z:
                selected_spoke += 1
                selected_dir_hits.append(int(math.copysign(1.0, selected_pred) == math.copysign(1.0, actual)))

        decisions.append(
            {
                "time": rows[i]["time"],
                "time_iso": rows[i]["time_iso"],
                "close": rows[i]["close"],
                "actual_next_return": actual,
                "chosen_observer": chosen,
                "all_feature_prediction": all_pred,
                "baseline_last_return": last_pred,
                "baseline_momentum_6": mom_pred,
                "spoke": abs(selected_pred / train_std) >= abstain_z,
            }
        )

    latest = decisions[-1]
    selected_hit = statistics.fmean(selected_dir_hits) if selected_dir_hits else None
    last_hit = statistics.fmean(last_dir_hits) if last_dir_hits else None
    mom_hit = statistics.fmean(mom_dir_hits) if mom_dir_hits else None
    selected_rmse = math.sqrt(statistics.fmean(selected_mse))
    zero_rmse = math.sqrt(statistics.fmean(zero_mse))
    all_rmse = math.sqrt(statistics.fmean(all_feature_mse))
    best_direction_baseline = max(v for v in [last_hit, mom_hit, 0.5] if v is not None)
    gate_passes = (
        selected_hit is not None
        and selected_rmse < zero_rmse
        and selected_hit > best_direction_baseline + 0.01
        and selected_spoke >= 20
    )
    gate_reason = (
        "rolling gate cleared zero-return RMSE and directional baselines"
        if gate_passes
        else "no validated live edge: rolling selected view did not clear fixed baselines"
    )

    metrics = {
        "test_points": len(decisions),
        "selected_rmse": selected_rmse,
        "all_feature_rmse": all_rmse,
        "zero_return_rmse": zero_rmse,
        "selected_direction_accuracy_when_speaking": selected_hit,
        "selected_coverage": selected_spoke / max(1, len(decisions)),
        "last_return_direction_accuracy": last_hit,
        "momentum_6_direction_accuracy": mom_hit,
        "best_direction_baseline": best_direction_baseline,
        "gate_passes": gate_passes,
        "gate_reason": gate_reason,
        "observer_wins": observer_wins,
        "abstain_z": abstain_z,
        "train_window_hours": train_window,
    }
    return {
        "generated_at": utc_iso(int(time.time())),
        "source": source,
        "candle_count": len(candles),
        "feature_rows": len(rows),
        "first_candle": utc_iso(candles[0].time),
        "last_candle": utc_iso(candles[-1].time),
        "last_close": candles[-1].close,
        "features": FEATURES,
        "observers": {k: list(v) for k, v in OBSERVERS.items()},
        "metrics": metrics,
        "latest_decision": latest,
        "recent_decisions": decisions[-48:],
        "scope": (
            "Experimental paper-trading observation gate. It uses only hourly OHLCV-derived "
            "features and local Gaussian prediction geometry. It is not financial advice."
        ),
    }


def write_outputs(result: dict[str, Any], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    rows = result["recent_decisions"]
    lines = [
        "time_iso,close,actual_next_return,observer,prediction,prediction_z,visible_r2,hidden_clock,spoke"
    ]
    for row in rows:
        chosen = row["chosen_observer"]
        lines.append(
            "{},{:.8f},{:.10f},{},{:.10f},{:.6f},{:.6f},{},{}".format(
                row["time_iso"],
                row["close"],
                row["actual_next_return"],
                chosen["name"],
                chosen["prediction"],
                chosen["prediction_z"],
                chosen["visible_r2"],
                chosen["hidden_clock"],
                row["spoke"],
            )
        )
    (out_dir / "recent_decisions.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
