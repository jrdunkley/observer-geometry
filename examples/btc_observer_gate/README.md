# BTC Observer Gate

Live paper-trading observation gate for hourly BTC-USD.

```bash
python -m examples.btc_observer_gate.run_main --live
python -m examples.btc_observer_gate.server
```

Then open:

```text
http://127.0.0.1:8765
```

## What It Does

The harness fetches recent public BTC-USD hourly candles from Coinbase Exchange,
builds fixed OHLCV-derived observer views, and runs a rolling one-hour-ahead
walk-forward check.

The UI reports:

- current last close
- selected observer view
- whether the gate is allowed to speak
- recent walk-forward coverage
- hit rate when speaking
- naive directional baselines
- zero-return and all-feature RMSE baselines

The gate refuses to present a directional call unless the selected view clears
fixed rolling validation checks against simple baselines.

## Data Source

Primary source:

```text
https://api.exchange.coinbase.com/products/BTC-USD/candles?granularity=3600
```

The Coinbase Exchange candle endpoint is public and unauthenticated. Coinbase
documents the candle schema, allowed granularities including `3600`, and a
single-request maximum of 300 candles:

```text
https://docs.cdp.coinbase.com/api-reference/exchange-api/rest-api/products/get-product-candles
```

Free public APIs can fail, rate-limit, or change behaviour. The harness includes
a synthetic offline mode for local validation:

```bash
python -m examples.btc_observer_gate.run_main
```

## Theory Surface

This is not a general trading model. It is an observer-geometry gate:

- estimate local Gaussian prediction geometry on a rolling train window
- compare fixed observer views of the hourly stream
- use `visible_precision` to compute the visible prediction geometry of each view
- use `hidden_load` to expose how much omitted feature structure still affects
  the chosen visible view
- test the selected view one hour ahead without lookahead

## Boundaries

- no financial advice
- no order execution
- no claim of durable market edge
- no hyperparameter search over the test window
- hourly OHLCV only
- local Gaussian/quadratic prediction geometry only
- Coinbase public API availability is not guaranteed

The useful failure mode is: "no validated live edge." That means the current
view does not justify a directional call under this harness.

See [RESULT.md](RESULT.md) for the latest checked live run.
