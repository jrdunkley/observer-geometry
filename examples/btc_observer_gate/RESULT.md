# Result

Live Coinbase run on 2026-04-11:

```text
mode: live_coinbase
candles: 300
window: 2026-03-30T06:00:00Z to 2026-04-11T17:00:00Z
last close: 73091.34
test points: 178
selected RMSE: 0.00451431
zero-return RMSE: 0.00419294
coverage: 53.4%
gate passes: false
```

Interpretation:

```text
No validated live edge. The selected observer view did not beat the fixed
zero-return and directional baselines on this rolling hourly window, so the UI
should refuse a directional BTC call.
```

That is the intended behaviour for a gate: it can say no.

Fixed-setting sweep on the same live candle pull:

```text
settings checked: 25
train windows: 48, 72, 96, 120, 144 hours
speak thresholds: 0.0, 0.15, 0.3, 0.5, 0.8
gate passes: 0
best selected/zero RMSE ratio: 1.0404 at 144 hours, threshold 0.0
best hit rate: 1.0 at 144 hours, threshold 0.8, but coverage was only 0.77%
```

Interpretation:

```text
No robust setting in this small fixed sweep cleared the baseline gate. The
apparently perfect hit-rate setting spoke on only one or about one case in the
test window, so it is rejected as a tiny-coverage artifact.
```
