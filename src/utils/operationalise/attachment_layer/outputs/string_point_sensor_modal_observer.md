# point displacement sensor over a three-mode fixed-end string truncation

- card: `string_point_sensor_modal_observer`
- attachment type: `exact_quadratic`
- derivation status: exact for the declared three-mode fixed-end truncation and point sensor position
- equation hook: `y = sum_{n=1}^3 sqrt(2 / L_string) sin(n pi x0 / L_string) q_n; U(q) = 1/2 q^T K q`
- coordinate note: This card is the first practical sensor-over-modal observer. The point sensor row is not full modal access; it induces a one-dimensional quotient stiffness for the measured displacement channel. Frequencies still require the separate pair (K, M), and hidden load still requires a declared ceiling/reference.
- allowed calls: visible_precision

## Scope

This card shows that a practical point sensor becomes a declared observer map over a finite modal object and exposes a one-dimensional quotient stiffness.

It does not infer full modal amplitudes, attach modes above n = 3, model sensor noise, attach damping or drive response, or license hidden-load structure.

## Sensor Surface

- ideal: point displacement y at x0 = L_string / 4 as a linear functional of q1, q2, q3
- practical: single optical displacement pickup, contact displacement probe, or microphone-like point channel calibrated at x0 = L_string / 4

## Readout

```json
{
  "visible_precision": [
    [
      612.5961352400288
    ]
  ]
}
```

## Refusals

```json
[]
```
