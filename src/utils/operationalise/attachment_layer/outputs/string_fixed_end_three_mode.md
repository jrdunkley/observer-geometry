# fixed-end string truncated to the first three sine modes

- card: `string_fixed_end_three_mode`
- attachment type: `exact_quadratic`
- derivation status: exact for the declared small-amplitude fixed-end string after basis choice and three-mode truncation
- equation hook: `U(q) = 1/2 sum_{n=1}^3 Tension (n pi / L_string)^2 q_n^2; M = mu I in the same orthonormal modal coordinates`
- coordinate note: This finite modal card separates stiffness Hessian K, mass metric M, and frequency readout. The runner reads the declared stiffness Hessian K; frequencies require the pair (K, M) and are not inferred from K alone.
- allowed calls: visible_precision

## Scope

This card converts a raw standing-wave formula into a finite positive modal stiffness object after boundary, basis, and truncation are declared.

It does not attach the infinite string, modes above n = 3, point-sensor inversion, damping, drive response, or hidden-load structure.

## Sensor Surface

- ideal: modal amplitudes q1, q2, q3
- practical: point displacement pickup, optical displacement tracking, or microphone-like channel after an additional modal observer map is declared

## Readout

```json
{
  "visible_precision": [
    [
      986.9604401089358,
      0.0,
      0.0
    ],
    [
      0.0,
      3947.8417604357433,
      0.0
    ],
    [
      0.0,
      0.0,
      8882.64396098042
    ]
  ]
}
```

## Refusals

```json
[]
```
