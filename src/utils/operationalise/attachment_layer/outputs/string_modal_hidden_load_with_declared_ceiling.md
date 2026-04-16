# point-displacement string modal observer with a declared reference-shape ceiling

- card: `string_modal_hidden_load_with_declared_ceiling`
- attachment type: `exact_quadratic_with_declared_ceiling`
- derivation status: exact for the declared three-mode fixed-end truncation, point sensor position, and reference-shape right inverse
- equation hook: `y = C q; r = C^T / (C C^T); T = r^T K r; Phi = (C K^{-1} C^T)^{-1}`
- coordinate note: This card is not a visible-block ceiling in the original modal coordinates. The ceiling is a declared reference-shape stiffness for the unit point displacement channel; hidden-load readout is therefore a reference-relaxation statement, not an intrinsic property of C alone.
- allowed calls: visible_precision, hidden_load

## Scope

This card shows that a point-sensor modal observer can support a hidden-load readout only after a reference shape and ceiling are declared.

It does not claim that the point sensor alone defines a canonical ceiling, infer full modal amplitudes, attach modes above n = 3, model noise, or attach damping and drive response.

## Sensor Surface

- ideal: point displacement y at x0 = L_string / 4 as a linear functional of q1, q2, q3
- practical: calibrated point displacement pickup plus a declared reference actuation or clamping protocol for the modal shape r

## Readout

```json
{
  "visible_precision": [
    [
      612.5961352400288
    ]
  ],
  "ceiling": [
    [
      1110.330495122552
    ]
  ],
  "ceiling_minus_visible_precision": [
    [
      497.7343598825232
    ]
  ],
  "ceiling_dominance_min_eigenvalue": 497.7343598825232,
  "ceiling_dominance_tolerance": 1.1103304951225521e-07,
  "hidden_load": [
    [
      0.8124999999999993
    ]
  ],
  "hidden_rank": 1,
  "clock": 0.5947071077466924
}
```

## Refusals

```json
[]
```
