# full simple pendulum potential over finite angle range

- card: `pendulum_full`
- attachment type: `not_yet_attachable`
- derivation status: not attachable as a single exact quadratic object
- equation hook: `U(theta) = m g ell (1 - cos theta)`
- coordinate note: No local coordinate chart or quadratic operating point is declared for module execution.
- allowed calls: none

## Scope

This card records the boundary: a famous physical formula is not automatically an observer-geometry object.

It cannot call theorem-local quadratic kernel surfaces until a local operating point, sampled family, or other declared object is supplied.

## Sensor Surface

- ideal: angle theta over a finite range
- practical: encoder, video tracking, or IMU angle estimate over finite amplitude

## Readout

```json
{}
```

## Refusals

```json
[
  {
    "call": "*",
    "reason": "card is not_yet_attachable; no theorem-local kernel calls are licensed"
  }
]
```
