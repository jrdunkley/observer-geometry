# single ideal mass-spring oscillator near equilibrium

- card: `mass_spring_single`
- attachment type: `exact_quadratic`
- derivation status: exact from the stated ideal spring equation
- equation hook: `U(x) = 1/2 k x^2`
- coordinate note: Single homogeneous displacement coordinate; the Hessian entry is directly interpretable as stiffness.
- allowed calls: visible_precision

## Scope

This card exposes physical stiffness as a direct one-dimensional quadratic module object.

It does not describe nonlinear springs, damping, driving, or measurement noise.

## Sensor Surface

- ideal: displacement x
- practical: linear displacement sensor, camera tracking, or encoder measuring x

## Readout

```json
{
  "visible_precision": [
    [
      11.999999999999995
    ]
  ]
}
```

## Refusals

```json
[]
```
