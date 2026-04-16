# simple pendulum at the stable equilibrium under small-angle reduction

- card: `pendulum_small_angle`
- attachment type: `local_quadratic_approximation`
- derivation status: exact for the declared local quadratic approximation at the operating point
- equation hook: `U(theta) = m g ell (1 - cos theta) approx 1/2 m g ell theta^2 at theta0 = 0`
- coordinate note: The angular coordinate is treated as the declared local coordinate at theta0 = 0; radians are handled as dimensionless for the Hessian.
- allowed calls: visible_precision

## Scope

This card exposes the local stiffness of the pendulum at a declared stable operating point.

It does not describe global large-angle pendulum motion.

## Sensor Surface

- ideal: angle theta near zero
- practical: encoder, video tracking, or IMU angle estimate near the stable equilibrium

## Readout

```json
{
  "visible_precision": [
    [
      9.809999999999999
    ]
  ]
}
```

## Refusals

```json
[]
```
