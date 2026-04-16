# translational and rotational kinetic energy in declared velocity coordinates

- card: `kinetic_and_rotational_energy_metrics`
- attachment type: `exact_quadratic`
- derivation status: exact from the stated nonrelativistic kinetic energy after the velocity coordinates and rotation axis are declared
- equation hook: `K(v, omega) = 1/2 m v^2 + 1/2 I_axis omega^2`
- coordinate note: This is kinetic metric geometry in velocity coordinates. The Hessian entries are mass and moment of inertia, not restoring stiffnesses; they must not be compared entry-by-entry with spring or LC storage Hessians.
- allowed calls: visible_precision

## Scope

This card admits kinetic energy as an exact finite positive quadratic object while keeping metric geometry distinct from stiffness and frequency geometry.

It does not describe potentials, oscillation frequencies, damping, relativistic kinetic energy, gyroscopic coupling, or hidden-load structure.

## Sensor Surface

- ideal: translational velocity v, angular velocity omega
- practical: motion tracking or encoder-derived velocity, tachometer, gyroscope, or rotary encoder angular velocity

## Readout

```json
{
  "visible_precision": [
    [
      2.0000000000000004,
      0.0
    ],
    [
      0.0,
      0.08
    ]
  ]
}
```

## Refusals

```json
[]
```
