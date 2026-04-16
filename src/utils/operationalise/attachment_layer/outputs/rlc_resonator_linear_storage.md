# single series RLC resonator with declared LC storage coordinates

- card: `rlc_resonator_linear_storage`
- attachment type: `exact_quadratic`
- derivation status: exact for the declared LC storage energy; damping is deliberately excluded from H
- equation hook: `U(q, flux_phi) = q^2/(2 Ccap) + flux_phi^2/(2 L); damping metadata R is outside U`
- coordinate note: Charge and flux have different units; resistance is a dissipative response parameter and must not be compared with Hessian entries.
- allowed calls: visible_precision

## Scope

This card exposes the storage part of a damped resonator as the same exact quadratic LC object while keeping damping out of the Hessian.

It does not describe dissipative linewidth, driven response, transient fitting, or energy loss as a static precision object.

## Sensor Surface

- ideal: charge q and flux_phi as storage coordinates
- practical: capacitor voltage V = q / Ccap, inductor current I = flux_phi / L, ring-down linewidth or Q-factor as damping metadata

## Readout

```json
{
  "visible_precision": [
    [
      10000.0,
      0.0
    ],
    [
      0.0,
      9.999999999999998
    ]
  ]
}
```

## Refusals

```json
[]
```
