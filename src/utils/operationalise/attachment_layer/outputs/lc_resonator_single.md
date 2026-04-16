# single ideal LC resonator in declared energy coordinates

- card: `lc_resonator_single`
- attachment type: `exact_quadratic`
- derivation status: exact after declaring charge and flux as the energy coordinates
- equation hook: `U(q, flux_phi) = q^2/(2 Ccap) + flux_phi^2/(2 L)`
- coordinate note: Charge and flux have different units; the Hessian is meaningful in the declared energy coordinates, not as a homogeneous stiffness table.
- allowed calls: visible_precision

## Scope

This card shows a circuit resonator as the same exact quadratic contact type as a mechanical oscillator.

It does not describe resistance, nonlinear components, drive, radiation, or measurement loading.

## Sensor Surface

- ideal: charge q and flux_phi as energy coordinates
- practical: capacitor voltage V = q / Ccap, inductor current I = flux_phi / L

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
