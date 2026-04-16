# single ideal capacitor in declared charge coordinate

- card: `capacitor_charge_coordinate`
- attachment type: `exact_quadratic`
- derivation status: exact after declaring charge as the capacitor energy coordinate
- equation hook: `U(q) = q^2/(2 Ccap)`
- coordinate note: This card uses charge q as the coordinate. If voltage V is used as the coordinate, U(V) = 1/2 Ccap V^2 gives a different Hessian with different units; the two coordinates must not be compared as the same stiffness table.
- allowed calls: visible_precision

## Scope

This card exposes capacitor storage as a scalar quadratic module object in a declared charge coordinate.

It does not describe voltage-coordinate geometry, leakage, nonlinear dielectrics, fringing fields, or dynamic RC/RLC response.

## Sensor Surface

- ideal: charge q
- practical: capacitor voltage V = q / Ccap measured with a high-impedance voltmeter

## Readout

```json
{
  "visible_precision": [
    [
      999999999.9999999
    ]
  ]
}
```

## Refusals

```json
[]
```
