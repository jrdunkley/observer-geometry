# ideal gas equation of state without ensemble or fluctuation model

- card: `ideal_gas_law_boundary`
- attachment type: `not_yet_attachable`
- derivation status: not attachable from the equation of state alone
- equation hook: `p V = n R_gas T`
- coordinate note: The equation of state relates state variables but does not by itself declare a positive Hessian, covariance, Fisher object, or finite observer map.
- allowed calls: none

## Scope

This card records the boundary that an equation of state is not automatically a module object.

It cannot call theorem-local quadratic or hidden-load surfaces without a declared thermodynamic potential, ensemble, fluctuation model, or local measurement route.

## Sensor Surface

- ideal: pressure p, volume V, temperature T, amount n
- practical: pressure transducer, volume mark or piston displacement, thermometer, mole estimate from sample preparation

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
