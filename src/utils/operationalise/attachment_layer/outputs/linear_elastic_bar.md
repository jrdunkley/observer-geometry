# small-strain linear elastic bar in axial extension coordinate

- card: `linear_elastic_bar`
- attachment type: `exact_quadratic`
- derivation status: exact from the stated small-strain uniform linear elastic energy
- equation hook: `U(delta_L) = 1/2 (E_young A_cross / L0) delta_L^2`
- coordinate note: The coordinate is total axial extension, not strain. If strain epsilon = delta_L / L0 is used as the coordinate, the Hessian changes to E_young A_cross L0.
- allowed calls: visible_precision

## Scope

This card exposes material stiffness as the same scalar quadratic contact type as an ideal spring after geometry is declared.

It does not describe plastic deformation, nonlinear elasticity, bending modes, shear, damping, or fracture.

## Sensor Surface

- ideal: axial extension delta_L
- practical: strain gauge, displacement probe, or load-frame extension measurement

## Readout

```json
{
  "visible_precision": [
    [
      199999.99999999994
    ]
  ]
}
```

## Refusals

```json
[]
```
