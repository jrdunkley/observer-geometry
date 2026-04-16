# closed-closed acoustic tube truncated to three displacement modes with a pressure microphone observer

- card: `acoustic_tube_three_mode_pressure_observer`
- attachment type: `exact_quadratic`
- derivation status: exact for the declared small-signal one-dimensional closed-tube model after basis choice and three-mode truncation
- equation hook: `U(xi) = 1/2 rho_air c_sound^2 A_tube integral_0^L (d xi / dx)^2 dx; p(x0) = -rho_air c_sound^2 sum_n xi_n d phi_n/dx at x0`
- coordinate note: This card uses longitudinal displacement modal amplitudes as the state coordinates and a pressure microphone as the observer row. The pressure quotient precision is meaningful only for this declared pressure channel; it is not acoustic impedance, intensity, reflection coefficient, or full field reconstruction.
- allowed calls: visible_precision

## Scope

This card shows that an acoustic resonator can expose a finite modal positive object through a realistic pressure sensor after boundary, basis, truncation, and sensor convention are declared.

It does not attach raw sound intensity, ultrasound reflection, acoustic impedance, damping, drive response, sensor noise, or full pressure-field reconstruction.

## Sensor Surface

- ideal: pressure p at x0 = L_tube / 3 as a linear functional of xi1, xi2, xi3
- practical: small microphone calibrated as a point pressure sensor at x0 = L_tube / 3 within the low-amplitude acoustic regime

## Readout

```json
{
  "visible_precision": [
    [
      2.3610721534205812e-08
    ]
  ]
}
```

## Refusals

```json
[]
```
