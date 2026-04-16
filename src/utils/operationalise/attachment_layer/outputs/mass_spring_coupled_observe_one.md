# two coupled ideal spring coordinates with one observed coordinate

- card: `mass_spring_coupled_observe_one`
- attachment type: `exact_quadratic_with_declared_ceiling`
- derivation status: exact from the stated coupled ideal spring energy
- equation hook: `U = 1/2 k1 x1^2 + 1/2 k2 x2^2 + 1/2 kc (x1 - x2)^2`
- coordinate note: Both coordinates are displacements in the same units; entries of H can be compared as stiffnesses in this declared coordinate chart.
- allowed calls: visible_precision, hidden_load

## Scope

This card exposes how an unobserved coupled coordinate renormalises the visible stiffness.

It does not describe nonlinear springs, dynamic resonance, damping, or active feedback.

## Sensor Surface

- ideal: displacement x1
- practical: encoder or camera tracking x1 while x2 is uninstrumented

## Readout

```json
{
  "visible_precision": [
    [
      13.750000000000004
    ]
  ],
  "ceiling": [
    [
      15.0
    ]
  ],
  "ceiling_minus_visible_precision": [
    [
      1.2499999999999964
    ]
  ],
  "ceiling_dominance_min_eigenvalue": 1.2499999999999964,
  "ceiling_dominance_tolerance": 1.5e-09,
  "hidden_load": [
    [
      0.09090909090909083
    ]
  ],
  "hidden_rank": 1,
  "clock": 0.08701137698962969
}
```

## Refusals

```json
[]
```
