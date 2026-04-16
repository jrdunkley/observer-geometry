# two coupled LC resonator flux coordinates with one observed node

- card: `coupled_lc_resonators_observe_one_node`
- attachment type: `exact_quadratic_with_declared_ceiling`
- derivation status: exact after declaring the linear coupled-flux energy model
- equation hook: `U = flux_phi1^2/(2 L1) + flux_phi2^2/(2 L2) + (flux_phi1 - flux_phi2)^2/(2 Lc)`
- coordinate note: Both active coordinates are flux coordinates in the same units; this card does not mix charge and flux blocks.
- allowed calls: visible_precision, hidden_load

## Scope

This card exposes the same hidden-renormalised visible stiffness story as coupled springs, now in a circuit substrate.

It does not describe nonlinear coupling, dissipative linewidths, drive response, or measurement back-action.

## Sensor Surface

- ideal: flux_phi1
- practical: node voltage on resonator 1 integrated or phase-referenced to infer flux_phi1

## Readout

```json
{
  "visible_precision": [
    [
      14.000000000000007
    ]
  ],
  "ceiling": [
    [
      30.0
    ]
  ],
  "ceiling_minus_visible_precision": [
    [
      15.999999999999993
    ]
  ],
  "ceiling_dominance_min_eigenvalue": 15.999999999999993,
  "ceiling_dominance_tolerance": 3e-09,
  "hidden_load": [
    [
      1.1428571428571423
    ]
  ],
  "hidden_rank": 1,
  "clock": 0.7621400520468965
}
```

## Refusals

```json
[]
```
