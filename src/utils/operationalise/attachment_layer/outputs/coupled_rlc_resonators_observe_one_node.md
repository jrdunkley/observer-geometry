# two coupled RLC resonators with one observed flux node and damping recorded outside storage

- card: `coupled_rlc_resonators_observe_one_node`
- attachment type: `exact_quadratic_with_declared_ceiling`
- derivation status: exact for the declared coupled LC storage model; exact damping response is not claimed
- equation hook: `U = flux_phi1^2/(2 L1) + flux_phi2^2/(2 L2) + (flux_phi1 - flux_phi2)^2/(2 Lc); R1,R2 are damping metadata outside U`
- coordinate note: Both active coordinates are flux coordinates in the same units; damping parameters are deliberately outside the Hessian and hidden-load call.
- allowed calls: visible_precision, hidden_load

## Scope

This card preserves the coupled-LC hidden-renormalisation story while making damping a recorded practical response feature rather than a static curvature term.

It does not describe dissipative linewidths, driven frequency response, nonlinear coupling, or measurement back-action.

## Sensor Surface

- ideal: flux_phi1
- practical: node voltage on resonator 1 integrated or phase-referenced to infer flux_phi1, linewidth or Q-factor as damping metadata

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
