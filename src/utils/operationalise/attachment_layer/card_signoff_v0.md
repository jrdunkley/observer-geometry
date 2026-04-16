# Card Sign-Off v0

Status: rigorous v0 sign-off for the current attachment card set.

This does not claim empirical validation of a physical apparatus. It signs off a card as a good attachment object: the conversion from formula-sheet language into declared finite observer-geometry input is explicit, runnable, validated, and honest about scope.

## Good Means

A v0 card is `good` only in this restricted sense:

- it declares physical coordinates, units, sensor surface, controls, regime, and failure mode;
- its equation hook actually yields the declared module object;
- its derivation status does not inflate the exactness claim;
- it licenses only the module calls it can support;
- the runner enforces those calls;
- hidden-load cards declare a ceiling and validate `T >= Phi` before calling `hidden_load`;
- generated output exposes what the module sees and what it does not see;
- at least one practical next measurement/control action is stated;
- refusal cards refuse theorem-local calls rather than pretending a surface formula is enough;
- validation passes with a negative/refusal check.

## Current Sign-Off Table

| Card | Sign-Off | Reason | Remaining Non-Blocking Work |
| --- | --- | --- | --- |
| `mass_spring_single` | `good_v0_baseline` | Exact quadratic stiffness card with declared displacement, stiffness units, sensor surface, controls, and licensed `visible_precision`. | Later add an experimental period/stiffness check if moving toward apparatus. |
| `lc_resonator_single` | `good_v0_baseline` | Exact quadratic storage card with declared charge/flux coordinates and explicit mixed-coordinate warning. | Later add a sibling capacitor charge/voltage coordinate-choice card. |
| `linear_elastic_bar` | `good_v0_storage` | Exact small-strain material stiffness card with declared extension coordinate, geometry-normalized stiffness, units, sensor surface, and licensed `visible_precision`. | Later add strain-coordinate sibling only if the coordinate transformation lesson is needed explicitly. |
| `capacitor_charge_coordinate` | `good_v0_storage` | Exact capacitor storage card in declared charge coordinate with explicit warning that voltage-coordinate Hessian has different units and meaning. | Later add voltage-coordinate sibling only if the sensor-surface ladder requires it. |
| `rlc_resonator_linear_storage` | `good_v0_storage` | Exact LC storage Hessian with resistance declared as damping metadata outside `H`, so the card broadens the circuit analogue without pretending damping is static curvature. | Later add driven-response or linewidth evidence card only after a dynamic convention is declared. |
| `kinetic_and_rotational_energy_metrics` | `good_v0_metric` | Exact kinetic metric card with declared velocity and angular-velocity coordinates, mass/inertia units, and an explicit boundary against treating the metric as stiffness or frequency geometry. | Use this as the vocabulary gate before promoting finite modal cards. |
| `string_fixed_end_three_mode` | `good_v0_finite_modal` | Exact three-mode fixed-end string stiffness card with boundary, sine basis, truncation, mass metric metadata, and explicit separation of `K`, `M`, and frequency readout. | Point-sensor observer and reference-ceiling hidden-load variants now exist; later add only non-reference variants with a declared ceiling/reference. |
| `string_point_sensor_modal_observer` | `good_v0_modal_observer` | Exact one-row point-displacement observer over the three-mode string card, exposing a quotient stiffness while refusing full modal access and hidden-load structure. | Later add noise/Fisher variants only with a declared sensor model. |
| `string_modal_hidden_load_with_declared_ceiling` | `good_v0_modal_reference` | Exact point-sensor modal hidden-load card with a declared reference-shape ceiling, making the stiffness gap a reference-relaxation statement rather than an intrinsic property of the sensor alone. | Later add a sweep over sensor position or reference-shape protocols; do not treat this as a canonical visible-block ceiling. |
| `acoustic_tube_three_mode_pressure_observer` | `good_v0_modal_observer` | Exact small-signal closed-tube acoustic modal card with rigid boundary, finite displacement basis, mass metric metadata, and a pressure microphone observer row. | Later add damping, drive response, impedance, or reflection only after a dynamic/sensor evidence route is declared. |
| `mass_spring_coupled_observe_one` | `good_v0_bridge` | Exact coupled quadratic with declared one-coordinate observer, visible-block ceiling, `T - Phi` readout, hidden load, clock, and coupling sweep. | Later add apparatus-level protocol for clamped versus relaxed visible stiffness. |
| `coupled_lc_resonators_observe_one_node` | `good_v0_bridge` | Exact coupled LC analogue of the spring bridge with declared flux coordinates, ceiling, `T - Phi` readout, hidden load, clock, and coupling sweep. | Later add RLC damping wrapper without putting damping into static `H`. |
| `coupled_rlc_resonators_observe_one_node` | `good_v0_bridge` | Exact coupled LC storage bridge with declared damping metadata outside `H`; the licensed hidden-load call uses only the storage Hessian and declared visible-block ceiling. | Later add frequency-response fitting only after the evidence layer is explicit. |
| `pendulum_small_angle` | `good_v0_regime` | Local quadratic approximation at declared operating point with tightened derivation status and no global overclaim. | Later add amplitude-bound convention if experimental data are attached. |
| `pendulum_full` | `good_v0_refusal` | Refuses theorem-local calls for the full nonlinear pendulum, preserving the rule that a formula is not yet an observer object. | None for v0; keep refusal strict. |
| `ideal_gas_law_boundary` | `good_v0_refusal` | Refuses theorem-local calls for `pV = nRT` by itself and records the missing ensemble, potential, fluctuation model, or measurement route. | Later revisit only through a thermodynamic or evidence attachment grammar. |

## Coupled Analogue Sweep

The generated sweep lives in:

- `outputs/sweeps/coupled_hidden_renormalisation_sweep.json`
- `outputs/sweeps/coupled_hidden_renormalisation_sweep.csv`
- `outputs/sweeps/coupled_hidden_renormalisation_sweep.md`

The sweep varies coupling in the coupled spring and coupled LC cards while preserving the one-visible, one-hidden structure.

Checked monotone fields:

- `ceiling`
- `ceiling_minus_visible_precision`
- `hidden_load`
- `clock`

Observed endpoints:

| Family | Coupling Range | `T - Phi` Range | Hidden Load Range | Clock Range |
| --- | ---: | ---: | ---: | ---: |
| coupled spring | `0.5 -> 20` | `0.0161290322581 -> 11.4285714286` | `0.00153846153846 -> 0.615384615385` | `0.00153727931889 -> 0.479573080262` |
| coupled LC | `1 -> 40` | `0.166666666667 -> 35.5555555556` | `0.0153846153846 -> 2.46153846154` | `0.0152674721308 -> 1.24171313231` |

Interpretation: both substrates realise the same qualitative observer-geometry story. A hidden passive coordinate relaxes the visible response below the visible-block ceiling, and the hidden-load signature strengthens with coupling.

This is the first `good_v0_bridge`: cross-substrate quadratic storage under partial observation.

## Validation Commands

```powershell
python validate_cards.py --nomogeo-src C:\observer_geometry_workspace_v0.3.2\observer_geometry\src
python run_card.py --all --markdown --nomogeo-src C:\observer_geometry_workspace_v0.3.2\observer_geometry\src
python analogue_sweep.py --nomogeo-src C:\observer_geometry_workspace_v0.3.2\observer_geometry\src
```

## Boundary

This sign-off does not yet cover:

- live sensor feeds;
- damping or driven response;
- damping is recorded on the RLC storage cards, but dynamic response itself is not yet attached;
- nonlinear large-amplitude regimes;
- empirical calibration;
- dynamic/evidence/Fisher cards;
- thermodynamic or chemistry fluctuation models;
- global observer discovery.

Those belong to the next diagnostic track, not to the current v0 good-card sign-off.
