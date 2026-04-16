# Non-Integration Diagnostics v0

Status: every non-carded hook in `candidate_backlog.json` has been routed by why it is not a v0 static attachment card.

The key result is that most failures are not theory failures. They are missing attachment declarations or misrouting relative to the first runner. The v0 runner is intentionally narrow: `visible_precision` and `hidden_load`. The broader module stack has other possible homes: evidence assembly, quotient descent, local quadratic ensembles, weighted-family frontiers, support-stratum transport, affine-hidden sectors, covariance/Fisher perturbations, and residual-margin diagnostics.

## Route Summary

| Route Class | Count | Meaning |
| --- | ---: | --- |
| `static_quadratic_next` | 0 | Ready or near-ready once coordinate/units details are declared. |
| `local_quadratic` | 3 | Needs operating point, sign/stability convention, or local expansion. |
| `dynamic_evidence` | 6 | Time-series, likelihood, Fisher, regression, or fit object rather than static Hessian. |
| `sensor_or_observer_map` | 6 | Measurement or observer map before module object. |
| `finite_weighted_family` | 1 | Needs finite basis/truncation/family declaration. |
| `thermodynamic_fluctuation` | 4 | Needs ensemble, potential, coordinates, and fluctuation/local Hessian convention. |
| `support_or_event` | 1 | Event/conservation route rather than standing quadratic object. |
| `exact_special_sector` | 2 | Requires special law sector; not generic v0. |
| `parameter_input` | 1 | Useful input, not standalone object. |
| `true_gap_or_defer` | 0 | No item is currently assigned to this bucket as its primary route, though relativity is a real SPD-boundary case. |

## Immediate Consequences

The apparent weak integrations divide cleanly:

- RC decay, chemistry kinetics, nuclear decay, blackbody fitting, and Hubble law are probably evidence/Fisher objects.
- Absorbance, optics, op-amps, and photon relations are sensor or observer-map objects.
- Acoustic hooks now have a narrow positive finite-modal pressure-observer card, but raw impedance, reflection, attenuation, and ultrasound imaging surfaces remain refused without additional declarations. Fluids still need finite control-volume declarations. Standing string waves have moved into the carded finite-modal set; the point-sensor modal observer and a reference-ceiling hidden-load variant are also carded. Non-reference hidden-load variants still need a declared ceiling/reference.
- Heat capacity, equilibrium constants, electrochemistry, and thermodynamic potentials need a declared thermodynamic/fluctuation schema. Ideal gas is now carded as an explicit refusal of the raw equation-of-state surface.
- Relativity should remain deferred unless a positive exact special sector is isolated.
- Quantum oscillator is promising but should wait until the classical oscillator/circuit vocabulary is stable.

## Strongest Near-Term Non-Static Route

The best next diagnostic target is probably `rc_charge_decay` or `absorbance_beer_lambert`, not because they are deeper than thermodynamics, but because they force the smallest useful evidence/sensor schema.

`rc_charge_decay` tests whether a time series can become a declared Fisher object.

`absorbance_beer_lambert` tests whether a sensor law can become a measurement Fisher object.

Electrochemistry remains the strongest chemistry bridge, but it has more hidden declarations: reaction coordinate, standard state, temperature, electron count, reaction quotient, and voltage noise.

## Current Boundary

Do not call these objects impossible. Also do not card them as static quadratic objects. The right next move is to define one minimal evidence/sensor attachment schema, then try one of the small hooks against it.

Note: `linear_elastic_bar`, `capacitor_charge_coordinate`, `rlc_resonator_linear_storage`, and `coupled_rlc_resonators_observe_one_node` have now moved out of this non-integration set into the carded storage-curvature set. `kinetic_and_rotational_energy_metrics` has moved into the carded metric-curvature set. `string_fixed_end_three_mode`, `string_point_sensor_modal_observer`, `string_modal_hidden_load_with_declared_ceiling`, and `acoustic_tube_three_mode_pressure_observer` have moved into the carded finite-modal set. `ideal_gas_law_boundary` has moved into the carded refusal set.

## Failure Probe Suite

The generated probe suite lives in:

- `outputs/failure_probes/failure_probe_suite.json`
- `outputs/failure_probes/failure_probe_suite.md`

Representative lessons:

| Probe | Lesson |
| --- | --- |
| `rc_charge_decay` | With sample times and Gaussian voltage noise, the decay law yields a scalar Fisher precision for `tau`; without those declarations, it is not a static card. |
| `absorbance_beer_lambert` | With absorbance noise and concentration parameter, the sensor law yields Fisher precision; without noise, it remains a measurement map. |
| `nuclear_decay_poisson` | With counting windows and Poisson likelihood, decay yields Fisher precision for `lambda`; without likelihood, it is just a law surface. |
| `standing_wave_string_modes` | Now carded as `string_fixed_end_three_mode`, `string_point_sensor_modal_observer`, and `string_modal_hidden_load_with_declared_ceiling`; non-reference hidden-load versions still require a declared ceiling/reference. |
| `gravity_inverse_square_local` | Raw radial gravitational curvature can have the wrong sign; effective-potential/stability declarations are decisive. |
| `ideal_gas_law_boundary` | A Helmholtz isothermal local Hessian can be produced only after ensemble and potential declaration; the gas law alone remains a refusal. |
| `electrochemistry_nernst_voltage` | Voltage is an excellent sensor surface, but it needs reaction coordinate, thermodynamic convention, and noise model. |
| `relativity_lorentz_metric` | The Lorentz metric is indefinite, so it remains a real boundary for current SPD calls. |
