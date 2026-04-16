# Next Schema Decisions v0

Status: planning artifact derived from the failure diagnostics, failure probe suite, and local-positive-curvature admissibility memo.

The failure probes show that the non-integrating cases are useful. Director correction: the schema is not the idea. The idea is declared local positive curvature admissibility. Schema work should follow only after the admissibility idea is tested.

## Pre-Decision: Local Positive Curvature Admissibility

Priority: immediate.

Why: the first successful layer is stored quadratic curvature. The failures suggest a broader class: declared local positive curvature, where the curvature may come from energy, likelihood, free energy, finite modes, or special sectors.

Already tested in:

- `rc_decay_local_curvature`
- `beer_lambert_local_curvature`

These prove-of-admissibility examples refuse the bare formula, declare the local coordinate/noise/sensor route, and produce a positive scalar Fisher object.

## Decision 1: Minimal Evidence/Fisher Attachment

Priority: high, but only after the admissibility proofs remain stable.

Why: `rc_charge_decay`, `absorbance_beer_lambert`, and `nuclear_decay_poisson` all become finite scalar Fisher objects after a small number of explicit declarations.

Candidate schema fields:

- `evidence_type`: `gaussian_time_series`, `gaussian_sensor_law`, or `poisson_counting_process`
- `estimated_parameters`
- `sensor_surface`
- `sampling_surface`
- `noise_model`
- `formula_hook`
- `fisher_object`
- `derivation_status`
- `allowed_module_route`: initially `evidence_only`, not `hidden_load`
- `refusal_without_noise_model`

First trial candidate after admissibility: `rc_charge_decay_fisher_tau`.

Reason: It is physically simple, uses circuit variables already in the bundle, and forces the dynamic/static separation cleanly.

Second trial candidate after admissibility: `absorbance_beer_lambert_fisher_concentration`.

Reason: It is the smallest chemistry sensor hook and avoids thermodynamic ambiguity.

## Decision 2: Finite Mode Family Attachment

Priority: high, after one evidence/Fisher trial.

Why: `standing_wave_string_modes` shows that a wave hook fails mostly because it is not finite until boundary, basis, and truncation are declared.

Candidate schema fields:

- `boundary_conditions`
- `basis`
- `truncation`
- `mass_metric`
- `stiffness_family`
- `sensor_surface`
- `control_variables`
- `allowed_nomogeo_calls`: likely no hidden-load unless an observer/ceiling is declared

First trial candidate: `string_fixed_end_three_mode`, now promoted as the first finite-modal card.

Reason: It is the cleanest route from AP wave formulas into finite quadratic geometry. The point-sensor modal observer and reference-ceiling hidden-load variant are now also promoted. The first acoustic analogue is now `acoustic_tube_three_mode_pressure_observer`, which admits only a closed-tube finite modal pressure observer while continuing to refuse raw impedance, attenuation, reflection, and ultrasound imaging surfaces.

## Decision 3: Thermodynamic/Fluctuation Attachment

Priority: medium; do not rush.

Why: `ideal_gas_law_boundary` and `electrochemistry_nernst_voltage` both show that thermodynamics can produce local Hessian or Fisher objects only after ensemble, potential, coordinate, and noise declarations.

Candidate schema fields:

- `ensemble`
- `thermodynamic_potential`
- `state_coordinate`
- `operating_point`
- `convexity_or_sign_convention`
- `fluctuation_or_noise_model`
- `sensor_surface`
- `refusal_without_ensemble`

First trial candidate: `ideal_gas_isothermal_helmholtz_boundary`.

Reason: It should remain a refusal-plus-conditional local Hessian lesson, not a broad thermodynamics claim.

Best chemistry candidate after that: `nernst_voltage_logQ_fisher`.

Reason: Voltage gives a practical sensor surface, but only after reaction-coordinate declarations.

## Decision 4: Local Field Patch Attachment

Priority: medium-low.

Why: `gravity_inverse_square_local` shows that raw field formulas can have the wrong local sign, while an effective-potential declaration can produce a positive local Hessian.

Candidate schema fields:

- `raw_potential`
- `operating_point`
- `local_coordinate`
- `stability_convention`
- `effective_potential_terms`
- `local_hessian`
- `refusal_without_stability`

First trial candidate: not yet. Capacitor and elastic cards are cleaner.

## Decision 5: Hard SPD Boundary

Priority: defer.

Why: `relativity_lorentz_metric` is not a missing-unit problem. The metric is indefinite. Current static `visible_precision`/`hidden_load` calls are built for SPD objects.

Rule:

Do not route Lorentzian metric material into current SPD calls unless a positive subproblem or exact special sector is explicitly declared.

## Recommended Next Work Order

1. Preserve v0 bridge stability:
   - keep `mass_spring_coupled_observe_one` and `coupled_lc_resonators_observe_one_node` as signed-off bridge cards.

2. Preserve the newly promoted direct static cards:
   - `linear_elastic_bar`
   - `capacitor_charge_coordinate`

3. Promote local positive curvature admissibility into one minimal evidence/Fisher prototype:
   - `rc_charge_decay_fisher_tau`

4. Add one minimal sensor-law Fisher prototype:
   - `absorbance_beer_lambert_fisher_concentration`

5. Extend the finite mode family:
   - `acoustic_modal_observer_with_declared_boundary`

6. Only then approach thermodynamics:
   - `ideal_gas_isothermal_helmholtz_boundary`
   - later `nernst_voltage_logQ_fisher`

This order preserves the v0 success while letting the failures teach the next attachment vocabulary.
