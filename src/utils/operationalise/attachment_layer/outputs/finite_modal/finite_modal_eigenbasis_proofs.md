# Finite Modal Eigenbasis Proofs v0

Status: generated proof packet for finite modal attachment.

These are not attachment cards. They test when a wave surface becomes a finite positive matrix object.

## string_fixed_end_three_mode

- route: `finite_modal_eigenbasis`
- status: `admissible_finite_modal_family`
- stiffness eigenvalues: `[986.9604401089358, 3947.8417604357433, 8882.643960980422]`
- frequencies Hz: `[50.0, 100.0, 150.0]`
- boundary: This proof admits only the declared three-mode truncation. A physical sensor still needs an observer map before visible precision or hidden-load claims.

## string_point_sensor_two_mode_projection

- route: `finite_modal_eigenbasis`
- status: `observer_map_recorded_not_hidden_load`
- observer map: `[1.0000000000000002, 1.4142135623730951]`
- lesson: A point sensor exposes a projection of modal coordinates. It is an observer map, but it is not enough for hidden-load without a declared ceiling/reference.

## acoustic_tube_three_mode_pressure_observer

- route: `finite_modal_eigenbasis`
- status: `admissible_finite_modal_pressure_observer`
- stiffness eigenvalues: `[13933.789058205142, 55735.15623282057, 125404.10152384629]`
- frequencies Hz: `[171.49999999999997, 342.99999999999994, 514.5]`
- boundary: This proof admits only the declared three-mode closed-tube pressure observer. It does not admit impedance, reflection, attenuation, damping, or full field reconstruction.

## Negative Controls

### raw_standing_wave_no_truncation_refusal

- status: `refused`
- reason: A standing-wave formula without boundary conditions, basis, and finite truncation is not a finite matrix object.
- refused route: `finite_modal_positive_curvature`

### string_zero_tension_stiffness_refusal

- status: `refused`
- reason: Zero tension gives zero stiffness eigenvalues, so the declared stiffness Hessian is not positive.
- refused route: `finite_modal_positive_curvature`
