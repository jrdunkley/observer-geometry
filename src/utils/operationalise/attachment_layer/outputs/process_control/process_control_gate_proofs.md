# Process-Control Gate Proofs v0

Status: generated proof packet for measurement-process gates.

These are not attachment cards. They test whether a measurement route may proceed at all. A passing gate emits a validity decision only; it does not license theorem-local kernel calls.

## check_standard_process_gate

- route: `process_control_gate`
- status: `gate_passed_no_module_call`
- emitted object: measurement-process validity decision only
- observed mean: `10.001125`
- max absolute z: `1.1`
- action limit z: `3`
- theorem-local calls licensed: `false`
- boundary: This gate says the measurement process is not obviously out of control under the declared check-standard rule. It does not create a module Hessian, Fisher object, or hidden-load object.

## Negative Controls

### check_standard_drift_refusal

- status: `refused`
- reason: At least one check-standard observation exceeds the declared three-sigma action limit.
- refused route: `gate_or_refusal_only`
- theorem-local calls licensed: `false`

### check_standard_missing_baseline_refusal

- status: `refused`
- reason: A process-control gate cannot be evaluated without a stable check standard, historical baseline, and variability declaration.
- refused route: `gate_or_refusal_only`
- theorem-local calls licensed: `false`
