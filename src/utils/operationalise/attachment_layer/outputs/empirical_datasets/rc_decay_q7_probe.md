# Q7 RC Decay Probe v0

- dataset: `datasets/Q7_CapacitorVoltageDecay.m`
- route: `observation_equation_fisher`
- status: `admissible_empirical_probe`
- observations: 9

## Admission

The MATLAB problem data compile into a small RC exponential observation-equation probe after declaring V0 = 100 V, positive tau, and residual-estimated Gaussian voltage noise.

Refused routes:

- `static storage Hessian`
- `hidden_load`
- `capacitance inference without separately declared resistance`
- `instrument-noise claim without instrument metadata`

## Readout

- fixed V0: `100.0 V`
- tau_hat: `0.99711584 s`
- residual sigma: `1.1221254 V`
- Fisher precision for tau: `3931.5081`
- local variance for tau: `0.00025435532 s^2`
- max abs standardized residual: `1.278458`

## Log-Linear Comparison

- V0_hat: `95.806342 V`
- tau_hat: `1.0023593 s`
- decay rate: `-0.9976463 s^-1`

## Warnings

- nonmonotone point: 3.5 s, 2.0 V -> 4.0 s, 3.0 V

## Boundary

This is a small textbook regression dataset. It supports a local Fisher precision for tau only under a residual-estimated voltage-noise convention. It does not identify capacitance without a declared resistance, does not expose static storage Hessian geometry, and does not license hidden-load calls.
