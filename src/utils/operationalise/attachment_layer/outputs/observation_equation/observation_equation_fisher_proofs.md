# Observation Equation Fisher Proofs v0

Status: generated proof packet for observation-equation local Fisher precision.

These are not attachment cards. They test when declared observations plus a noise model expose local positive curvature.

## rc_decay_tau_observation_equation

- route: `observation_equation_fisher`
- status: `admissible_local_precision`
- equation: `V_i = V0 exp(-t_i / tau) + E_i, E_i ~ N(0, sigma^2)`
- object: `scalar Fisher precision`
- local precision: `2392.75348727`
- residual gate passed: `True`
- boundary: This compiles a sampled voltage trace into local information about tau; it is not capacitor storage or a hidden-load object.

## repeated_current_mean_observation_equation

- route: `observation_equation_fisher`
- status: `admissible_local_precision_with_estimated_sigma`
- equation: `I_j = I0 + E_j, E_j ~ N(0, sigma^2)`
- object: `scalar Fisher precision`
- local precision: `11148.2720178`
- residual gate passed: `True`
- boundary: This is a statistical observation model for a repeated measurement, not a physical storage law.

## load_cell_linear_calibration_observation_equation

- route: `observation_equation_fisher`
- status: `admissible_local_precision_matrix`
- equation: `D_i = alpha + beta F_i + E_i, E_i ~ N(0, sigma^2)`
- object: `2x2 Fisher precision matrix`
- SPD eigenvalues: `[4769.178718899588, 13760230.8212811]`
- residual gate passed: `True`
- boundary: This proves calibration-parameter precision for a declared linear response; using it as an inverse force-measurement function is a later analysis-function step.

## Negative Controls

### rc_decay_zero_sensitivity_grid_refusal

- status: `refused`
- reason: All samples are at t = 0, so dV/dtau is zero and no local precision for tau is exposed.
- refused route: `local_positive_curvature`

### load_cell_rank_deficient_design_refusal

- status: `refused`
- reason: All force standards have the same force value, so intercept and slope are not jointly identifiable.
- refused route: `local_positive_curvature_matrix`
