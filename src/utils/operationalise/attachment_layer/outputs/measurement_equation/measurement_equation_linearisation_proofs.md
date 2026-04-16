# Measurement Equation Linearisation Proofs v0

Status: generated proof packet for measurement-equation local linearisation.

These are not attachment cards. They test when a measurement equation can emit a local positive curvature object.

## beer_lambert_concentration_measurement_equation

- route: `measurement_equation_local_linearisation`
- status: `admissible_local_precision`
- equation: `c = A_abs / (epsilon path_length)`
- estimate: `4.8e-05`
- local standard uncertainty: `1.1931657238e-06`
- local precision: `702422577423`
- max sensitivity error: `1.20003000004e-09`
- boundary: This is a local concentration precision, not a storage Hessian; it must remain inside the calibrated optical range.

## cadmium_calibration_standard_measurement_equation

- route: `measurement_equation_local_linearisation`
- status: `admissible_local_precision`
- equation: `c_Cd = 1000 mass_cd_mg purity / volume_ml`
- estimate: `1002.69972`
- local standard uncertainty: `0.863702590151`
- local precision: `1.34051462347`
- max sensitivity error: `4.91323030793e-06`
- boundary: This is an uncertainty-propagation object for a calibration solution, not a physical hidden-load object.

## bromine_abundance_linearisation_negative_control

- route: `measurement_equation_local_linearisation`
- status: `first_order_refused`
- equation: `x160 = 2 x79 (1 - x79)`
- first-order variance: `0`
- second-order standard uncertainty: `7.07106781187e-05`
- allowed route: `higher_order_or_monte_carlo_only`
- lesson: A zero first derivative can make first-order uncertainty propagation report zero variance even when the nonlinear model has nonzero uncertainty.
