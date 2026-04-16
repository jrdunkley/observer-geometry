# Calibration Analysis Function Proofs v0

Status: generated proof packet for calibration-analysis functions.

These are not attachment cards. They test when a calibrated sensor can expose an analysis function with local uncertainty.

## load_cell_linear_analysis_function

- route: `calibration_analysis_function`
- status: `admissible_analysis_function_record`
- calibration function: `D = alpha + beta F`
- analysis function: `F = (D - alpha) / beta`
- force estimate: `25.0006668724`
- local standard uncertainty: `0.432182739186`
- local precision: `5.35383724362`
- boundary: This proves inversion of a declared linear calibration inside its calibration range. It does not make the load cell a storage Hessian or hidden-load object.

## Negative Controls

### load_cell_analysis_extrapolation_refusal

- status: `refused`
- reason: A new indication outside the calibrated indication span cannot use the empirical analysis function without a separate extrapolation argument.
- refused route: `analysis_function_record`

### load_cell_analysis_rank_deficient_refusal

- status: `refused`
- reason: Repeated use of one force standard cannot identify both intercept and slope for the analysis function.
- refused route: `analysis_function_record`
