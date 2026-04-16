# Fifteen-Point Operationalisation Plan

Status: active plan derived from the five-step spine.

## Banked Five-Step Spine

1. Write the route admissibility contract.
2. Validate it against the current proof packets and cards.
3. Build one `calibration_analysis_function` proof, load-cell first or thermistor if richer inversion is needed.
4. Only then consider promoting a JSON schema for measurement attachments.
5. Defer thermodynamics until this measurement route is stable.

## Expanded Fifteen-Point Plan

1. Freeze the route vocabulary around actual successful routes, not speculative domains.
2. Separate storage curvature, information curvature, finite modal curvature, calibration analysis, and process gates.
3. Define the emitted object for each route before defining any card schema.
4. Define the allowed module surfaces for each route.
5. Define mandatory refusal conditions for each route.
6. Record the route contract in machine-readable JSON and human-readable Markdown.
7. Validate that current cards occupy only the storage routes they can justify.
8. Validate that measurement-equation proof packets occupy the local-linearisation route only.
9. Validate that observation-equation proof packets occupy the Fisher route only.
10. Validate that finite-mode proof packets require boundary, basis, truncation, and positive matrices.
11. Add a first calibration-analysis proof using a linear load-cell model and inverse analysis function.
12. Include calibration-analysis refusals for extrapolation and rank-deficient design.
13. Keep check-standard/process-control material as a gate, not a module object.
14. Keep thermodynamics and chemistry equilibrium deferred until ensemble/free-energy/fluctuation conventions are explicit.
15. Promote a measurement attachment schema only after the contract and route-specific proof packets remain boring under validation.

## Current Execution Target

This pass executes points 1-12 and records points 13-15 as binding constraints for the next pass.
