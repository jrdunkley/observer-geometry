# Stack Validation Note

This note records the v0.2 validation and soak pass across the current stack:

- `nomogeo`
- `nomodescent`
- `evidence`
- current public examples
- micro-real bundles

## Worker Modes

Exercised worker modes:

- serial
- thread
- process

Agreement on the full seeded soak run was exact:

- serial vs thread
  - ordering match: `true`
  - classification mismatch count: `0`
  - max numeric difference: `0.0`
- serial vs process
  - ordering match: `true`
  - classification mismatch count: `0`
  - max numeric difference: `0.0`

Failure propagation was also identical across all three backends. A deliberate bad task at index `1` raised:

`batch task 1 failed for task={'label': 'bad_1', 'fail': True}: intentional failure for bad_1`

The forced process-blockage path fell back to the deterministic thread implementation and matched serial exactly.

## Long-Run Dashboards

The current machine-readable soak summary is:

- [stack_soak_summary.json](../tools/outputs/stack_soak_summary.json)

High-signal kernel metrics from that run:

- max hidden-rank gap: `0.0`
- max inverse round-trip residual: `1.42e-12`
- median inverse round-trip residual: `7.74e-16`
- max projector residual: `1.19e-12`
- median projector residual: `4.25e-15`
- max two-step transport residual: `2.51e-14`
- max long-chain clock residual: `3.20e-14`
- max zero-support residual: `0.0`

High-signal descent metrics:

- max tower residual: `3.06e-13`
- median tower residual: `4.00e-15`
- exact factorisation and relation outputs were stable across all worker modes
- incompatible completion cases retained their explicit linear residual certificates

Example stability metrics:

- entanglement max `tau` residual: `1.33e-15`
- Bell phase classifications were stable across reruns
- Arrow clock ordering remained stable

Kernel validation sweep highlights:

- max lift/projector residual: `7.90e-11`
- max `H`-projector symmetry residual: `4.28e-12`
- max tower-law residual: `9.20e-14`
- max inverse `Lambda -> X -> Lambda` relative error: `2.50e-15`
- max inverse `X -> Lambda -> X` relative error: `1.92e-15`
- max hidden clock residual: `6.22e-15`
- max downstairs associativity residual: `3.30e-15`
- max divergence contraction violation: `0.0`

The asymptotic/local tests remained the largest residuals, as expected:

- max curvature split residual: `4.07e-2`
- max local `O(t^3)` normalized residual: `1.35e-1`
- max local hidden-birth residual: `1.95e-2`
- max DV hidden normalized residual: `6.76e-2`

These are not theorem failures. They are the known finite-precision remainder regime of the local expansions.

## Bug And Fixes

One issue was exposed in the soak harness, not in the kernel:

- the boundary kernel task originally chose a smallest active ceiling eigenvalue too close to the support cutoff and then asked for reduced-coordinate inverse reconstruction
- this could cross the numerical support threshold and make the harness ask for the wrong reduced dimension
- fix:
  - raise the smallest active test eigenvalue from `1e-8` to `1e-4`
  - compute the hidden rank-gap dashboard entry with the hidden-load metadata tolerance rather than NumPy's default rank heuristic

No theorem-layer defect was exposed in `nomogeo`.

## Package And Install Surface

The current package-soak summary is:

- [install_soak_summary.json](../tools/outputs/install_soak_summary.json)

Completed successfully:

- editable install of root, `nomodescent`, and `evidence`
- root test suite: `48 passed`
- `nomodescent` test suite: `14 passed`
- `evidence` test suite: `18 passed`
- example runs
- `nomodescent` worked examples
- `evidence` worked examples
- micro-real bundle validation runs
- wheel builds for all three packages
- clean wheel-environment install and minimal import smoke

Operational note:

- `pip wheel` inside the sandbox failed with Windows temp-permission errors
- the wheel-build step was rerun outside the sandbox with explicit approval
- this was an environment issue, not a package defect

## Micro-Real Bundle Status

All three micro-real bundles ran cleanly:

- Bell counts bundle
  - raw exact counts -> `incompatible_by_linear_inconsistency`
  - projected no-signalling Gaussian surrogate -> `approximate_common_descent`
- Iris protocol mismatch
  - relation -> `non_nested_observers`
  - completion -> `underdetermined_affine_family`
  - refinement winner -> `full_iris`
- leaderboard benchmark slice
  - initial assembly -> `underdetermined_evidence`
  - selected relation -> `non_nested_observers`
  - completion -> `underdetermined_affine_family`

## Overall Judgement

The current stack survives hard local use within its declared domain.

What held exactly:

- multiworker determinism
- batch ordering
- explicit failure propagation
- identity-level geometry checks
- inverse theorem round-trips
- factorisation and relation classification
- packaged editable installs and installed-surface imports

What remains explicitly approximate:

- the audited low-dimensional affine PSD search in `nomodescent`
- the Gaussian-surrogate layer when real non-Gaussian evidence is deliberately projected into the stack

That boundary is now exposed in the evidence bundles and their audits rather than hidden.


