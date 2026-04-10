# Release Scope 0.30

This note freezes the incremental exact public claim for `nomogeo 0.30.0`.

## Exact Additions Since 0.25.0

The `0.30.0` release adds a narrow real-valued observation-field layer:

- precision-side hidden-load conversion `Pi <-> Lambda`
- support-stable transport right-hand sides for `Pi`, `Lambda`, and `tau`
- support-stratum diagnostics that separate generator existence from PSD-valid positive transport
- forced finite restart maps for support birth and support death in explicit nested support bases
- finite-order kernel Schur-complement event-jet classification
- semisimple event-block pole, clock, and desingularisation diagnostics
- explicit-derivative local coupled birth extraction for `A_cpl`
- sampled interval-family leakage, stationarity, closure, and Hessian diagnostics

## Explicitly Out Of Scope

The `0.30.0` release does not claim:

- a global smooth hidden-load PDE across support changes
- a public stratified integrator
- a public noncommuting observer optimiser
- sampled branch-continuation as a solver
- automatic finite-difference extraction of path derivatives
- replacement of `A_cpl` by anchored `Q`
- full moving-observer connection geometry
- Hermitian packet / complex positive-cone APIs
- raw single-amplitude field ontology

## Release Reading Rule

Read `0.30.0` as:

- exact support-stratified field diagnostics in finite real symmetric coordinates
- explicit local coupled extraction when derivative data are supplied
- exact event-jet diagnostics for supplied Taylor coefficients
- sampled interval-family diagnostics, not optimisation

Do not read `0.30.0` as:

- a complete field simulator
- an autonomous event-detecting integrator
- a general observer-design optimiser
- a complex wave/coherence bridge
