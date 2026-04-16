# Metric And Modal Route Plan v0

Status: scientific plan for the next attachment move, not a new theory claim.

The next exact route should be metric/modal curvature. It is closer to the module's positive-cone core than the measurement pipeline, and it addresses a current weakpoint directly: many AP formulas involving speed, angular speed, frequency, and modes are not stiffness Hessians, but they are often quadratic after a metric and coordinate convention are declared.

## Core Judgement

The existing good cards prove the storage route:

```text
potential or field energy -> stiffness/storage Hessian -> visible precision
```

The metric route is different:

```text
kinetic energy -> mass/inertia metric -> finite positive quadratic object
```

The modal route then combines the two:

```text
boundary + basis + truncation -> (M, K) -> modal frequencies and finite observer maps
```

This distinction matters because frequency formulas are not module objects by themselves. A modal frequency is an eigenvalue readout after both a stiffness matrix `K` and mass metric `M` have been declared.

## Promotions Completed

`kinetic_and_rotational_energy_metrics` was promoted first.

Reason: it is exact, finite, and forces the vocabulary distinction before wave cards depend on it. The card should license `visible_precision` only over the declared kinetic metric. It must state that this is metric geometry, not a restoring stiffness or frequency geometry.

`string_fixed_end_three_mode` was then promoted once the sign-off language said clearly:

- which coordinates are modal amplitudes;
- which object is the mass metric `M`;
- which object is the stiffness Hessian `K`;
- whether the runner is reading `K`, `M`, or a mass-whitened operator;
- why hidden-load is refused unless a ceiling/reference is declared.

`string_point_sensor_modal_observer` was then promoted as the first practical observer-map variant over the finite modal object.

`string_modal_hidden_load_with_declared_ceiling` was then promoted as the first modal reference-ceiling hidden-load variant. It does not create a canonical visible block for a point sensor; it declares a reference shape `r` and ceiling `T = r^T K r`, then compares that ceiling with the relaxed quotient precision.

The existing finite-modal proof packet is valuable and remains the staging proof. It establishes that the raw standing-wave formula becomes finite only after boundary conditions, basis, and truncation.

## Scientific Payoff

If this route holds, the operational layer gains a second exact bridge:

```text
storage curvature: spring / LC / elastic bar / capacitor
metric curvature: mass / inertia / finite modal mass metric
modal sensor route: finite modal K/M object -> point-sensor observer map -> quotient visible stiffness -> declared reference-ceiling hidden-load readout
```

That prepares the later wave/acoustic bridge without pretending that every AP formula is a storage Hessian.

## Refusal Discipline

The following remain refusals:

- raw kinetic formulas without a velocity or angular-velocity coordinate;
- comparing mass and stiffness entries as if they had the same units;
- deriving modal frequency from `K` alone;
- calling hidden load on a metric or modal object without a declared ceiling/reference;
- treating a reference-shape ceiling as a canonical visible-block ceiling for the point sensor;
- treating an infinite wave equation as a finite matrix object without boundary, basis, and truncation.
