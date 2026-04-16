# Local Positive Curvature Admissibility v0

Status: theoretical operational memo, not a theorem added to `nomogeo`.

Director correction: do not jump directly to an Evidence Attachment schema. First extract the idea.

The idea is:

> A physical surface becomes admissible for observer-geometry attachment only after it yields a declared finite local positive curvature object.

The source of curvature may differ. It can come from stored energy, local approximation, sensor likelihood, thermodynamic free energy, a finite mode truncation, or a supplied exact special sector. But the module-facing object must be explicit and positive on the active support before theorem-local calls are licensed.

## Why This Matters

The working cards showed cross-substrate quadratic storage:

```text
spring energy -> stiffness Hessian
LC energy -> storage Hessian
coupled spring -> visible precision and hidden load
coupled LC -> visible precision and hidden load
```

The failure probes suggest a wider but still disciplined class:

```text
declared local positive curvature
```

That includes energy curvature, information curvature, free-energy curvature, finite modal curvature, and possibly special-sector curvature. It does not include bare formulas.

## Admissibility Test

A physical hook is admissible only if all of the following are declared:

1. `surface`
   - the physical law, formula, energy, sensor relation, likelihood, or potential being used.

2. `coordinate`
   - the finite coordinate or parameter in which curvature is being taken.

3. `locality`
   - the point, regime, support, basis, sample grid, or truncation that makes the object finite and local.

4. `curvature_source`
   - one of: energy, likelihood, free energy, finite modal truncation, local effective potential, exact special sector.

5. `positive_object`
   - a finite SPD/PSD Hessian, precision, Fisher, covariance inverse, or declared positive sector.

6. `units_and_scale`
   - units and coordinate scaling, especially where cross-substrate comparison is tempting.

7. `observer_or_sensor_surface`
   - the physical measurement channel or observer map that exposes the object.

8. `allowed_module_route`
   - which module surfaces are allowed, and which remain refused.

9. `boundary`
   - what is still not captured.

If any one of these is missing, the correct result is refusal or deferment.

## Curvature Sources

| Source | Example | What Becomes Positive | Main Risk |
| --- | --- | --- | --- |
| Stored energy | spring, LC, capacitor, elastic bar | energy Hessian | confusing coordinate units or dynamic response with storage |
| Local approximation | small-angle pendulum, local field patch | local Hessian | global overclaim |
| Likelihood / evidence | RC decay, Beer-Lambert, nuclear decay | Fisher precision | hiding noise/sampling assumptions |
| Thermodynamic free energy | ideal gas under declared ensemble, Nernst voltage with reaction coordinate | free-energy Hessian or Fisher | ensemble/potential ambiguity |
| Finite modal truncation | string modes, acoustic modes | finite modal stiffness/mass metric | infinite field object smuggled in without basis/truncation |
| Exact special sector | quantum oscillator, affine-hidden law sector | declared sector curvature | pretending special-sector exactness is generic |
| Indefinite metric | Lorentzian relativity | no SPD object without positive sector | forcing SPD calls on an indefinite object |

## Core Principle

The attachment unit is not:

```text
physical system
```

and not:

```text
formula
```

It is:

```text
declared local positive curvature object with observer/sensor route
```

This is the conceptual bridge from the successful storage cards to the failure cases.

## Consequence For Cross-Substrate Claims

The strongest defensible claim is not that different sciences are the same.

The stronger and more disciplined claim is:

> Different substrates can instantiate the same observer-facing local positive curvature role after their attachment declarations are made explicit.

The spring/LC pair proves this for stored energy curvature in miniature.

RC decay and Beer-Lambert test whether the same discipline works for information curvature.

## Why RC Decay And Beer-Lambert Are The Right Tests

They are deliberately simple and ontologically different.

`rc_charge_decay`:

- source: dynamic physical relaxation law;
- coordinate: parameter `tau`;
- locality: declared finite sample times;
- curvature source: Gaussian likelihood;
- positive object: scalar Fisher precision for `tau`;
- sensor surface: voltage over time.

`absorbance_beer_lambert`:

- source: static sensor law;
- coordinate: concentration `c`;
- locality: declared measurement protocol;
- curvature source: Gaussian measurement likelihood;
- positive object: scalar Fisher precision for `c`;
- sensor surface: absorbance.

If both pass, the result is not merely an evidence schema. The result is a broader admissibility principle:

```text
stored energy curvature and information curvature can both be admitted,
but only after the local positive object and its route are declared.
```

## What Would Count As Failure

The test fails if:

- the positive object depends on undeclared noise;
- the coordinate is not physically interpretable;
- the Fisher object is singular without acknowledging support;
- the readout is treated as hidden load without a ceiling/reference;
- the result is described as intrinsic to the formula alone;
- the dynamic/sensor route is collapsed into the static storage route.

## Current Recommendation

Do not build a general Evidence Attachment schema yet.

First write two proof-of-admissibility examples:

1. `rc_decay_local_curvature`
2. `beer_lambert_local_curvature`

Each should explicitly show:

- refusal of the bare formula;
- declarations required for admissibility;
- derivation of the scalar Fisher curvature;
- positivity check;
- allowed route as `local_positive_curvature`, not `hidden_load`;
- boundary statement.

Only after those two work should the schema be promoted.

## Stress Test Requirements

The admissibility examples must also survive basic invariance and refusal checks before schema promotion:

- reparameterization covariance: Fisher curvature must transform as `I_new = J^T I_old J`;
- noise monotonicity: increasing declared sensor noise should weaken Fisher curvature;
- zero sensitivity: a sensor with zero derivative should produce zero curvature, not a false positive object;
- missing noise refusal: no likelihood curvature is admitted without a positive noise/likelihood scale;
- route refusal: information-curvature examples must still refuse `hidden_load` without a declared ceiling and must not be described as static storage Hessians.

The generated checks live in:

- `outputs/local_curvature/local_curvature_stress_tests.json`
- `outputs/local_curvature/local_curvature_stress_tests.md`
