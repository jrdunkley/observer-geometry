# Two-Track Operational Programme v0

Status: working memo for the attachment layer.

The current attachment layer has two different jobs:

1. understand why some AP/A-level physical hooks do not integrate cleanly into the first static card runner;
2. push the hooks that do integrate cleanly to a sign-off standard strong enough to call them good.

The two tracks should run in parallel. Do not broaden the runnable card set just because a formula is famous. Do not leave the clean quadratic cards at demo quality when they are the first real bridge to physical use.

## Track A: Diagnose Non-Integrating Hooks

A non-integrating hook is not automatically a failure. It may be:

- misrouted to the wrong module surface;
- missing an attachment declaration;
- outside the current v0 runner but inside the broader module stack;
- a real theoretical or implementation gap.

The current runner mostly exercises `visible_precision` and `hidden_load`. The module stack is broader: local visible calculus, local quadratic ensembles, weighted-family frontiers, evidence assembly, quotient descent, support-stratum transport, affine-hidden reduction, covariance/Fisher perturbation, and residual-margin diagnostics all exist as possible homes.

## Diagnosis Classes

| Class | Meaning | Likely Module Route | Example Hooks |
| --- | --- | --- | --- |
| Static quadratic | Already a finite Hessian/precision after coordinate declaration | `visible_precision`, `hidden_load` | springs, LC storage, elastic bar, capacitor energy |
| Local quadratic | Becomes quadratic only at an operating point or regime | `visible_precision`, `local_visible_calculus`, local ensembles | small-angle pendulum, local gravity/Coulomb expansion |
| Dynamic law | Time evolution rather than static energy geometry | evidence/Fisher layer, state-space wrapper, local calculus only after derivative data | RC decay, kinetics, nuclear decay |
| Sensor law | Measurement map rather than module object | evidence assembly, observer-map grammar, Fisher after noise model | absorbance, lens equations, diffraction readouts |
| Finite family | Requires several declared quadratics/observers | weighted-family frontier, declared ladder, quotient descent | normal modes, observer comparisons, wave-mode families |
| Thermodynamic/fluctuation | Needs ensemble, potential, coordinates, and local expansion | local quadratic ensemble, evidence, declared fluctuation model | ideal gas, equilibrium constants, heat capacity |
| Support/event | About rank changes, branch changes, or near-zero modes | support-stratum transport, kernel Schur jets, residual margins | buckling-like transitions, support birth/death, threshold events |
| Exact special sector | Not generic quadratic v0, but may have a supplied exact law sector | affine-hidden, covariance/Fisher perturbation, exact special sectors | quantum oscillator, selected non-Gaussian affine-hidden cases |
| True gap | No honest current route without new theory/API | record as gap; do not card as if solved | arbitrary nonlinear full-law systems, global noncommuting optimisers |

## Immediate Diagnostic Priorities

1. `rc_charge_decay`
   - Current issue: dynamic law, not static Hessian.
   - Likely route: time-series evidence/Fisher card, or state-space wrapper.
   - Gap question: do we want a lightweight dynamic attachment schema, or should this wait for evidence integration?

2. `chemistry_kinetics_rate_laws`
   - Current issue: parameter-estimation law over concentration time series.
   - Likely route: evidence/Fisher layer.
   - Gap question: how do we declare observation noise and sampling without overbuilding?

3. `electrochemistry_nernst_voltage`
   - Current issue: voltage is a strong sensor surface, but free-energy coordinate is undeclared.
   - Likely route: declared thermodynamic potential plus evidence/Fisher voltage readout.
   - Gap question: what is the smallest reaction-coordinate schema that is not fake?

4. `ideal_gas_law_boundary`
   - Current issue: equation of state alone is not a module object.
   - Likely route: refusal first; later declared fluctuation or thermodynamic Hessian model.
   - Gap question: which ensemble/potential should be the first honest gas-law attachment?

5. `acoustic_ultrasound_impedance` and non-reference hidden-load string variants
   - Current issue: acoustic fields still need finite variables, boundary, observer map, and sensor convention; the first string hidden-load card now exists only as a declared reference-ceiling variant.
   - Likely route: finite modal family, then observer comparison or sensor/evidence route.
   - Gap question: what is the smallest acoustic modal observer or non-reference string ceiling/reference that remains physically legible without claiming hidden load prematurely?

6. `paraxial_optics_lens` and `interference_diffraction`
   - Current issue: observer/sensor maps, not SPD objects by themselves.
   - Likely route: evidence or observer-map task mode.
   - Gap question: should optics enter through ray-transfer observer maps or through wave-field/sensor likelihoods?

7. `relativity_lorentz_metric`
   - Current issue: indefinite metric and frame geometry do not belong in current SPD calls.
   - Likely route: exact special sector only after a positive declared subproblem is isolated.
   - Gap question: is this strategically useful soon, or a distraction?

## Track B: Sign Off The Integrating Hooks

The clean objects should not remain loose demos. A signed-off good card should satisfy a stricter standard:

1. The physical regime is declared and not vague.
2. The state coordinates are explicit, with units.
3. The sensor surface is split into ideal and practical.
4. The control variables are explicit and actually controllable.
5. The equation hook derives the declared matrix, not just inspires it.
6. The module object is named precisely: Hessian, precision, covariance, weighted family, or refusal.
7. The derivation status is exact/local/declared/inferred with no inflation.
8. The allowed module calls are minimal and enforced by the runner.
9. Hidden-load cards declare a ceiling/reference and validate `T >= Phi` before calling `hidden_load`.
10. The output exposes the human interpretation: `T`, `Phi`, `T - Phi`, hidden load, and clock where applicable.
11. The positive scope and boundary are both present and physically meaningful.
12. There is at least one practical measurement or control move that would stress the card.
13. There is a sibling analogue or a declared reason why no analogue exists yet.
14. The card runs from default numerical parameters.
15. The card survives validation and a negative/refusal test where relevant.

## Good v0 Candidates

| Card | Current Sign-Off Status | Remaining Work |
| --- | --- | --- |
| `mass_spring_single` | near good | add explicit experimental check and maybe a units sanity assertion |
| `lc_resonator_single` | near good | preserve coordinate warning; add charge/voltage sibling later |
| `mass_spring_coupled_observe_one` | good v0 bridge | coupling sweep now shows monotone ceiling gap, hidden load, and clock |
| `coupled_lc_resonators_observe_one_node` | good v0 bridge | coupling sweep now shows the same monotone hidden-renormalisation story as coupled springs |
| `pendulum_small_angle` | good as regime card | add amplitude boundary note or refusal-pair link |
| `pendulum_full` | good as refusal card | keep refusal strict; do not sneak in local calls |

## What This Means Broadly

The project is not discovering that all physics is already a nomogeo object. It is discovering that many physical hooks fall into a small set of attachment routes:

```text
stored quadratic -> static observer geometry
local approximation -> local observer geometry
time series -> evidence/Fisher route
sensor law -> observer/evidence route
finite modes -> weighted family/frontier route
thermodynamics -> declared fluctuation route
support transition -> support-stratum route
special law sector -> exact special-sector route
```

That is a better result than indiscriminate integration. It means the attachment layer is beginning to classify the difficulty of physical contact.

## Current Superintelligent Summation

The first true physical bridgehead is cross-substrate quadratic storage under partial observation. Springs and LC resonators are not just examples; they are the first evidence that different material systems can present the same observer-facing geometry.

The non-integrating objects are not noise. They are a map of the next schemas the operational layer needs: dynamic evidence, sensor maps, finite mode families, thermodynamic fluctuations, and support events.

The immediate discipline should be:

- keep the signed-off spring/LC/pendulum v0 objects stable;
- add the elastic-bar/capacitor/RLC analogue cards only after that;
- reclassify the non-integrating backlog by module route before calling anything a theoretical gap.
