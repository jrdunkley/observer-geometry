# AP Hook Digest v0

Status: working notes for attachment candidates, not cards.

Working mode:

`AP theory digestion -> candidate physical hook -> attachment test -> card or refusal -> lesson for the next hook`

## Current Source Impressions

The AP/A-level sheets are useful because they are already organized around operational primitives:

- measured quantities: displacement, velocity, force, voltage, current, pressure, temperature, intensity
- controls: spring constant, mass, length, capacitance, inductance, resistance, charge, field strength
- assumptions: inertial frame, negligible air resistance, ideal springs/strings, ideal fluids
- equations that are either already quadratic, become quadratic after a declared operating point, or require an additional fluctuation model

## Strong Candidate Hooks

### RLC Linear Response

- Formula surface: circuit energy plus resistance and driven response.
- Candidate attachment: exact quadratic storage object plus dissipative/dynamic wrapper.
- Module object: LC energy Hessian for storage; response/frequency data as a later sensor task.
- Attachment test: pass for storage Hessian; defer damping response until a dynamic convention is declared.
- Lesson: RLC should be the next analogue after static LC, but it must not blur static energy geometry with dissipative time evolution.

### Normal Modes / Coupled Oscillators

- Formula surface: springs, masses, oscillation frequency, period.
- Candidate attachment: exact quadratic Hessian plus observer map over coordinates or modal coordinates.
- Module object: stiffness Hessian; possibly mass-weighted Hessian if dynamic frequencies are the target.
- Attachment test: pass after coordinate and mass convention are declared.
- Lesson: mass matrix matters; do not compare stiffness geometry and frequency geometry without declaring the metric.

### Small-Oscillation Pendulum Family

- Formula surface: `T = 2 pi sqrt(ell / g)` and small-angle assumptions.
- Candidate attachment: local quadratic approximation at stable equilibrium.
- Module object: local Hessian `m g ell`.
- Attachment test: pass only at declared operating point.
- Lesson: textbook period formulas often hide a local reduction step.

### Capacitor / Electric Field Energy

- Formula surface: `U = 1/2 C V^2`, `C = epsilon A / d`.
- Candidate attachment: exact quadratic in declared voltage or charge coordinate.
- Module object: Hessian differs by coordinate choice: charge coordinate gives `1/C`; voltage coordinate gives `C`.
- Attachment test: pass only after coordinate choice.
- Lesson: the same physical device can produce different-looking Hessians depending on whether charge or voltage is the state coordinate.

### Hooke / Young Modulus

- Formula surface: `F = k Delta L`, stress, strain, Young modulus.
- Candidate attachment: exact quadratic elastic energy in small-strain linear regime.
- Module object: stiffness/Hessian after geometry and strain coordinate declaration.
- Attachment test: pass in linear elastic regime.
- Lesson: good route from school mechanics into material stiffness, but geometry normalization must be explicit.

## Boundary Hooks

### Ideal Gas Law

- Formula surface: `pV = nRT`.
- Attachment test: not enough by itself.
- Missing declaration: ensemble, fluctuation model, thermodynamic potential, or local expansion.
- Lesson: famous formula is not automatically a Hessian/covariance object.

### Chemistry Equilibrium Constants

- Formula surface: `Kc`, `Kp`, `Ka`, `Kb`, `Delta G = -RT ln K`.
- Attachment test: promising but not v0.
- Missing declaration: reaction coordinate, thermodynamic potential, concentration covariance or Hessian convention.
- Lesson: strong future attachment route, but it should enter as a declared fluctuation/thermodynamic model rather than an immediate quadratic card.

### Kinetics

- Formula surface: zero/first/second-order rate laws.
- Attachment test: not a static quadratic object.
- Missing declaration: observation time series, noise model, Fisher/Hessian of parameter estimation, or local log-likelihood.
- Lesson: likely an evidence/estimation attachment rather than a direct `visible_precision` card.

## Next Candidate After Hardening

The RLC storage pair has now been added:

- `rlc_resonator_linear_storage`
- `coupled_rlc_resonators_observe_one_node`

Both cards keep resistance out of the quadratic storage object and mark damping/linewidth as practical sensor features or deferred dynamic wrappers.

`kinetic_and_rotational_energy_metrics` has now been promoted as the first metric-curvature card. `string_fixed_end_three_mode` has now been promoted as the first finite-modal card after explicitly separating mass metric, stiffness Hessian, and frequency readout.
