# Core Formula-Sheet Coverage Map v0

Status: breadth pass over the supplied AP Physics, A-level Physics, AP Chemistry, and broader physics formulary PDFs.

This is not a card set. It is a routing map: each core topic gets an attachment test and a status. A topic can be "approached" by becoming a runnable card, a candidate card, or an explicit refusal/defer note.

## Status Key

- `carded`: already represented in v0 runnable cards.
- `next`: strong next card candidate after v0 hardening.
- `candidate`: promising, but needs a declared convention first.
- `defer`: do not card until a richer layer exists.
- `refusal`: formula-sheet surface is not a module object by itself.

## Mechanics

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Hooke spring | `F = k Delta L`, `U = 1/2 kx^2` | Exact Hessian after displacement coordinate is declared | `carded` | `mass_spring_single` is the cleanest stiffness-to-Hessian card. |
| Coupled springs | spring energy with coupling | Exact Hessian plus observer map and visible-block ceiling | `carded` | `mass_spring_coupled_observe_one` is the first hidden-renormalisation card. |
| Kinematics | `v = u + at`, `s = ut + 1/2 at^2` | State evolution, not a static Hessian by itself | `defer` | Could become a trajectory/estimation card with a noise/Fisher model. |
| Work and kinetic energy | `W = F s cos theta`, `K = 1/2 mv^2` | Kinetic part is quadratic in velocity after coordinate/metric declaration | `carded` | `kinetic_and_rotational_energy_metrics` admits the mass metric while refusing to treat it as stiffness. |
| Circular motion | `a = v^2/r`, `F = mv^2/r` | Constraint dynamics, not immediately an SPD Hessian | `candidate` | Local radial stiffness cards are possible only after a potential or constraint model is declared. |
| Rotational dynamics | `E_k = 1/2 I omega^2`, `tau = I alpha` | Exact quadratic in angular velocity with declared inertia metric | `carded` | Covered by `kinetic_and_rotational_energy_metrics` after axis and angular-velocity coordinate declaration. |
| Pendulum small angle | `T = 2pi sqrt(ell/g)`, local expansion of `m g ell(1 - cos theta)` | Exact for declared local quadratic approximation at operating point | `carded` | `pendulum_small_angle`. |
| Full pendulum | `U = m g ell(1 - cos theta)` | Not one global quadratic object | `carded` as refusal | `pendulum_full` teaches boundary discipline. |
| Collisions/momentum | `p = mv`, impulse | Conservation/event relation, not a direct quadratic module object | `defer` | Could become evidence/constraint cards later. |

## Materials And Elasticity

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Young modulus | stress/strain, `F = k Delta L`, stored energy | Exact local quadratic elastic energy after geometry normalization | `carded` | `linear_elastic_bar` declares extension coordinate, area, length, and Young modulus. |
| Density | `rho = m/V` | Parameter, not module object | `candidate` | Useful as mass matrix input for wave/oscillation cards. |
| Tensile energy | `E = 1/2 F Delta L` | Quadratic after linear-elastic coordinate declaration | `carded` | Covered by `linear_elastic_bar` under the small-strain uniform axial regime. |

## Circuits And Electronics

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Single LC storage | `U = q^2/(2C) + phi^2/(2L)` | Exact Hessian in declared charge/flux coordinates | `carded` | `lc_resonator_single`; includes mixed-coordinate note. |
| Coupled LC storage | coupled flux or charge energy | Exact Hessian plus observer map and ceiling | `carded` | `coupled_lc_resonators_observe_one_node`. |
| RC charge decay | `Q = Q0 exp(-t/RC)` | Dynamic first-order system, not static Hessian | `candidate` | Good live hook later; needs time-series/noise or state-space convention. |
| RLC resonance | `f0 = 1/(2pi sqrt(LC))`, Q-factor | Storage Hessian plus dissipative dynamic wrapper | `carded` | `rlc_resonator_linear_storage` and `coupled_rlc_resonators_observe_one_node` keep resistance outside static `H`. |
| Resistance/Ohm law | `V = IR`, `P = IV` | Constitutive/dissipative relation | `candidate` | Could attach through Fisher/estimation or response, not hidden load alone. |
| Op-amp formulas | gain equations | Linear observer/control map, not SPD geometry by itself | `candidate` | Promising for sensor surfaces and observer maps, but not v0. |

## Waves, Acoustics, And Optics

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Harmonic oscillator | sinusoidal solution, frequency | Exact only after storage energy and coordinate metric are declared | `carded` via spring/LC | The equation of motion alone is not the object; the energy/Hessian is. |
| Standing waves/string | wave speed, harmonics | Normal-mode quadratic family after boundary, basis, truncation, observer declaration, and optional reference ceiling | `carded` | `string_fixed_end_three_mode` admits the first three fixed-end modes; `string_point_sensor_modal_observer` adds a point-sensor observer map; `string_modal_hidden_load_with_declared_ceiling` adds a reference-ceiling hidden-load readout. |
| Sound/acoustic waves | pressure wave speed, intensity | Finite closed-tube acoustic modal stiffness with declared pressure observer; raw impedance/intensity/reflection still refused | `carded_narrow` | `acoustic_tube_three_mode_pressure_observer` gives the narrow positive route; live microphone/speaker data remain later work. |
| Lenses/paraxial optics | lens equations, magnification | Linear ray-transfer map, not SPD precision by itself | `defer` | Could become observer-map/control card, but needs a separate optical geometry convention. |
| Interference/diffraction | fringe spacing, grating equation | Measurement equation rather than Hessian | `defer` | Requires sensor/noise model or wave-field quadratic energy. |
| Intensity/absorption | inverse-square, exponential attenuation | Observation/sensor law | `candidate` | Potential evidence encoding hook, not immediate visible precision. |

## Fields And Gravitation

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Gravitational inverse-square | `F = Gm1m2/r^2`, potential | Local Hessian only after operating point and coordinate split | `candidate` | Useful local quadratic approximation; global potential is not SPD everywhere. |
| Electric inverse-square | Coulomb force/potential | Same as gravitation, plus charge sign issues | `candidate` | Local Hessian or fluctuation model needed. |
| Capacitor energy | `U = 1/2 QV = Q^2/(2C) = 1/2 C V^2` | Exact quadratic after charge or voltage coordinate declaration | `carded` | `capacitor_charge_coordinate` cards the charge coordinate and warns that voltage coordinate is a distinct Hessian. |
| Magnetic Lorentz force | `F = qvB`, flux | Force law/dynamic constraint, not static Hessian | `candidate` | Trapped/linearized systems may attach later. |
| Induction/transformer | flux linkage, induced emf | Dynamic/observer relation | `candidate` | Strong sensor/control surface, not v0 static kernel. |

## Thermal Physics And Statistical Mechanics

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Heat capacity | `Q = mc Delta T` | Needs thermodynamic potential or fluctuation model | `candidate` | Can become local Hessian in entropy/free-energy coordinates. |
| Ideal gas | `pV = nRT` | Not enough by itself | `carded` as refusal | `ideal_gas_law_boundary` refuses theorem-local calls until ensemble, fluctuation model, local expansion, or measurement route is declared. |
| Kinetic theory | `pV = NkT`, rms speed | Distributional/statistical object | `candidate` | Potential covariance/Fisher route after variable declaration. |
| Thermodynamic potentials | `G = H - TS`, `Delta G = -RT ln K` | Hessian possible after potential and coordinates declared | `candidate` | Good future chemistry bridge. |
| Heat engines | efficiency formulas | Constraint/energy accounting, not direct module object | `defer` | Needs process model and state variables. |

## Chemistry

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Atomic photon relations | `E = hf = hc/lambda` | Measurement conversion, not quadratic object | `candidate` | Useful sensor/energy hook; not v0. |
| Gases/liquids/solutions | gas-law ratios, concentration | Needs fluctuation or thermodynamic model | `candidate` | Same boundary as ideal gas. |
| Kinetics | zero/first/second-order rate laws | Time-series parameter model, not static Hessian | `candidate` | Likely Fisher/log-likelihood card later. |
| Equilibrium constants | `Kc`, `Kp`, `Ka`, `Kb` | Thermodynamic potential bridge after reaction coordinate declaration | `candidate` | Strong, but needs careful chemistry schema. |
| Electrochemistry | `Delta G = -n F E`, Nernst form | Potential/chemical-coordinate bridge | `candidate` | Good future attachment to measured voltage sensors. |
| Calorimetry | `q = mc Delta T` | Needs uncertainty/fluctuation model for module contact | `candidate` | Useful lab hook, but not direct hidden-load card. |

## Modern Physics, Nuclear, Astrophysics

| Topic | Formula-Sheet Surface | Attachment Test | Status | Notes |
| --- | --- | --- | --- | --- |
| Photoelectric effect | `hf = phi + Kmax` | Linear energy accounting, not SPD object | `candidate` | Could attach via measurement/Fisher model. |
| de Broglie relation | `lambda = h/p` | Coordinate transformation/measurement equation | `candidate` | Not a quadratic object by itself. |
| Nuclear decay | `N = N0 exp(-lambda t)` | Time-series/Poisson likelihood object | `candidate` | Strong evidence/Fisher candidate. |
| Nuclear radius | `R = R0 A^(1/3)` | Scaling relation | `defer` | Needs model fitting context. |
| Blackbody/Wien/Stefan | radiation laws | Distribution/field model | `candidate` | Good future fluctuation/fit layer, not v0. |
| Hubble law | `v = H d` | Linear regression/observer evidence object | `candidate` | Candidate for evidence layer, not immediate Hessian. |
| Relativity | Lorentz transforms, metric | Metric/coordinate geometry, not current SPD observer kernel by default | `defer` | Only local Euclideanized/quadratic sectors should attach. |
| Quantum oscillator | harmonic oscillator Hamiltonian | Exact quadratic operator/phase-space sector after convention | `candidate` | Promising but not AP v0. |

## Strong Next Cards

The clean direct storage additions, `linear_elastic_bar`, `capacitor_charge_coordinate`, `rlc_resonator_linear_storage`, and `coupled_rlc_resonators_observe_one_node`, have been promoted. `kinetic_and_rotational_energy_metrics` has been promoted as the first metric-curvature card. `string_fixed_end_three_mode` has been promoted as the first finite-modal card, `string_point_sensor_modal_observer` adds the first practical modal observer map, and `string_modal_hidden_load_with_declared_ceiling` adds the first modal reference-ceiling hidden-load readout. `ideal_gas_law_boundary` has also been promoted as a refusal card.

The remaining breadth-preserving but disciplined next moves should now be:

1. `rc_charge_decay_fisher_tau`
   - admitted as `rc_decay_q7_probe`; the next improvement is a cleaner apparatus dataset with declared resistance/noise if capacitance is to be inferred.

2. `absorbance_beer_lambert_fisher_concentration`
   - admitted as `beer_lambert_uvvis_probe`; the next improvement is declared path length and instrument/replicate uncertainty.

3. `acoustic_modal_observer_with_declared_boundary`
   - now represented narrowly by `acoustic_tube_three_mode_pressure_observer`;
   - raw impedance, intensity, attenuation, and reflection still need additional dynamic/sensor declarations.

4. `thermistor_calibration_analysis_function`
   - still a candidate analysis-function route with calibration range and uncertainty propagation.

5. thermodynamic fluctuation route
   - still the next deep theoretical route; do not promote raw gas laws or equilibrium constants without ensemble/potential declarations.
