# Source Exhaustion Log v0

Status: deletion-safe routing log for the supplied physics and chemistry PDFs.

This file does not preserve the PDFs. It preserves the operational result of reading them: which core formula-sheet families have been approached, what attachment test they face, and where they sit in the attachment backlog.

## Sources Read

| Source | Pages | Role In This Pass | Coverage Captured |
| --- | ---: | --- | --- |
| `ap-physics-1-equations-sheet.pdf` | 3 | 2026 AP Physics 1 reference sheet | Conventions, mechanics and fluids variables, geometry, ideal spring/string/fluid assumptions. |
| `20200327 ADVANCED PLACEMENT PHYSICS 1 EQUATIONS.pdf` | 5 | Older AP Physics 1 equation sheet with usage notes | Mechanics, electricity, waves, geometry; useful explicit reminders that equations are usage surfaces, not module objects. |
| `physics_equation_tables.pdf` | 6 | 2008/2009 AP Physics B/C equation tables | Mechanics, E&M, fluids/thermal, waves/optics, atomic/nuclear, calculus Physics C extensions. |
| `Data Sheet.pdf` | 8 | AQA A-level physics formula/data sheet | Materials, electricity, fields/capacitors, SHM, thermal physics, magnetic fields, nuclear, astrophysics, medical physics, engineering physics, electronics. |
| `ap-chemistry-equations-sheet.pdf` | 4 | 2026 AP Chemistry reference sheet | Atomic photon/Coulomb hooks, gases/solutions, kinetics, equilibrium, thermodynamics, electrochemistry, calorimetry. |
| `physics-formulary.pdf` | 108 | Broad advanced physics formulary | Breadth check across mechanics, Hamilton/Lagrange mechanics, E&M, relativity, oscillations, waves, optics, statistical physics, thermodynamics, transport, quantum, plasma, solid state, nuclear, QFT/particle, astrophysics. |
| `hypothetical_objectives.txt` | n/a | Director objective note | Cross-substrate and shared-geometry orientation; used as a filter for analogue-pair priority. |
| `JCGM_GUM_6_2020.pdf` | 103 | Measurement model and uncertainty source | Measurand declaration, measurement equations, observation equations, statistical models, local linearisation, Monte Carlo adequacy checks, random variation, shared effects, calibration/analysis functions, process-model adequacy. |
| `NIST.TN.1900.pdf` | 105 | Practical uncertainty and inference source | Measurement equation versus observation equation workflow, uncertainty propagation, Monte Carlo, thermistor calibration, cadmium calibration, random effects, likelihood/Bayesian examples, model adequacy. |
| `mpc.pdf` | 426 | Measurement process characterization source | Check standards, bias, long-term variability, repeatability/reproducibility, calibration designs, gauge R&R, uncertainty budgets, control charts, process validity gates. |
| `pmd.pdf` | 60 | Process modeling source | Deterministic function plus random component, design principles, least squares, weighted least squares, model selection/fitting/validation, residual diagnostics, calibration use cases. |
| `5-Eigenvalue and eigenvector.ipynb` | n/a | Finite-mode/eigenbasis source | Diagonalisation, eigenvalue stability, matrix powers, matrix exponentials, Hermitian/symmetric structure, similarity transforms; used to discipline finite modal route. |

## Exhaustion Standard

For this pass, a source item is considered approached when it has one of these outcomes:

- `carded`: already represented by a runnable/refusal JSON card.
- `next`: strong next card candidate after v0 hardening.
- `candidate`: promising, but missing a declaration such as coordinate, metric, operating point, ensemble, noise model, or dynamic wrapper.
- `defer`: outside the narrow static attachment grammar, but recorded for later layers.
- `refusal`: the formula-sheet surface is explicitly not a module object by itself.

The current breadth map is `core_coverage_map.md`. The machine-readable backlog is `candidate_backlog.json`.

The newer measurement/inference digest is `measurement_inference_source_digest.md`; its machine-readable non-promoted backlog is `measurement_route_backlog.json`.

The first manual-dataset intake protocol is `empirical_dataset_protocol_v0.md`, with templates under `dataset_templates/`. It currently targets RC decay traces and Beer-Lambert calibration tables because those are the smallest useful tests of observation-equation Fisher curvature and measurement-equation local linearisation.

## Core Families Approached

### Mechanics

Hooke springs, coupled springs, small-angle pendulum, full pendulum, kinematics, kinetic energy, rotational inertia, circular motion, collisions/impulse, gravitation, work/power, and density have all been routed.

The strongest direct attachments are exact or local quadratic energy objects. Kinetic and rotational energy now add a metric-curvature attachment, which is exact but not stiffness. The remaining non-quadratic mechanics formulas are mostly dynamics, constraints, conservation laws, or local-linearization candidates.

### Materials

Hooke/Young modulus, tensile energy, density, and geometry-normalized stiffness have been routed.

`linear_elastic_bar` is now carded: it is the same story as the spring card, but with area, length, extension coordinate, and Young modulus declared explicitly.

### Circuits And Electronics

Single LC, coupled LC, capacitor energy, RC decay, RLC resonance, resistance/Ohm law, Kirchhoff circuits, transformers/induction, and op-amp gain formulas have been routed.

The clean direct objects are storage quadratics. Resistance, damping, and op-amp gain are operationally useful, but belong first as sensor/control or dynamic wrappers around declared storage objects.

### Waves, Acoustics, And Optics

Sinusoidal oscillator formulas, standing waves, wave speed, string modes, acoustic/ultrasound impedance, interference, diffraction, paraxial lens formulas, attenuation, and intensity laws have been routed.

The main lesson is that wave formulas become module objects only after boundary conditions, a basis, and a finite modal truncation are declared. The first fixed-end string modal card is now promoted. Optics formulas are mostly observer maps or sensor laws, not static SPD objects by themselves.

### Fields

Newtonian gravity, Coulomb force, electric field/potential, capacitor fields, Lorentz force, magnetic flux, induction, and transformer equations have been routed.

Inverse-square forces are local candidates, not global SPD objects. Capacitor storage is the clean field-to-circuit bridge because the quadratic energy is exact after a coordinate choice.

### Thermal And Statistical Physics

Heat capacity, ideal gas, kinetic theory, thermodynamic potentials, heat engines, blackbody laws, and statistical-distribution formulas have been routed.

The direct lesson is refusal discipline: equations of state and process formulas need a thermodynamic potential, ensemble, fluctuation model, or local expansion before the module can touch them.

### Chemistry

Photon relations, Coulombic atomic structure, gas laws, concentration/absorbance, kinetics, equilibrium constants, acid/base formulas, calorimetry, Gibbs relations, Nernst/electrochemistry, and current/charge relations have been routed.

The most promising future chemistry hook is electrochemistry, because measured voltage gives a practical sensor surface. It still needs a reaction coordinate and thermodynamic convention before becoming a valid attachment card.

### Modern, Nuclear, Astrophysics

Photoelectric effect, de Broglie relation, nuclear decay, nuclear radius scaling, mass-energy, Wien/Stefan laws, redshift, Hubble law, relativity, and quantum oscillator content have been routed.

Decay laws and Hubble law look like evidence/Fisher cards, not static hidden-load cards. Quantum oscillator and local field modes are promising special-sector candidates after the classical attachment vocabulary is stable.

## High-Value Analogue Pairs

| Pair | Shared Story | Status |
| --- | --- | --- |
| Mass-spring oscillator / LC resonator | Exact quadratic storage in different substrates | `carded` |
| Coupled springs / coupled LC resonators | Hidden coordinate changes visible precision under one-coordinate observation | `carded` |
| Elastic bar / capacitor energy | Geometry-normalized material parameter or field storage becomes a finite Hessian after coordinate choice | `carded` |
| Damped oscillator / RLC response | Storage Hessian plus dissipative dynamic wrapper | `carded` for storage only; dynamic response deferred |
| String standing wave / acoustic mode | Boundary-conditioned modal quadratic family plus observer map and declared reference ceiling | string modal basis, point-sensor observer, reference-ceiling hidden-load, and closed-tube acoustic pressure observer carded; raw acoustic impedance/reflection still refused |
| Nuclear decay / chemistry kinetics | Time-series likelihood and parameter precision rather than static energy Hessian | `candidate` for evidence layer |
| Electrochemistry / capacitor voltage sensing | Measured voltage as practical sensor surface over stored/free-energy coordinates | `candidate` |

## Current Narrow Next Move

The existing runner and validator now remain boring with the clean direct storage, metric, finite-modal, point-sensor modal, reference-ceiling modal, acoustic pressure-observer, process-gate, and empirical information-curvature additions promoted. The next breadth-preserving additions should not repeat those.

1. `rc_charge_decay_fisher_tau`

   Admitted as `rc_decay_q7_probe`; next useful work is a cleaner apparatus dataset with declared resistance/noise if capacitance is to be inferred.
2. `absorbance_beer_lambert_fisher_concentration`

   Admitted as `beer_lambert_uvvis_probe`; next useful work is declared path length and instrument/replicate uncertainty.
3. `acoustic_modal_observer_with_declared_boundary`

   Now partially complete as `acoustic_tube_three_mode_pressure_observer`. The next acoustic move is not raw impedance or ultrasound imaging; it is apparatus-level microphone calibration, damping/drive response, or a frequency-response dataset.
4. `thermistor_calibration_analysis_function`

   Still a candidate analysis-function route.
5. thermodynamic fluctuation route

   Still the next deep theoretical route; do not promote raw gas laws or equilibrium constants without ensemble/potential declarations.

This keeps the AP material active while preserving the director's rule: attachment test first, card or refusal second.
