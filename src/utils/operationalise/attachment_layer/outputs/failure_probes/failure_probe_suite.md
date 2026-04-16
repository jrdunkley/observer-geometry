# Failure Probe Suite v0

These are not attachment cards. They interrogate representative failures to determine what extra declaration would be needed before a module-facing object exists.

| probe | route | declared object after extra declarations | lesson |
| --- | --- | --- | --- |
| `rc_charge_decay` | `dynamic_evidence` | scalar Fisher precision for tau | The failure points to an evidence/Fisher attachment schema, not to a static hidden-load card. |
| `absorbance_beer_lambert` | `sensor_or_observer_map` | scalar measurement Fisher precision for concentration | Sensor laws can become module-facing evidence only after a noise model and estimated parameter are declared. |
| `nuclear_decay_poisson` | `dynamic_evidence` | scalar Poisson Fisher precision for lambda | This is a strong evidence hook; the missing object is likelihood structure, not more formula digestion. |
| `standing_wave_string_modes` | `finite_weighted_family` | finite diagonal modal stiffness family with mass metric | The obstruction is infinite/basis-free description; a finite mode declaration converts it into a clean quadratic family. |
| `gravity_inverse_square_local` | `local_quadratic` | positive local radial Hessian of an effective potential at a stable circular orbit | The failure is informative: local field attachment depends on stability and effective-potential declarations. |
| `ideal_gas_law_boundary` | `thermodynamic_fluctuation` | local scalar Hessian d2F/dV2 in declared isothermal Helmholtz model | Thermodynamic formulas attach only after ensemble and potential are declared; the refusal remains correct. |
| `electrochemistry_nernst_voltage` | `thermodynamic_fluctuation` | scalar voltage-measurement Fisher precision for log Q | Electrochemistry is a strong bridge, but it should enter through a declared thermodynamic/evidence schema. |
| `relativity_lorentz_metric` | `exact_special_sector` | none in v0 static attachment runner | This is a real boundary for the current SPD route, not merely missing units. |

## Readouts

### `rc_charge_decay`

- surface formula: `V(t) = V0 exp(-t/tau)`
- refusal: The exponential law is not a static Hessian.
- minimal extra declarations: sample times, Gaussian voltage noise sigma, fit parameter tau

```json
{
  "tau": 2.0,
  "v0": 5.0,
  "sigma": 0.05,
  "sample_count": 16,
  "fisher_tau": 2392.753487271695
}
```

### `absorbance_beer_lambert`

- surface formula: `A = epsilon * b * c`
- refusal: The absorbance equation is a sensor law, not a storage Hessian.
- minimal extra declarations: molar absorptivity, path length, Gaussian absorbance noise, concentration parameter

```json
{
  "molar_absorptivity": 1.2,
  "path_length": 1.0,
  "sigma_absorbance": 0.01,
  "fisher_concentration": 14400.0
}
```

### `nuclear_decay_poisson`

- surface formula: `N(t) = N0 exp(-lambda t)`
- refusal: The decay law is a stochastic counting process, not a static quadratic storage object.
- minimal extra declarations: counting windows, Poisson likelihood, initial rate, decay parameter lambda

```json
{
  "lambda": 0.18,
  "initial_rate": 200.0,
  "sample_count": 8,
  "fisher_lambda": 13622.170030659636
}
```

### `standing_wave_string_modes`

- surface formula: `string harmonic frequencies and wave speed`
- refusal: The wave formula is not finite until boundary conditions, basis, and truncation are declared.
- minimal extra declarations: fixed-end boundary, sine mode basis, three-mode truncation, tension, linear density

```json
{
  "modes": [
    1,
    2,
    3
  ],
  "stiffness_diagonal": [
    197.39208802178717,
    789.5683520871487,
    1776.5287921960844
  ],
  "mass_metric_diagonal": [
    0.025,
    0.025,
    0.025
  ],
  "omega_squared": [
    7895.6835208714865,
    31582.734083485946,
    71061.15168784338
  ],
  "min_stiffness": 197.39208802178717
}
```

### `gravity_inverse_square_local`

- surface formula: `U(r) = -mu/r`
- refusal: The raw inverse-square potential is not globally positive quadratic; local radial curvature has the wrong sign in this simple coordinate.
- minimal extra declarations: operating radius, effective potential, angular momentum/constraint convention, local stability sector

```json
{
  "mu": 1.0,
  "r0": 2.0,
  "raw_potential_second_derivative": -0.25,
  "effective_potential_second_derivative": 0.125
}
```

### `ideal_gas_law_boundary`

- surface formula: `p V = n R T`
- refusal: The equation of state alone is not a module object.
- minimal extra declarations: fixed n,T ensemble, Helmholtz free-energy model, volume coordinate, operating point

```json
{
  "n_moles": 1.0,
  "temperature_K": 300.0,
  "volume_m3": 0.024,
  "helmholtz_volume_hessian": 4330449.280208333
}
```

### `electrochemistry_nernst_voltage`

- surface formula: `E = E0 - (RT/nF) ln Q`
- refusal: Voltage is measurable, but reaction coordinate and thermodynamic convention are missing.
- minimal extra declarations: reaction quotient coordinate log Q, electron count, temperature, voltage noise

```json
{
  "temperature_K": 298.15,
  "electrons": 1.0,
  "sigma_voltage": 0.002,
  "dE_dlogQ": -0.025692579121493725,
  "fisher_logQ": 165.02715547855382
}
```

### `relativity_lorentz_metric`

- surface formula: `Minkowski metric diag(-1,1,1,1)`
- refusal: The current static kernel calls require SPD objects; this metric is indefinite.
- minimal extra declarations: positive subproblem or special sector, signature convention, observable map

```json
{
  "metric_diagonal": [
    -1.0,
    1.0,
    1.0,
    1.0
  ],
  "eigenvalues": [
    -1.0,
    1.0,
    1.0,
    1.0
  ],
  "min_eigenvalue": -1.0
}
```
