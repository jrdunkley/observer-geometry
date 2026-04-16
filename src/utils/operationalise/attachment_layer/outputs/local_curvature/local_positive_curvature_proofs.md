# Local Positive Curvature Proofs v0

These examples test the admissibility idea before promoting any evidence schema.

A bare formula is refused. A declared local positive curvature object is admitted.

## `rc_decay_local_curvature`

- source hook: `V(t; tau) = V0 exp(-t/tau)`
- bare formula refusal: The decay formula alone is a time-evolution surface, not a module object.
- allowed route: `local_positive_curvature`
- refused routes: hidden_load without a ceiling, static storage Hessian
- positive object: `scalar Fisher precision` on `tau` = 2392.75348727
- boundary: This proof does not infer tau from data and does not model resistor/capacitor storage; it only shows admissibility of declared local information curvature.

### Declarations

```json
{
  "coordinate": "tau",
  "coordinate_meaning": "RC time constant",
  "locality": "finite declared sample times",
  "curvature_source": "Gaussian voltage likelihood",
  "sensor_surface": "voltage samples V(t_i)",
  "noise_model": "independent Gaussian voltage noise with standard deviation sigma",
  "known_controls": {
    "V0": 5.0,
    "sample_times": [
      0.25,
      0.7666666666666667,
      1.2833333333333334,
      1.8000000000000003,
      2.316666666666667,
      2.8333333333333335,
      3.3500000000000005,
      3.866666666666667,
      4.383333333333334,
      4.9,
      5.416666666666667,
      5.933333333333334,
      6.450000000000001,
      6.966666666666668,
      7.483333333333334,
      8.0
    ],
    "sigma": 0.05
  }
}
```

### Derivation

- For independent Gaussian observations with fixed sigma, the local Fisher curvature is sum_i (1/sigma^2) (dV(t_i;tau)/dtau)^2.
- For V(t;tau) = V0 exp(-t/tau), dV/dtau = V(t;tau) * t / tau^2.
- Substituting the declared sample grid gives a positive scalar Fisher object.

## `beer_lambert_local_curvature`

- source hook: `A(c) = epsilon b c`
- bare formula refusal: The absorbance formula alone is a sensor law, not a module object.
- allowed route: `local_positive_curvature`
- refused routes: hidden_load without a ceiling, static storage Hessian
- positive object: `scalar Fisher precision` on `c` = 14400
- boundary: This proof does not model chemistry thermodynamics or concentration dynamics; it only shows admissibility of declared local information curvature from a sensor law.

### Declarations

```json
{
  "coordinate": "c",
  "coordinate_meaning": "solute concentration",
  "locality": "single declared absorbance measurement protocol; linear law is global in this coordinate under the stated regime",
  "curvature_source": "Gaussian absorbance likelihood",
  "sensor_surface": "absorbance A",
  "noise_model": "Gaussian absorbance noise with standard deviation sigma_A",
  "known_controls": {
    "molar_absorptivity": 1.2,
    "path_length": 1.0,
    "sigma_absorbance": 0.01
  }
}
```

### Derivation

- For one Gaussian absorbance observation with fixed sigma_A, the Fisher curvature is (1/sigma_A^2) (dA/dc)^2.
- For A(c) = epsilon b c, dA/dc = epsilon b.
- Substituting the declared sensor constants gives a positive scalar Fisher object.
