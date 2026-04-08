# Citations

## How To Cite This Workspace

If you use the code, tests, demonstrations, or evidence bundles from this
workspace, cite: https://doi.org/10.5281/zenodo.19474899

- the `nomogeo` kernel 
- the `nomodescent` layer if you use descent / compatibility results
- the `evidence` layer if you use evidence encoding, micro-real bundles, or
  curated ingestion

If you discuss a micro case study, also cite the source paper used for that
case study. The case studies are not standalone scientific sources.

## Internal Project Sources

Core mathematical and validation sources live in the local `papers/` tree and
the corresponding docs under `docs/`.

Canonical DOI list for the local paper set:
- [papers/DOIs.txt](papers/DOIs.txt)

Current DOI entries:
- `10.5281/zenodo.19102609`
- `10.5281/zenodo.19154445`
- `10.5281/zenodo.19361501`
- `10.5281/zenodo.19473828`
- `10.5281/zenodo.19473822`

Key workspace docs:

- [theorem_map.md](docs/theorem_map.md)
- [inverse_theorem.md](docs/inverse_theorem.md)
- [validation_note.md](docs/validation_note.md)
- [stack_validation_note.md](docs/stack_validation_note.md)

## Demonstrations

If you refer to a release demonstration, cite the demonstration itself and the
standard external quantity or method it is compared against where relevant.

- entanglement demonstration
  - standard Gaussian mutual-information formulas
- Bell / common-gluing demonstration
  - correlator-level CHSH / arcsine compatibility methods
- Arrow observer-rank demonstration
  - synthetic worked example; no external empirical source is implied

## Micro Case Studies

These are tiny cleaned application probes. They should always be cited together
with their source paper.

### 2604.06078

- Title: *A proximal approach to the Schrodinger bridge problem with incomplete information and application to contamination tracking in water networks*
- Use in this repo: observer/coarsening translation of the paper's explicit
  partial-observation / non-uniqueness examples
- Case study path:
  [tests/micro_case_studies/2604_06078_observability](tests/micro_case_studies/2604_06078_observability)

### 2604.06108

- Title: *Investigating ACS/WFC Amp-to-Amp Sensitivities*
- Use in this repo: small protocol-mismatch probe on the paper's exact red/blue
  table slice
- Case study path:
  [tests/micro_case_studies/2604_06108_amp_protocol](tests/micro_case_studies/2604_06108_amp_protocol)

### 2604.06064

- Title: *Star-planet magnetic interactions in photoevaporating exoplanets*
- Use in this repo: coarsening vs non-nested observer split from the published
  summary table
- Case study path:
  [tests/micro_case_studies/2604_06064_spmi_degeneracy](tests/micro_case_studies/2604_06064_spmi_degeneracy)

### 2604.06099

- Title: *Extending to Robust Medical Imaging: Corruption and Adversarial Stress Testing in Low-Data Regimes*
- Use in this repo: benchmark-summary observer split between realistic-shift
  and adversarial-stress views
- Case study path:
  [tests/micro_case_studies/2604_06099_robustness_views](tests/micro_case_studies/2604_06099_robustness_views)

## Evidence Bundles And Provenance

For any micro-real bundle under `evidence/micro_real_bundles`, follow the
bundle-local files as the authoritative citation / provenance surface:

- `PROVENANCE.md`
- `LICENSE.md`
- `README.md`

That local surface takes precedence over any short summary here.

