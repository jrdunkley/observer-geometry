# Molecular Observer Geometry Atlas v0 Status

Status: first acquisition and atlas pass completed inside `src/utils/operationalise/attachment_layer`.

## Data Acquired

- Dataset: Hessian QM9, Figshare article `26363959`, version 4.
- Archive: `incoming_data/hessian_qm9/hessian_qm9_DatasetDict.zip`.
- Official archive MD5: `f3e36130e5cc47021ab403767a19ddf7`.
- Local archive MD5: `f3e36130e5cc47021ab403767a19ddf7`.
- License: CC0 in Figshare metadata.
- Extracted for this pass:
  - all five `vacuum` Arrow shards;
  - first Arrow shard for `thf`, `toluene`, and `water`.

## Declared Compiler Route

The v0 route is:

`molecular record -> Cartesian Hessian -> mass-weighted Hessian -> rigid translation/rotation removal -> vibrational SPD object -> transformed atom-selection observer -> visible_precision`

Important declaration:

Cartesian atom-selection observers are not passed directly to the kernel. They are transformed into the mass-weighted vibrational coordinate chart before `visible_precision` is called.

Boundary:

No `hidden_load` is emitted in this v0 atlas, because atom-selection observers after vibrational projection do not yet carry a declared ceiling convention.

## Full Vacuum Pass

- Vacuum records scanned: `41,645`.
- Molecules admitted as projected vibrational SPD objects: `40,583`.
- Molecule-level refusals: `1,062`.
- Observer-level refusals: `8`.
- Observer readings emitted: `243,479`.
- Admission rate: `0.9745`.

Vibrational Hessian conditioning over admitted molecules:

- condition min: `5.2828`.
- condition median: `602.2289`.
- condition q90: `1786.4310`.
- condition q99: `11159.9846`.
- condition max: `4550552.9331`.

Minimum vibrational eigenvalue over admitted molecules:

- min: `1.1221e-05`.
- q01: `0.0044170`.
- q10: `0.0269048`.
- median: `0.0796517`.
- max: `11.1153`.

The molecule-level refusals are currently dominated by negative projected vibrational eigenvalues. That is a useful gate result, not a nuisance: the atlas is refusing records that fail the declared local positive-curvature route.

## Spectrum Validation

A 250-row vacuum sample was used to check the compiler against the dataset-provided frequency record.

Method:

- compile the projected vibrational Hessian with the v0 compiler;
- take the square roots of its positive eigenvalues;
- compare the sorted log spectrum against the largest `vibrational_dim` dataset frequency magnitudes, leaving the six smallest rigid/near-rigid entries outside the comparison for nonlinear molecules.

Result:

- admitted and compared: `248`.
- refused during compile: `2`.
- log-spectrum correlation min: `0.906685`.
- q10: `0.997187`.
- median: `0.999911`.

Interpretation:

The mass-weighted rigid-projection compiler is reading the Hessian field in the right convention up to unit scaling. The low outliers should be inspected, but the median and q10 are strong enough to treat the v0 compiler as a valid first-pass extraction route.

## Observer Signatures

Median visible geometric mean eigenvalue by observer type:

| Observer | Readings | Median rank | q10 | median | q90 | median condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heavy atoms | 40,583 | 27 | 37.1501 | 48.6561 | 64.1354 | 1636.6193 |
| hydrogens | 40,575 | 27 | 4.7934 | 5.3374 | 6.0682 | 86.8917 |
| local heavy site plus attached H | 162,321 | 6 | 7.2993 | 11.1839 | 25.7937 | 55.6400 |

Heavy-atom versus hydrogen observer separation:

- heavy/hydrogen visible geometric-mean ratio q10: `7.3350`.
- median: `8.9929`.
- q90: `11.5019`.
- max: `150.7778`.

Interpretation:

The heavy skeleton, hydrogen surface, and local heavy-site observers are already separated by the geometric reader at scale. This does not yet claim chemical classification, but it is the first strong signal that simple atom-selection observers produce nontrivial molecular access geometry.

Composition sanity check:

For molecules with both heavy-atom and hydrogen observers, the heavy/hydrogen visible geometric-mean ratio shows:

- correlation with heavy-atom count: `-0.3005`;
- correlation with hydrogen count: `-0.8487`;
- correlation with heteroatom count: `0.5593`.

Median heavy/hydrogen ratio by heteroatom count:

| Heteroatom count | Molecules | Median ratio |
| ---: | ---: | ---: |
| 0 | 1,014 | 7.399 |
| 1 | 6,113 | 7.931 |
| 2 | 13,702 | 8.600 |
| 3 | 13,048 | 9.469 |
| 4 | 5,295 | 10.730 |
| 5 | 1,220 | 12.443 |
| 6 | 168 | 14.619 |

Interpretation:

This is a real atlas signal but not yet a chemistry result. It may reflect heteroatom stiffness, hydrogen count, formula class, or all of these together. The next pass should separate these effects instead of claiming functional-group recovery too early.

## Solvent Perturbation Probe

Probe declaration:

For first-shard same-label records, each solvent Hessian is mass-weighted, mass-weighted-Kabsch aligned to the vacuum geometry, and pulled back to the vacuum vibrational chart before forming solvent perturbations. The perturbation family is then read with `closure_scores`.

First-shard common labels across vacuum, THF, toluene, and water: `54`.

Rows emitted: `216`.

Alignment RMSD:

- median: `0.0093269`.
- max: `0.177603`.

Median solvent leakage fraction `eta` by observer:

| Observer | Readings | eta q10 | eta median | eta q90 |
| --- | ---: | ---: | ---: | ---: |
| heavy atoms | 54 | 0.0301 | 0.1756 | 0.3569 |
| hydrogens | 54 | 0.0453 | 0.0907 | 0.1986 |
| local heavy site plus attached H | 108 | 0.0731 | 0.2762 | 0.6075 |

Interpretation:

The solvent route is promising but more convention-dependent than the vacuum atlas. Under the declared pullback convention, local heavy-site observers leak more of the solvent perturbation family than whole heavy-atom or hydrogen observers. This should be stress-tested before it becomes chemistry-facing evidence.

## Current Scientific Read

This is no longer a toy attachment case. The module has now read tens of thousands of molecular curvature objects under declared atom observers.

The strongest current result is the contact grammar itself:

`finite molecular curvature -> declared observer -> visible molecular access geometry -> refusal when positive curvature fails`

What is integrating well:

- equilibrium vibrational Hessians as local positive curvature objects;
- atom-selection observers after mass-weighted vibrational projection;
- full-vacuum batch execution with stable admission/refusal behaviour;
- solvent perturbation as a possible closure/leakage family after explicit alignment and pullback.

What remains weak:

- no ceiling convention yet for molecular hidden load;
- no detailed validation yet against dataset-provided normal-mode vectors;
- no chemistry labels or isomer grouping yet, so semantic recovery is not tested;
- solvent comparisons need stronger coordinate-chart and alignment stress tests;
- local-site observers use a simple covalent-radius heuristic, not a chemistry-certified graph.

## Next Best Moves

1. Validate the vibrational compiler against Hessian QM9 frequencies and normal modes.
2. Separate numerical near-zero refusals from chemically meaningful indefinite records.
3. Add molecule fingerprints: formula, heavy-atom count, heteroatom count, hydrogen count, eigenvalue summaries, observer ratios.
4. Test whether observer signatures cluster by formula or inferred composition before adding any external chemistry labels.
5. Add a declared ceiling convention only after deciding what visible-block reference means in the vibrational observer chart.
6. Stress-test the solvent pullback convention across all shards only after the normal-mode validation passes.
