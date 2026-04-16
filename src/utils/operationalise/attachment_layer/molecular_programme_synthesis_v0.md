# Molecular Programme Synthesis v0

Status: synthesis after the full Hessian QM9 visible-atlas pass, four-environment support pass, solvent perturbation probe, and within-formula isomer geometry analysis.

This note is not a new theorem layer. It records what the molecular route is currently teaching the operationalisation programme, with hard claims separated from leads.

## Core Contact

The molecular route now has a working compiler:

`Cartesian Hessian -> mass-weighted rigid projection -> declared vibrational SPD chart -> atom-selection observer -> observer-geometry reading`

This is the cleanest chemistry contact point found so far. It is not chemistry labels first; it is curvature first. The entry object is the molecular vibrational Hessian at a local minimum, and the first observers are simple atom-selection windows.

The doctrine remains strict:

- a projected vibrational Hessian admitted as SPD can enter the visible-precision reader;
- a projected vibrational Hessian that is indefinite or singular is a refusal for that route;
- rigid translations and rotations are support provenance, not discarded junk;
- raw zero-mode or soft-mode facts are not hidden-load facts without a declared ceiling, path, or event coefficient.

## Full-Batch Facts

The visible vacuum atlas admitted `40,583` of `41,645` records, for an admission rate of `0.9745`. It produced `243,479` observer readings over heavy-atom, hydrogen, and local-heavy-site windows.

The support atlas processed all four environments, `166,580` rows total, with zero parser failures:

| Environment | Rows | SPD | Indefinite/singular | SPD rate | softest positive eig median | condition median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| vacuum | 41,645 | 40,583 | 1,062 | 0.9745 | 0.0784205 | 610.833 |
| THF | 41,645 | 40,391 | 1,254 | 0.9699 | 0.0284726 | 1603.928 |
| toluene | 41,645 | 41,264 | 381 | 0.9909 | 0.0341577 | 1345.811 |
| water | 41,645 | 39,811 | 1,834 | 0.9560 | 0.0213982 | 2172.276 |

The dominant cross-environment support pattern is `SSSS` with `37,799` molecules. The strongest neighboring refusal pattern is `SSSI` with `1,446` molecules: admitted in vacuum, THF, and toluene, but refused in water under the declared projection/tolerance convention.

## Signals

Heavy-atom windows see a much stiffer visible geometry than hydrogen windows in the vacuum atlas. The median heavy/hydrogen visible-geometric-mean ratio is about `8.99`, and it rises with heteroatom count from about `7.40` at zero heteroatoms to about `14.62` at six heteroatoms.

Local-site windows are sensitive to local saturation. The median local-site visible geometric mean decreases with attached hydrogens:

| Attached H count | Median local visible gmean |
| ---: | ---: |
| 0 | 20.729 |
| 1 | 12.318 |
| 2 | 9.783 |
| 3 | 7.741 |

Environment shifts are strong in the support layer. Relative to vacuum, the median log softest-eigenvalue ratios are:

| Environment | Median log ratio |
| --- | ---: |
| THF | -0.866 |
| toluene | -0.724 |
| water | -1.152 |

Water is both the softest and the least often admitted environment in this pass. Toluene is the cleanest by projected-SPD admission rate.

Soft-mode participation also shifts. Median hydrogen participation in the softest positive mode rises from `0.179` in vacuum to `0.215` in THF, `0.207` in toluene, and `0.222` in water. The most frequent top soft-mode element shifts from carbon in vacuum to oxygen in the solvent splits.

## Within-Formula Structure

The isomer analysis asks a narrow empirical question: if formula is held fixed, does observer/support geometry still vary enough to generate useful leads?

Answer: yes.

Across `40,583` vacuum-admitted molecules, there are `141` formula groups with at least 20 records. Several same-formula groups have multiple support patterns and substantial observer-fingerprint spread. High-entropy groups include:

| Formula | n | support patterns | entropy |
| --- | ---: | ---: | ---: |
| C4H9N3O2 | 27 | 5 | 0.903 |
| C4H8N2O3 | 35 | 5 | 0.866 |
| C5H11NO3 | 51 | 6 | 0.783 |
| C7H15NO | 64 | 6 | 0.774 |
| C4H6N4O | 177 | 8 | 0.764 |

This does not prove chemical novelty by itself. It does show that the reader is not merely reproducing elemental composition. There is within-formula geometry left to explain.

## Best Leads

1. `SSSI` water-only support refusals. These are the cleanest environment-specific refusal leads because the same label is admitted in the other three environments.

2. High-entropy formula families. These are the best places to inspect isomer or conformer structure because composition is fixed while support and observer geometry move.

3. Coarse-collapse candidate pairs. These are same-formula pairs that are almost identical under a coarse heavy-atom fingerprint but separated by richer hydrogen/local-site/water-softening coordinates. The first graph inspection found that `0` of the top `25` pairs share the same simple bond/local-environment/cyclomatic signature, so the strongest current pair leads mostly look like ordinary graph-isomer separation rather than deeper same-graph curvature separation.

4. Same-graph separations. A sharper graph-stratified pass over the `40,583` vacuum-admitted fingerprints found `20,958` heuristic graph groups, including `7,300` multi-record groups. This creates a better lead class: same graph signature, different observer/support fingerprint. The strongest current examples have large rich-fingerprint separation while the simple graph signature is held fixed, but they still require manual structure/conformer inspection.

5. Solvent softening. The water and THF shifts are large enough to deserve attention, but they must be checked against geometry alignment, dataset convention, and frequency metadata before being read chemically.

6. Heteroatom-dependent heavy/hydrogen separation. This is a simple and robust first sign that atom-selection observer geometry is tracking chemically meaningful curvature structure without labels.

## Boundaries

- The atom-selection atlas currently licenses visible precision only.
- No molecular hidden-load claim is made yet; no ceiling convention has been declared.
- The support analysis is static across separate optimized environments. It is not a continuous reaction path or support-birth/death clock.
- Formula is not a linear module observer. Same-formula comparisons are useful lead generation, not theorem-level false-collapse.
- The solvent perturbation probe uses a declared Kabsch pullback into the vacuum vibrational chart. That makes it interpretable as a convention, not as a coordinate-free solvent theorem.

## Scientific Read

The molecular route is now genuinely productive. It gives the module a high-volume physical substrate where the positive object is native, the refusal gate is strict, and simple observers already produce nontrivial structure.

The most important point is not that molecules can be forced into the module. It is that the same discipline developed for springs, LC systems, modal truncations, and local curvature now scales to `41,645` real molecular curvature objects without changing the kernel idea.

The gold is likely at the boundary between admitted geometry and structured refusal. Water-only refusals, within-formula observer splits, and soft-mode environment shifts are telling us where the current declared chart is informative and where it may be hiding chemistry, numerical convention, or both.

## Next Move

The graph-aware inspection layer is now started. The next high-leverage refinement is to move from heuristic graph stratification to chemistry-facing graph audit:

`molecular graph -> atom environment classes -> observer grouping -> lead-pair inspection`

That should be done before adding hidden load. A molecular ceiling convention will be much more credible once the observer families are chemically legible and the strongest support-refusal families have been inspected against graph structure.
