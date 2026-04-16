# Source Law on QM9 Molecular Data

**Date:** 15 April 2026
**Code:** `0_4_0_qm9_source_law.py`
**Data:** QM9 raw xyz (35,693 molecules, first 500 used)

---

## What was done

Applied the source law and conservation law to real molecular data from the QM9 dataset. Constructed diagonal Hessians from vibrational frequencies and computed observer geometry for isomer pairs (molecules with the same formula but different structures).

## Results

### Conservation on 72 isomer pairs

The first-order conservation law vis\_rate + hid\_rate = ambient\_rate was verified on 72 isomer pairs with up to 54 vibrational modes. **Maximum conservation error: 3.55e-15.** The law holds at machine precision on real molecular data.

### Information budget for molecular isomers

| Statistic | Value |
|-----------|-------|
| Pairs analysed | 72 |
| V > 0 (source law active) | 3/72 (4%) |
| Mean visible rate | 0.14 |
| Mean hidden rate | 2.81 |
| Mean visible fraction | 9.7% |

**Interpretation:** Observing only the 3 softest vibrational modes captures ~10% of the total information change between isomers. The remaining 90% lives in the hidden (higher-frequency) modes. This is the information budget for molecular observation: the softest modes are the most physically accessible but carry a small share of the total.

### V > 0 is rare on real data

Only 3 of 72 isomer pairs have V > 0 (all visible eigenvalues positive). For most pairs, the perturbation from one isomer to another increases some soft-mode frequencies while decreasing others, giving mixed-sign V. The source law's support-stability condition is restrictive on real molecular data.

This is the **information-destruction locus in chemistry**: molecular isomerisation generically destroys information in some visible directions while creating it in others.

### Single-molecule observer sweep (methane, CH4)

For methane (9 vibrational modes, frequencies 1341-3152 cm^{-1}), the conservation law was verified for 8 different observers (m = 1 to 8):

| m | Visible fraction | Conservation error |
|---|-----------------|-------------------|
| 1 | -14% | 0.0 |
| 2 | 31% | 0.0 |
| 3 | -5% | 0.0 |
| 4 | -51% | 5.5e-17 |
| 5 | 32% | 0.0 |
| 6 | 90% | 0.0 |
| 7 | 84% | 0.0 |
| 8 | 99% | 0.0 |

At m=4, the observer captures -51% of the ambient rate (working against the geometry). At m=8, it captures 99%. The observer is a splitter — the total is always -0.41.

## What the Hessian data will add

The raw QM9 xyz data gives diagonal Hessians (from frequencies only). The HessianQM9 Arrow dataset provides full 3N x 3N Hessian matrices in four solvent environments (vacuum, THF, toluene, water). This will enable:

1. **Solvent perturbation paths:** H(solvent) as a path through SPD space. The source law tracks how visible precision changes across solvents — a genuine chemical prediction.

2. **Off-diagonal Hessian structure:** Full matrices have off-diagonal coupling that diagonal frequency-based Hessians miss. This is where the hidden defect Q\_hat and the curvature F\_alpha become nontrivial.

3. **Mass-weighted rigid-motion projection:** The full Hessian includes rigid-body modes that must be projected out. The existing `molecular_vibrational_atlas.py` handles this.

4. **Observer design for molecular spectroscopy:** Which vibrational modes should a spectroscopist observe to maximise information about solvent effects? The conservation law + optimal observer alignment (R5) answers this directly.
