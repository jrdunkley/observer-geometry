# Claim Hierarchy

This workspace now contains several epistemic layers. They should not be read as if they make the same kind of claim.

## 1. Exact Kernel Claims

These live in `nomogeo`.

Meaning:

- theorem-facing matrix identities inside the declared linear / Gaussian domain
- exact fixed-ceiling hidden-load inverse theorem
- exact quotient / tower identities
- exact deterministic batch semantics

Status:

- theorem-grade within the documented domain and tolerances

Examples:

- visible-precision identities
- hidden clock identity
- inverse round-trips beneath fixed ceiling `T`
- exact two-step hidden-load transport law

## 2. Exact Descent Claims

These live in `nomodescent`.

Meaning:

- exact linear observer relation tests
- exact factorisation and non-factorisation
- exact tower checks
- exact incompatibility when linear constraints already fail
- exact incompatibility when a unique affine completion is non-PSD
- exact underdetermination when the current exact engine deliberately refuses to guess

Status:

- theorem-grade within the declared linear / Gaussian domain

Examples:

- `exact_factorisation`
- `non_nested_observers`
- `incompatible_by_linear_inconsistency`
- `underdetermined_affine_family`

## 3. Audited Approximate Claims

These currently appear mainly in `nomodescent` and in some micro-real bundles.

Meaning:

- the underlying exact object is known
- the final positive / negative conclusion depends on a finite deterministic approximation layer
- that approximation layer must be exposed in the audit

Status:

- not theorem-grade
- acceptable only when explicitly labelled and accompanied by audit and falsification routes

Examples:

- `approximate_common_descent`
- `incompatible_by_approximate_psd_search`
- Gaussian-surrogate layers built from real non-Gaussian evidence

## 4. Synthetic Worked Examples

These are examples chosen for clarity of geometry, not because they are direct empirical datasets.

Meaning:

- explicit, executable examples
- often mathematically cleaner than real-data cases
- useful for pedagogy, architecture validation, and scope control

Status:

- can contain exact theorem demonstrations
- can also contain explicitly synthetic interpretation layers
- must never be presented as empirical evidence

Examples:

- entanglement example: exact identity demonstration in a clean Gaussian family
- Arrow example: explicitly synthetic observer-geometry example
- `evidence` benchmark worked example: explicitly synthetic encoding example

## 5. Micro-Real Bundles

These live in `evidence/micro_real_bundles`.

Meaning:

- tiny legally clean real or manually reconstructed evidence objects
- exact extracted facts remain separate from encoded inference
- downstream conclusions inherit the epistemic status of the assembled problem

Status:

- stronger than purely synthetic examples as evidence of practical contact with real material
- still not a licence to overclaim beyond the tiny bundle and its encoding choices

Examples:

- Bell count table transcription plus explicit no-signalling projection
- tiny Iris panel mismatch bundle
- tiny public leaderboard score slice with explicit benchmark-observer ambiguity

## Reading Rule

When reading any result in this workspace, ask:

1. What was exact input?
2. What was encoded inference?
3. What remains ambiguous?
4. Is the conclusion exact, audited approximate, synthetic, or micro-real?

If those answers are not visible, the surface needs tightening.

