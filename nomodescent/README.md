# nomodescent

`nomodescent` is the first exact scientific superstructure built on top of the
frozen `nomogeo` kernel.

It is a narrow workbench for disciplined observer-geometry questions in the
linear / Gaussian regime:

- factorisation and nestedness of observer families
- staged quotient descent and tower checks
- common Gaussian completion and compatibility classification
- explicit obstruction certificates
- deterministic finite refinement search
- audit-first outputs with exact vs approximate boundaries made explicit

## Claim Hierarchy

`nomodescent` produces two different kinds of scientific output and they should
not be read interchangeably.

Exact descent claims:

- factorisation / non-factorisation
- nestedness / non-nestedness
- exact tower agreement
- exact linear incompatibility
- exact unique-PSD obstruction
- exact underdetermination when the affine family is not resolved

Audited approximate claims:

- deterministic low-dimensional affine PSD search when explicitly enabled
- classifications such as `approximate_common_descent`
- any result whose audit layer says the conclusion depends on approximate PSD
  search rather than an exact theorem closure

If a result is not exact, the audit object is part of the claim.

## Exact Domain

Supported in `v0.2`:

- linear observers
- Gaussian / quadratic visible evidence
- exact matrix-based factorisation and tower logic
- exact covariance completion when linear constraints already determine a
  compatible or incompatible outcome
- deterministic low-dimensional PSD-feasibility search only when explicitly
  requested and clearly audited as approximate

Deferred in `v0.2`:

- nonlinear inference
- generic heuristic search
- automated PDF extraction
- arbitrary symbolic science automation
- broad application wrappers

## Install

`nomogeo v0.1.0` must already be installed or available in the environment.

```bash
pip install -e nomodescent
```

## Quick Start

```python
from nomodescent import GoalSpec, ObserverSpec, ProblemSpec, VisibleEvidenceSpec, common_descent_test

problem = ProblemSpec(
    name="toy",
    latent_dim=2,
    observers=[
        ObserverSpec(name="full", matrix=[[1.0, 0.0], [0.0, 1.0]]),
        ObserverSpec(name="coarse", matrix=[[1.0, 0.0]]),
    ],
    evidence=[
        VisibleEvidenceSpec(name="coarse_cov", observer="coarse", kind="covariance", matrix=[[1.0]]),
    ],
    goals=[GoalSpec(kind="common_completion")],
)

result = common_descent_test(problem)
print(result.classification)
print(result.audit.theorem_layer)
```

## Worked Examples

- [worked_examples/bell_descent](worked_examples/bell_descent)
- [worked_examples/free_gaussian_rg](worked_examples/free_gaussian_rg)
- [worked_examples/replication_fragility](worked_examples/replication_fragility)

## QD Core

The narrow operational Quotient Descent core is documented in:

- [docs/QD_CORE_SCOPE.md](docs/QD_CORE_SCOPE.md)
- [external_comparisons/README.md](external_comparisons/README.md)

## Docs

- [docs/DOMAIN.md](docs/DOMAIN.md)
- [docs/HOW_TO_USE_THIS_ON_VAGUE_SCIENTIFIC_MATERIAL_WITHOUT_LYING_TO_YOURSELF.md](docs/HOW_TO_USE_THIS_ON_VAGUE_SCIENTIFIC_MATERIAL_WITHOUT_LYING_TO_YOURSELF.md)
- [docs/HOW_TO_COMPARE_QUOTIENT_DESCENT_AGAINST_EXTERNAL_METHODS_HONESTLY.md](docs/HOW_TO_COMPARE_QUOTIENT_DESCENT_AGAINST_EXTERNAL_METHODS_HONESTLY.md)
- [docs/QD_CERTIFICATES.md](docs/QD_CERTIFICATES.md)
- [docs/QD_CORE_SCOPE.md](docs/QD_CORE_SCOPE.md)
- [docs/RELEASE_NOTE.md](docs/RELEASE_NOTE.md)

