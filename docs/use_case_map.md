# Use Case Map

This map routes common project uses to the exact surfaces that currently exist.
It is intentionally conservative: every route starts from declared finite data,
not from an unconstrained prose claim, raw paper, or raw dataset.

## Reading Rule

Use `nomogeo` when the task can be reduced to one of these declared objects:

- a covariance, precision, Hessian, or Fisher matrix
- a linear observer or observer family
- a finite weighted quadratic family
- an evidence bundle with exact / inferred / ambiguous status
- an explicitly supplied affine-hidden law sector
- a support-stratum or restart datum

If the object is not yet in one of those forms, the first step is encoding, not
kernel inference.

For agent execution standards and workflow gates, see
[docs/agent_requirements.md](agent_requirements.md).

## Task Routes

| Goal | Start With | Use | Example Surface | Boundary |
| --- | --- | --- | --- | --- |
| Understand a theory, paper, or scientific idea more clearly | Claimed observers, hidden variables, Hessian/Fisher objects, or evidence bundles | `visible_precision`, `canonical_lift`, `hidden_load`, local calculus, evidence assemblies | `examples/minimal_visible_precision.py`, `examples/local_calculus.py`, `evidence/worked_examples` | The module does not read or judge raw papers by itself. |
| Check whether a claim is really supported | Exact inputs plus declared inferred or ambiguous steps | `docs/claim_hierarchy.md`, `evidence`, residual-margin diagnostics, declared frontier certificates | `tests/micro_case_studies/README.md`, `examples/*/CLAIMS.md` | A supported kernel result is not automatically a supported empirical claim. |
| Find the right approach to a hard problem | Several candidate observers, quotients, branches, or residual bounds | closure scores, hidden-load coordinates, `same_rank_observer_comparison`, `weighted_family_frontier_scores`, `declared_frontier_local_certificate` | `examples/graph_frontier_declared_certificate` | These tools compare supplied candidates; they are not global observer discovery. |
| Push research forward | A precise theorem boundary or suspected breakpoint | non-Gaussian stress tests, affine-hidden reductions, branch Hessians, graph-frontier Hessians | `papers/0.3.2_Technical_Note_1.tex` | Local quadratic survival is not full non-Gaussian law exactness. |
| Build apps and live tools | Stable finite kernel calls and deterministic inputs | `batch_map`, visible precision, hidden load, local calculus, interval diagnostics, evidence encoders | public demos in `examples/` | Apps should expose provenance and theorem scope, not hide it behind a generic answer. |
| Apply it to datasets | A covariance, Hessian/Fisher estimate, weighted family, or encoded evidence bundle | `visible_precision`, local quadratic ensembles, `rank_k_covariance_perturbation`, evidence problem assembly | `evidence/micro_real_bundles`, `tests/micro_case_studies` | Raw data ingestion and statistical estimation are upstream of the kernel. |
| Compare different observers | A finite observer list, same-rank family, or declared ladder | same-rank comparison, leakage / visibility scores, `declared_ladder_dimension_cost_intervals`, branch and frontier Hessian diagnostics | `examples/declared_ladder_cost.py`, `examples/graph_frontier_declared_certificate` | Ladder intervals rank only the supplied ladder, not all possible observers. |
| Diagnose and understand failures | Noncommuting families, residual gaps, conditioning issues, or support events | closure certificates, leakage channels, residual margins, support-stratum transport, `guarded_fibre_dominance` | `examples/affine_hidden_branch_reversal.py`, support-related tests | Matrix support strata are PSD rank/kernel events, not probability atoms. |
| Build simpler models that keep important structure | A ceiling, quotient, hidden space, or composition chain | fixed-ceiling inverse theorem, hidden load, quotient precision, hidden contraction factors | `examples/hidden_load_transport.py`, `docs/inverse_theorem.md` | The ceiling is part of the theorem data; it is not recovered globally from visibility alone. |
| Track changes, branches, and critical events | A path, finite branch family, support basis, or affine-hidden fibre | local visible calculus, kernel Schur jets, support transport/restarts, affine-hidden staged elimination, branch-reversal diagnostics | `examples/affine_hidden_branch_reversal.py`, `examples/graph_frontier_declared_certificate` | Branch probabilities and global event laws need additional law data beyond Hessians. |

## Strong Current Uses

- finite observer comparison under supplied quadratic geometry
- exact Gaussian-law demonstrations
- exact local Hessian/Fisher geometry outside Gaussian law mode
- exact affine-hidden Gaussian-fibre reductions
- declared branch/frontier local certificates
- explicit evidence encoding and claim-status separation

## Do Not Use It As

- a generic non-Gaussian marginalisation engine
- a global Grassmannian optimizer
- a PDF reader or scientific-claim oracle
- a raw dataset analysis pipeline before estimation/encoding
- a branch-probability calculator from Hessians alone
- a support-event detector for probability atoms without law-level support data
