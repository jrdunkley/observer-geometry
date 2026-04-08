# Theorem Map

This map uses exact theorem / proposition / definition labels from the TeX sources where they exist. Where the paper statement is an unlabeled displayed formula or a code-relevant internal formula rather than a theorem-like environment, the fallback cites the exact equation label or section formula.

## Geometry Core

- `visible_precision`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `[Visible precision]`, label `prop:visible-precision`
  - formulas: `eq:visible-variational`, `eq:tower-law`, `eq:schur-form`
  - quotient confirmation: `Geometric_Consciousness_MASTER.tex`, Definition `def:induced-visible-object`, Proposition `prop:variational`, Proposition `prop:composition`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `canonical_lift`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:lift-projector`
  - theorem using it: `thm:local-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `hidden_projector`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:lift-projector`
  - theorem using it: `thm:local-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

- `local_visible_calculus`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Local visible calculus]`, label `thm:local-calculus`
  - formulas: `eq:local-expansion`, `eq:VQ-def`, `eq:Q-gram`, `eq:det-split`
  - quotient descent confirmation: `Quotient_Descent_MASTER.tex`, Proposition `prop:hilbert-quotient-calculus`
  - tests: [`tests/test_core.py`](../tests/test_core.py), [`tests/test_release_gate_core.py`](../tests/test_release_gate_core.py)

## Hidden-Load Layer

- `hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Hidden-load geometry]`, label `thm:hidden-load`
  - formulas: `eq:hidden-load-def`, `eq:X-from-Lambda`, `eq:rank-and-clock`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load`, Corollary `cor:hidden-load-rank-volume`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Definition `def:latent-return`, Proposition `prop:hidden-load`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `visible_from_hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `thm:hidden-load`
  - formula: `eq:X-from-Lambda`
  - inverse-theorem role: fixed-ceiling bijection `Lambda <-> X` on active support `S = Ran(T)`
  - domain note: `Lambda` must be positive semidefinite on active support; full-rank ceilings require explicit coordinate interpretation via `lambda_representation='reduced'` or `'ambient'`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_inverse_theorem.py`](../tests/test_inverse_theorem.py)

- `inverse_visible_class`
  - source: same fixed-ceiling inverse theorem as `visible_from_hidden_load`, promoted from `eq:X-from-Lambda` in `thm:hidden-load`
  - role: theorem-facing alias stating the ceiling-conditioned inverse map explicitly
  - non-claim: this is not an inverse of the global fibre map `(H, C) -> Phi_C(H)`
  - tests: [`tests/test_inverse_theorem.py`](../tests/test_inverse_theorem.py)

- `canonical_hidden_realisation`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:canonical-hidden`
  - theorem context: `thm:hidden-load`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

- `minimal_hidden_realisation`
  - source fallback: `Geometric_Observation_MASTER.tex`, formula `eq:min-hidden`
  - rank / gauge confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:min-hidden-dim`, Proposition `prop:bridge-hidden-factor-gauge`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

- `transport_hidden_load`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `[Transport law and determinant clock]`, label `prop:transport`
  - formulas: `eq:transport-law`, `eq:clock-add`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport`
  - domain note: theorem-domain kernel only; both inputs must satisfy `Lambda >= 0` on active support
  - structural note: this is the two-step Gram-shadow transport law, not the associative downstairs object
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `hidden_contraction`
  - source fallback: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport` together with the preceding remark `[The visible law is not the associative object]`
  - formulas: `Pi = (I_S + Lambda)^(-1)` from `thm:hidden-load-transport`, exposed as the canonical factor `K = Pi^(1/2) = (I + Lambda)^(-1/2)`
  - role: exposes the actual associative downstairs object used for long-chain composition
  - tests: [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `load_from_hidden_contraction`
  - source fallback: `Quotient_Observation_MASTER.tex`, Theorem `thm:hidden-load-transport` with `Pi = K^T K`
  - formulas: `Lambda = Pi^(-1) - I`
  - role: converts the associative downstairs object back into hidden-load coordinates
  - tests: [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- `clock`
  - source: `Geometric_Observation_MASTER.tex`, Proposition `prop:transport`
  - formula: `eq:clock-def`
  - quotient confirmation: `Quotient_Observation_MASTER.tex`, Corollary `cor:det-clock`, Theorem `thm:clock-uniqueness`
  - domain note: theorem-domain kernel only; `Lambda` must satisfy `Lambda >= 0` on active support
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py), [`tests/test_release_gate_hidden.py`](../tests/test_release_gate_hidden.py)

- local hidden birth smoke
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Local hidden return as hidden-load birth]`, label `thm:birth`
  - formula: `eq:birth-law`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:local-birth`
  - tests: [`tests/test_hidden.py`](../tests/test_hidden.py)

## Finite Bridge

- `dv_bridge`
  - source: `Geometric_Observation_MASTER.tex`, Theorem `[Donsker-Varadhan bridge]`, label `thm:dv-bridge`
  - formulas: `eq:dv-bridge`, `eq:dv-gram`
  - finite source confirmation: `Finite_Observation_MASTER.tex`, Theorem `[Quadraticity of the DV correction]`, label `thm:exact_quadraticity`
  - quotient descent confirmation: `Quotient_Descent_MASTER.tex`, Proposition `prop:dv-hilbert-lift`
  - tests: [`tests/test_bridge.py`](../tests/test_bridge.py), [`tests/test_release_gate_bridge.py`](../tests/test_release_gate_bridge.py)

- DV quadratic onset and quartic hidden defect smoke
  - source: `Geometric_Observation_MASTER.tex`, Corollary `[Quadratic onset and quartic defect]`, label `cor:quartic-defect`
  - source: `Geometric_Observation_MASTER.tex`, Corollary `[The quartic defect is the first hidden birth]`, label `cor:final`
  - consciousness confirmation: `Geometric_Consciousness_MASTER.tex`, Corollary `cor:quadratic-dv`
  - tests: [`tests/test_bridge.py`](../tests/test_bridge.py)

## Quotient Adapter Layer

- `observed_covariance` and Gaussian divergences
  - source context: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:gaussian-dp`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_quotient.py`](../tests/test_release_gate_quotient.py)

- `observer_collapse_descends`
  - source: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:collapse-asymmetry`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py), [`tests/test_release_gate_quotient.py`](../tests/test_release_gate_quotient.py)

- coarse observer misses a real change smoke
  - source: `Geometric_Consciousness_MASTER.tex`, Proposition `prop:toy-example`
  - tests: [`tests/test_quotient.py`](../tests/test_quotient.py)

## Batch Layer

- `batch_map`
  - source fallback: implementation requirement from `PLAN.txt`, section `Multicore support`
  - no theorem label exists for this engineering layer
  - tests: [`tests/test_batch.py`](../tests/test_batch.py)

