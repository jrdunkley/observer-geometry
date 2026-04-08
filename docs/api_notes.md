# API Notes

## Public Entry Points

- `visible_precision(H, C)` returns the visible precision matrix.
- `canonical_lift(H, C)` returns `L_{C,H}`.
- `hidden_projector(H, C)` returns `P_{C,H}`.
- `visible_geometry(H, C)` returns all three together as a `VisibleGeometryResult`.
- `local_visible_calculus(H, C, Delta)` returns `phi`, `V`, `Q`, the lift/projector, and the determinant-curvature split.
- `dv_bridge(H0, Jhat)` returns `H_DV`, `Delta_DV`, and the Gram factor.
- `hidden_load(T, X)` returns the support-aware hidden load and metadata.
- `visible_from_hidden_load(T, Lambda, ..., lambda_representation=...)` reconstructs `X`.
- `inverse_visible_class(T, Lambda, ..., lambda_representation=...)` is the theorem-facing alias for the same fixed-ceiling inverse map.
- `hidden_contraction(Lambda)` returns the canonical contraction factor `K = (I + Lambda)^(-1/2)`.
- `load_from_hidden_contraction(K)` recovers the hidden load from `K^T K`.
- `canonical_hidden_realisation(T, X)` and `minimal_hidden_realisation(T, X)` work on active-support coordinates and return explicit block matrices.
- `transport_hidden_load(Lambda, M)` returns the theorem-native two-step transported load.
- `clock(Lambda)` returns `log det(I + Lambda)`.
- `observed_covariance`, `gaussian_forward_kl`, `gaussian_reverse_kl`, `gaussian_hellinger_squared`, and `gaussian_data_processing_contraction` provide the thin quotient adapter layer.
- `batch_map(func, tasks, ...)` is the thin deterministic batch wrapper.

## Support-Aware Convention

The hidden-load layer follows the papers: the real operator lives on `S = Ran(T)`.

- `hidden_load` returns:
  - `lambda_`: ambient embedding with zeros off support
  - `reduced_lambda`: the active-support operator
  - `support_basis`: the basis used for the active support
- `visible_from_hidden_load` accepts either:
- `inverse_visible_class` accepts the same coordinate choices and exists to make the fixed-ceiling inverse theorem explicit at the API surface
  - an ambient `n x n` load with `lambda_representation="ambient"`, or
  - a reduced `rank(T) x rank(T)` load with `lambda_representation="reduced"`
- when `rank(T) = n`, shape alone is ambiguous and the function requires the representation to be stated explicitly
- this inverse surface is exact only after fixing a ceiling/reference `T`; it does not invert the global forward map `(H, C) -> Phi_C(H)`
- `clock` and `transport_hidden_load` are theorem-domain functions and reject indefinite hidden loads
- the associative downstairs object is the contraction factor `K = (I + Lambda)^(-1/2)`, not the raw visible-effect or raw load coordinate
- long chains should therefore be composed in contraction-factor coordinates and converted back with `load_from_hidden_contraction`

`canonical_hidden_realisation` and `minimal_hidden_realisation` return support-coordinate block matrices, because that is the canonical domain of the theorem.

## Result Objects

The package uses small frozen dataclasses instead of mutable object graphs:

- [`src/nomogeo/types.py`](../src/nomogeo/types.py)

Each result carries `LinearAlgebraMetadata` with tolerance, support rank, and method information.

