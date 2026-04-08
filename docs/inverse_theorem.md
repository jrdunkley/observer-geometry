# Fixed-Ceiling Inverse Theorem

The hidden-load layer has an exact inverse theorem, but only after fixing a
ceiling/reference visible object `T`.

## Exact Statement

Let `T >= 0` with active support `S = Ran(T)`. Write `T_S` for the restriction
of `T` to `S`, so `T_S` is strictly positive definite. Then:

```text
Lambda -> X = T_S^(1/2) (I + Lambda)^(-1) T_S^(1/2)
X -> Lambda = T_S^(1/2) X_S^(-1) T_S^(1/2) - I
```

are mutually inverse bijections between:

- positive semidefinite hidden loads on `S`
- support-preserving visible objects `X` beneath `T`

This is the theorem-level inverse surface implemented by:

- `hidden_load(T, X)`
- `visible_from_hidden_load(T, Lambda, ...)`
- `inverse_visible_class(T, Lambda, ...)`

## What It Does Not Claim

This is not a global inverse of the forward map:

```text
(H, C) -> Phi_C(H)
```

That global latent-fibre problem remains highly non-unique. The inverse theorem
only recovers the exact hidden-load compatibility class relative to a fixed
ceiling/reference `T`.

## Structural Consequences

- order reversal:
  `Lambda_2 >= Lambda_1` implies `X_2 <= X_1`
- rank theorem:
  `rank(Lambda) = rank(T - X)`
- clock identity:
  `log pdet(T) - log det(X|S) = log det(I + Lambda)`

## Coordinate Rules

- the hidden operator lives on active support `S = Ran(T)`
- `hidden_load` returns both reduced-support and ambient coordinates
- `visible_from_hidden_load` and `inverse_visible_class` accept either:
  - reduced-support coordinates with `lambda_representation="reduced"`
  - ambient coordinates with `lambda_representation="ambient"`
- when `rank(T) = n`, reduced and ambient shapes coincide, so the
  representation must be stated explicitly

## Numerical Boundary Notes

- zero-support ceilings are handled exactly and return the zero visible object
- near the ceiling, tiny hidden loads make rank-sensitive statements threshold
  dependent at the machine floor; that is a tolerance issue, not a theorem
  failure
- large hidden loads remain within theorem domain, but the visible object on
  support becomes correspondingly ill-conditioned

## Validation Surface

Direct inverse-theorem tests live in:

- [tests/test_inverse_theorem.py](../tests/test_inverse_theorem.py)
- [tests/test_hidden.py](../tests/test_hidden.py)
- [tests/test_release_gate_hidden.py](../tests/test_release_gate_hidden.py)

The seeded sweep also records direct inverse residuals in:

- [tools/validation_sweep.py](../tools/validation_sweep.py)

