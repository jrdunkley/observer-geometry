# 0.3.2 Declared Frontier Local Certificate Note

Date: 2026-04-11

Status: pure research note. This is not an implementation plan and does not
modify module source.

## Executive Result

The conservative near-branch certificate can be widened conceptually. Once
the general graph-chart first variation and Hessian are available, the same
abstract stability lemma gives a declared-observer local certificate for any
observer with:

```text
small stationarity residual,
definite graph Hessian margin,
Hessian-Lipschitz control on the chart.
```

This does not require the observer to be an exact branch. Exact branches are
just the simplest stationarity sector.

## Certificate Shape

Let `F_mu` be the weighted frontier score in a graph chart around a declared
observer. Suppose on `||X||_F <= rho`:

```text
||grad F_mu(0)|| <= eps,
D^2F_mu(0) <= -lambda I,        lambda > 0,
||D^2F_mu(X)-D^2F_mu(0)|| <= L ||X||.
```

Let

```text
r0 = min(rho, lambda/(2L)).
```

If

```text
4 eps / lambda < r0,
```

then the declared observer has a unique nearby local maximizer in that chart,
with displacement bounded by

```text
||X_*|| <= 2 eps / lambda.
```

The local-minimum version is obtained by applying the same statement to
`-F_mu`.

## Lipschitz Bound

The conservative rho-dependent bound from
`audit/0_3_2_projector_bounds_note.md` can be used:

```text
L_cert(rho) =
  G3 C1(rho)^3
  + 3 G2 C2(rho) C1(rho)
  + G1 C3(rho),
```

with the norm conventions and conservatism caveats from that note.

The important point is that this bound is a graph-chart functional bound. It
does not require exact-branch invariance. It only requires the supplied
weighted family, observer rank, penalty `mu`, chart radius, and norm contract.

## Stress Checks

The script

```text
audit/0_3_2_declared_frontier_certificate_check.py
```

writes

```text
audit/outputs/0_3_2_declared_frontier_certificate_check.json
```

It checks three sectors:

- stationary non-exact symmetric-pair cases;
- near-stationary imbalanced symmetric-pair cases;
- random declared near-branch cases.

For the symmetric non-exact family

```text
A_+ = [[3,e],[ e,1]],
A_- = [[3,-e],[-e,1]],
mu = 0,
```

the max certificate passes for `e=1` and `e=1.7`, fails at the degenerate
threshold `e=sqrt(3)`, and then the min certificate passes for `e=1.75` and
`e=2`. This is the intended behavior: the certificate follows the true graph
Hessian sign, not the invalid exact-branch proxy.

For small imbalances away from perfect cancellation, the certificate becomes
controlled by the ratio `eps/lambda`. It passes at very small imbalance and
eventually becomes vacuous as the stationarity residual grows relative to the
conservative Lipschitz radius.

## Implication

The right future API, if one is ever approved, should not be named as an
approximate exact-branch Hessian. The theory points instead toward a declared
local frontier certificate:

```text
declared observer
-> stationarity residual
-> graph Hessian spectrum
-> conservative local certificate if margins beat residuals.
```

That would sit cleanly in exact local quadratic geometry. It would still not
be a global Grassmannian optimizer and not a full non-Gaussian law selector.

## Final Position

This is a real 0.3.2 research win. It turns the non-exact stationary examples
from pathologies into evidence for a broader exact local object. The limiting
factor is no longer exact-branch status; it is certificate non-vacuity under a
conservative Lipschitz bound.
