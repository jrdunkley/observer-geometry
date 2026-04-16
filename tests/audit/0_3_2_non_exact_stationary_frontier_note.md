# 0.3.2 Non-Exact Stationary Frontier Note

Date: 2026-04-11

Status: pure research note. This is not an implementation plan and does not
modify module source.

## Executive Result

Exact branches are sufficient for weighted-frontier stationarity, but they are
not necessary. A weighted family can have a stationary declared observer even
when individual family members do not preserve that observer. The mechanism is
first-variation cancellation.

This is a real theory boundary:

- `exact_branch_hessian` is correctly strict and should stay strict;
- a future "near branch" diagnostic cannot safely reuse the exact-branch
  Hessian after merely measuring a small off-block;
- the missing broader object is a general graph-chart frontier Hessian at
  stationary, not necessarily invariant, observers.

## 1. Two-Dimensional Test Family

Take a one-dimensional visible observer in `R^2`, represented by the graph
coordinate `x` with unit vector

```text
u(x) = (1,x) / sqrt(1+x^2),
P(x) = u(x)u(x)^T.
```

For one symmetric family member

```text
A = [[a,e],
     [e,d]],
```

the weighted frontier score is

```text
F_mu(P) = ||PAP||_F^2 - mu * 1/2 ||[A,P]||_F^2.
```

At `x=0`, direct differentiation gives

```text
D F_mu(0)
= 2 e ((2+mu)a - mu d).
```

The directional second derivative is

```text
D^2 F_mu(0)
= 2 [ (1+mu)(4e^2 + 2ad - 2a^2)
       - mu(d^2-a^2) ].
```

For a finite weighted family, average these expressions with the declared
weights.

## 2. Consequence

The exact-branch condition is `e_i=0` for every family member. It implies
stationarity. But stationarity only requires

```text
sum_i w_i e_i ((2+mu)a_i - mu d_i) = 0.
```

This can hold with nonzero off-blocks.

Therefore:

```text
exact branch  => stationary,
stationary    does not imply exact branch.
```

This was easy to miss because the exact-branch sector is clean and already
module-backed.

## 3. Sharp Symmetric-Pair Breakpoint

Let

```text
A_+ = [[3,e],
       [e,1]],

A_- = [[3,-e],
       [-e,1]],
```

with equal weights and `mu=0`.

The off-blocks cancel in the first variation, so `x=0` is stationary for every
`e`. But the second derivative is

```text
D^2 F_0(0) = -24 + 8e^2.
```

If one incorrectly drops the off-block terms and applies the exact-branch
proxy, one gets the constant verdict

```text
D^2 F_proxy(0) = -24.
```

The true local sign changes at

```text
|e| = sqrt(3).
```

Thus the exact-branch proxy predicts a strict local maximum even after the
true stationary point has become locally unstable.

This is a clean breakpoint. It is not a numerical artifact and not a generic
"non-Gaussian" warning. It is a local quadratic frontier effect inside the
weighted-family geometry itself.

## 4. Numerical Check

The script

```text
audit/0_3_2_non_exact_stationary_frontier_check.py
```

writes

```text
audit/outputs/0_3_2_non_exact_stationary_frontier_check.json
```

It checks:

- the symmetric-pair threshold above;
- single-operator leakage-penalty cancellations;
- random two-operator stationary but non-exact families.

The symmetric sweep confirms the analytic threshold:

```text
e = 1.70 < sqrt(3):  true Hessian < 0, local maximum at x=0
e = sqrt(3):         true Hessian = 0, degeneracy
e = 1.75 > sqrt(3):  true Hessian > 0, x=0 no longer local maximum
```

The exact-branch proxy remains `-24` throughout, demonstrating why it must not
be used outside its exact invariance hypothesis.

## 5. Implication For 0.3.2

This note changes the near-branch research target.

Before this pass, the missing object looked like:

```text
exact-branch Hessian + small off-block residual => approximate certificate.
```

That is incomplete. The better target is:

```text
general graph-chart first variation
+ general graph-chart Hessian
+ residual/Lipschitz certificate
=> stationary or near-stationary frontier certificate.
```

The exact-branch Hessian remains a valuable closed sector, but it is not the
right parent object for all frontier critical points.

## 6. Module Boundary

No immediate module change follows from this note.

The current strictness of `exact_branch_hessian` is correct. A future API would
need a different name and different semantics, for example a declared-observer
graph Hessian or stationary-frontier diagnostic. It should not be presented as
an approximate exact-branch Hessian unless a theorem explicitly controls the
off-block second-order terms.

## Final Position

The biggest new finding is that weighted-frontier criticality is broader than
exact branch geometry even at order two. That is good news for the theory, but
it is also a warning: the project should not extend exact-branch results by
tolerance alone.

The next rigorous object is the general graph-chart Hessian at a declared
observer, with exact-branch Hessian recovered as the invariant special case.
That object is now made explicit in
`audit/0_3_2_general_graph_hessian_note.md`.
