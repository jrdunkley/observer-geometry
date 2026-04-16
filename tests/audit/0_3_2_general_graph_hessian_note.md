# 0.3.2 General Graph-Chart Frontier Hessian Note

Date: 2026-04-11

Status: pure research note. This is not an implementation plan and does not
modify module source.

## Executive Result

The right order-two parent object for weighted-frontier criticality is the
general graph-chart Hessian at a declared observer. The existing
`exact_branch_hessian` is recovered exactly when all off-blocks vanish, but it
does not cover stationary non-exact observers.

This note records the theorem target and the numerical checks in

```text
audit/0_3_2_general_graph_hessian_check.py
audit/outputs/0_3_2_general_graph_hessian_check.json
audit/0_3_2_graph_hessian_invariance_check.py
audit/outputs/0_3_2_graph_hessian_invariance_check.json
```

## Setup

Use a fixed splitting

```text
R^n = U plus W,
dim U = m,
```

and graph-chart coordinates

```text
Y(X) = [I; X],
P(X) = Y(X)(I + X^T X)^{-1}Y(X)^T,
X in Hom(U,W).
```

In an arbitrary orthonormal observer frame `B` for `U` and complement frame
`N` for `W`, the same chart is

```text
span(B + N X).
```

Changing frames by

```text
B' = B R_U,
N' = N R_W
```

changes coordinates by

```text
X' = R_W^T X R_U,
```

but not the represented tangent direction. This frame contract is important:
Hessian matrices are only directly comparable after the complement and visible
bases are fixed.

For each symmetric family member write the block decomposition in the chosen
observer/complement frame:

```text
K_i = [[U_i, E_i^T],
       [E_i, D_i  ]].
```

The weighted frontier is

```text
F_mu(P) = sum_i w_i [
  ||P K_i P||_F^2
  - mu * 1/2 ||[K_i,P]||_F^2
].
```

## First Variation

The exact first variation at the declared observer is

```text
DF_mu(0)[X] = <G, X>_F,
```

with gradient

```text
G = 2 sum_i w_i [
      (2+mu) E_i U_i
      - mu D_i E_i
    ].
```

Thus exact branch, `E_i=0` for all `i`, implies stationarity, but
stationarity can also occur by cancellation in the weighted sum.

## Directional Second Variation

For a tangent direction `X`, the exact directional second variation at the
declared observer `X=0` is

```text
D^2 F_mu(0)[X,X]
= 2 sum_i w_i [
    (1+mu) S_{i,2}(X)
    - mu T_{i,2}(X)
  ],
```

where

```text
B_i(X) = E_i^T X + X^T E_i,

S_{i,2}(X)
= Tr(B_i(X)^2)
 + 2 Tr(U_i X^T D_i X)
 - 2 Tr((X^T X) U_i^2),

T_{i,2}(X)
= Tr(X^T (K_i^2)_{WW} X)
 - Tr((X^T X) (K_i^2)_{UU}).
```

Equivalently, using the chosen block notation,

```text
(K_i^2)_{UU} = U_i^2 + E_i^T E_i,
(K_i^2)_{WW} = D_i^2 + E_i E_i^T.
```

This formula is exact whether or not the observer is an exact branch.

The associated Hessian bilinear form is obtained by polarization:

```text
H(X,Y)
= 1/4 [
  D^2F(0)[X+Y,X+Y]
  - D^2F(0)[X-Y,X-Y]
].
```

This is the clean theorem target for a future general declared-observer
frontier Hessian.

## Frame Invariance Checks

The script

```text
audit/0_3_2_graph_hessian_invariance_check.py
```

checks that the formulas are geometric rather than coordinate artifacts.

For a random embedded observer and complement frame, the first variation
matched finite differences with

```text
max_abs_gradient_formula_minus_finite = 9.52e-09,
gradient_relative_fro_error           = 6.55e-10.
```

The bilinear Hessian matched finite differences with

```text
max_abs_formula_minus_finite = 3.64e-05,
max_relative_error           = 4.62e-06.
```

Changing visible and hidden orthonormal frames gave

```text
max_abs_bilinear_invariance_error = 5.68e-13.
```

For arbitrary-frame exact-branch reduction, comparing the raw Hessian matrix
in a declared complement frame against the module's result can be misleading,
because the module computes its own complement basis. The raw matrix mismatch
in the check was large:

```text
raw_matrix_error_in_different_complement_frames = 18.92.
```

But after expressing the formula in the module's complement frame, the
agreement returned to roundoff:

```text
max_abs_formula_minus_module_in_module_complement_frame = 1.42e-14,
max_abs_eigenvalue_error                                = 4.26e-14.
```

This is the right interpretation: the bilinear form and spectrum are
coordinate-invariant; matrix entries are not.

## Reduction To Exact-Branch Hessian

When `E_i=0` for every family member, the formula reduces to the existing
exact-branch sector. The check script compared the polarized graph Hessian
against the module's `exact_branch_hessian(...).second_variation_operator` in
a random block-diagonal family and found

```text
max_abs_formula_minus_module = 3.55e-15.
```

So the general graph Hessian is not a competing object. It contains the
existing exact-branch Hessian as the invariant/off-block-zero specialization.

## Non-Exact Check

For a random non-exact family, the polarized formula was compared against a
finite-difference Hessian in graph coordinates. The check found

```text
max_abs_formula_minus_finite = 1.80e-05,
relative_fro_error           = 3.35e-07.
```

This is consistent with the finite-difference step size and confirms that the
formula remains valid outside the exact-branch sector.

## Sharp Sign Breakpoint Recovered

For the symmetric pair

```text
A_+ = [[3,e],[ e,1]],
A_- = [[3,-e],[-e,1]],
mu = 0,
```

the general graph Hessian gives

```text
e = 1:       H = -16
e = sqrt(3): H = 0
e = 2:       H = 8
```

while the invalid exact-branch proxy remains

```text
H_proxy = -24.
```

This recovers the breakpoint from
`audit/0_3_2_non_exact_stationary_frontier_note.md`.

## Implication

The near-branch research program should now be reframed:

```text
general graph Hessian
+ stationarity residual
+ Hessian-Lipschitz certificate
=> local declared-observer verdict
```

This reframing is developed in
`audit/0_3_2_declared_frontier_certificate_note.md`.

Exact-branch Hessian is a closed, strict, useful sector inside that theory.
It should not be softened by tolerance into a generic approximate branch
Hessian.

## Final Position

This is the most concrete next theoretical object for 0.3.2. It is exact
local quadratic geometry, not a full non-Gaussian law result and not a
Grassmannian optimizer. It is nevertheless broader than exact-branch geometry
and explains the stationary non-exact edge cases.
