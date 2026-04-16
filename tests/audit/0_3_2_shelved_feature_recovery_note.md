# 0.3.2 Shelved Feature Recovery Note

Date: 2026-04-11

Status: pure research note. This is not an implementation plan and does not
modify module source.

## Executive Result

Several ideas that were correctly withheld in older builds now have safe
semantics, provided they are named narrowly:

- rank-k covariance/Fisher perturbation classification;
- declared-ladder observer comparison with dimension cost;
- guarded fibre-dominance diagnostics;
- affine-hidden stage-sign and branch-reversal diagnostics;
- residual/ensemble finite-candidate margin wrappers;
- declared local frontier certificates from the general graph Hessian.

The shared rule is:

```text
declared finite object + exact local/special-sector theorem + explicit scope
```

not:

```text
generic non-Gaussian branch selector or global observer optimizer.
```

The companion script is:

```text
audit/0_3_2_shelved_feature_recovery_check.py
audit/outputs/0_3_2_shelved_feature_recovery_check.json
```

## Readiness Table

| Candidate | Status | Reason |
|---|---|---|
| Rank-k covariance/Fisher perturbation | Ready after theorem note | Exact Woodbury formula; hidden-gap increment is difference of two rank-k PSD terms, rank at most `2k`. |
| Declared-ladder dimension-cost intervals | Ready after theorem note | Exact finite comparison of supplied scores `s_j - c d_j`; no global optimizer claim. |
| General graph-chart frontier Hessian | Ready for theorem note | Exact first/Hessian formulas; frame invariance checked; exact-branch Hessian recovered as invariant sector. |
| Declared local frontier certificate | Short theorem note, then maybe API planning | Exact abstract certificate, but `L_cert` is conservative and can be vacuous. |
| Guarded fibre-dominance diagnostic | Diagnostic-only, low risk | Useful only as raw centered norms plus optional ratio under denominator floor and declared sample measure. |
| Affine-hidden stage sign / branch reversal | Ready after caveated theorem note | Exact in affine-hidden Gaussian-fibre sector; needs measure/coordinate convention. |
| Residual-margin finite candidate wrapper | Ready | Exact supplied-bound theorem; already partially surfaced. |
| Ensemble finite-candidate wrapper | Short theorem note | Exact once sampling/confidence bounds are supplied or validly computed; not a cumulant engine. |
| Generic non-Gaussian branch probabilities | Still blocked | Needs full-law/higher-order/support data. |
| Global noncommuting Grassmannian optimizer | Still blocked | Needs optimizer policy and proof. |
| Support/restart probability-event engine | Still blocked as probability law | Current support strata are PSD rank/kernel events, not probability atoms. |

## 1. Rank-k Perturbation

For

```text
Sigma_1 = Sigma_0 + F F^T,   F in R^{n x k},
```

Woodbury gives

```text
H_1 = H_0 - H_0 F (I + F^T H_0 F)^{-1} F^T H_0.
```

On a coordinate visible block, the visible covariance block also updates by
`F_V F_V^T`, so

```text
Phi_1 = Phi_0
  - Phi_0 F_V (I + F_V^T Phi_0 F_V)^{-1} F_V^T Phi_0.
```

Therefore the visible hidden-gap increment is

```text
Delta gap =
  Phi_0 F_V Gamma F_V^T Phi_0
  - (H_0 F)_V Beta (H_0 F)_V^T,
```

where

```text
Gamma = (I + F_V^T Phi_0 F_V)^{-1},
Beta  = (I + F^T H_0 F)^{-1}.
```

This is a difference of two rank-k PSD terms, hence rank at most `2k`.

The probe checked `k=1,2,3` with residuals:

```text
2.36e-16, 1.53e-16, 6.94e-17.
```

Ranks were bounded by `2k` in every case. This is theorem-grade.

## 2. Declared-Ladder Dimension Cost

For a finite declared ladder with scores `s_j` and dimensions `d_j`, compare

```text
score_j(c) = s_j - c d_j.
```

Pairwise crossings are exact:

```text
c_{ab} = (s_a - s_b)/(d_a - d_b),
```

when `d_a != d_b`.

The probe used:

```text
small_clean: score 4.0, dim 1
medium:      score 5.8, dim 2
large:       score 6.4, dim 4
```

The winner changed:

```text
c < 0.3:        large
0.3 < c < 1.8:  medium
c > 1.8:        small_clean
```

The extra crossing at `0.8` is real but not winner-changing. This shows why a
phase diagram should report all pairwise crossings and the upper envelope.

This is safe if and only if the ladder is declared. It is not global
Grassmannian optimization.

## 3. Guarded Fibre Dominance

A naked ratio

```text
||fibre|| / ||variational||
```

is unsafe when the denominator is small. The safe diagnostic is:

```text
fibre_centered_norm,
variational_centered_norm,
ratio if variational_centered_norm >= floor else undefined,
declared sample measure / grid.
```

The probe showed the intended behavior. With flat variational action and zero
or tiny coupling:

```text
ratio_defined = false
```

because the denominator was `0` or `3.07e-08` under a `1e-6` floor. With
larger coupling the ratio became defined, and with a quartic base action it
was defined throughout.

This should remain diagnostic-only unless the sample measure and denominator
guard are part of the theorem/API contract.

## 4. Affine-Hidden Sign And Branch Reversal

The affine-hidden sector is exact:

```text
p(v,h) proportional to exp(-A(v) - 1/2 h^T D(v)h - J(v)^T h).
```

The visible action is

```text
A(v) + 1/2 log det D(v) - 1/2 J(v)^T D(v)^{-1}J(v),
```

up to a visible-independent additive constant.

The probe reproduced a branch flip where the variational action is flat:

```text
variational_action = [0, 0],
fibre_volume       = [-1.1513, 1.1513],
visible_action     = [-1.1513, 1.1513].
```

It also checked stage signs for `D_ee=4`:

```text
j0=0: action_shift =  0.6931
j0=1: action_shift =  0.5681
j0=3: action_shift = -0.4319
```

matching

```text
1/2 log det D_ee - 1/2 J_e^T D_ee^{-1}J_e.
```

This is ready for caveated theorem writing, not generic marginalization.

## 5. Ensemble And Residual Margins

The residual-margin theorem is safe:

```text
quadratic_gap > R_a + R_b
```

certifies the ordering when residuals are bounded by `R_a,R_b`.

For ensemble means, sampling error is just another residual once the sampling
bound is valid. The probe used Hoeffding bounds on bounded scores:

- `plausible_but_not_certified`: empirical winner `A`, but gaps did not beat
  the Hoeffding residual, so no robust certificate.
- `large_gap_certified`: empirical winner `A`, large sample and large gaps;
  all pairwise comparisons certified.

This is exactly the desired behavior. It refuses seductive but underpowered
empirical gaps.

## 6. Declared Local Frontier Certificate

The general graph Hessian work changes the feature boundary. A future surface
should not be an approximate exact-branch Hessian. It should be a declared
local frontier certificate:

```text
declared observer
stationarity residual
general graph Hessian spectrum
conservative Lipschitz radius
pass/fail sufficient certificate
```

This is broader than exact branches because stationary non-exact observers
exist. It is still local and sufficient-only.

## 7. Stack Handoff Implications

A quick scan of `nomodescent` and `evidence` did not reveal a reason to widen
their epistemic claims. They are already disciplined:

- `nomodescent` separates exact descent from audited approximate PSD search;
- `evidence` preserves exact, inferred, and ambiguous source material before
  downstream assembly;
- finite-family observer suggestion in `evidence` is already labelled as a
  deterministic ranking aid, not a scientific selector.

The recovered features could still buy something in the broader stack later:

- declared-ladder dimension-cost intervals can improve finite observer-family
  ranking without claiming global optimality;
- residual/ensemble finite-candidate margins can turn some evidence-level
  score comparisons into honest sufficient certificates;
- guarded fibre-dominance can be a reportable diagnostic when affine-hidden
  fibre terms are explicitly supplied;
- declared local frontier certificates could help `nomodescent` distinguish
  "local declared observer verdict" from "global common-descent result."

But these are integration opportunities, not immediate claim expansions. The
claim hierarchy should remain unchanged.

## Planning Recommendation

For the next planning phase, split recovered work into:

1. theorem/writeup first:
   rank-k perturbation, declared-ladder intervals, affine-hidden sign,
   general graph Hessian;
2. possible module candidates after theorem text:
   declared local frontier certificate, guarded fibre dominance, ensemble
   finite-candidate certification;
3. stack integration candidates after kernel semantics are stable:
   declared-ladder ranking handoff to `evidence`, residual-margin handoff from
   evidence bundles to declared finite comparisons, and local frontier
   certificate handoff to `nomodescent` audit objects;
4. keep blocked:
   full non-Gaussian branch probabilities, global noncommuting optimizer,
   probability support/restart event engine from Hessians alone.

## Final Position

There is real low-hanging fruit. Most of it is not new mathematics but safer
naming and theorem packaging. The strongest immediate module-shaped candidate
is the declared local frontier certificate, but only after the general graph
Hessian theorem is written and the conservative `L_cert` semantics are
accepted.
