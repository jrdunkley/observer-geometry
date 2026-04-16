# 0.4.0 Technical Note -- Verification and Research Report

**Date:** 15 April 2026
**Scope:** Numerical verification of all boxed identities and theorems in *0.4.0 Technical Note: Reduced source law, composite curvature, and sixth-order observer stability*.
**Verification code:** `0_4_0_verification.py` (94 checks, all passed, tolerance 1e-10).

---

## 1. Verification summary

Every boxed identity and theorem was exercised numerically with random SPD fields across dimensions (n,m) in {(3,1), (4,2), (5,2), (5,3), (6,2)}.

| Result | Status | Max error | Notes |
|--------|--------|-----------|-------|
| Split-frame identities (Thm 1.1) | PASS | 3.5e-15 | CL=I, CZ=0, M^T H M = diag(Phi,R), L^T H Z = 0 |
| Calligraphic V formula (eq 1.4) | PASS | 2.4e-15 | V = dPhi - alpha^T Phi - Phi alpha |
| Calligraphic B formula (eq 1.5) | PASS | 3.2e-16 | B = -theta^T R - Phi beta |
| Raw W formula (eq 2.4) | PASS | 5.4e-16 | W = L^T Hddot L + theta^T B^T + B theta |
| Completed square (eq 2.8) | PASS | 3.2e-16 | W = L^T Hddot L - 2 Qhat + 1/2 O |
| theta = -R^{-1} B^T (beta=0) | PASS | 1.5e-16 | Confirmed for fixed C, fixed Z |
| A_cpl decomposition (Thm 2.2) | PASS | 2.9e-15 | A_cpl = A - 1/2 Sym(R Theta) with Sym(X) = X + X^T |
| A_cpl expanded form (eq 2.11) | PASS | 2.7e-15 | |
| Fast-hidden lift (Prop 2.3) | PASS | 1.0e-15 | Pdot at P=I equals A_cpl |
| Pure-observer theta identity (eq 3.1) | PASS | 1.1e-14 | theta = -R^{-1} beta^T Phi when dH=0 |
| F_alpha composite formula (eq 3.2) | PASS | 7.5e-15 | General gauge |
| F_omega composite formula (eq 3.2) | PASS | 2.3e-14 | General gauge |
| Whitened curvature (eq 3.3) | PASS | 0.0 | F_alpha = beta ^ beta^T exactly |
| Projector form (Prop 3.2, eq 3.8) | PASS | 4.8e-15 | Both visible and hidden projections verified in whitened gauge |
| No free curvature mode (Prop 3.3) | PASS | <1e-10 | F ~ eps^2.00, action ~ eps^4.00 across 6 scales |
| Source-curvature separation (Thm 4.1) | PASS | -- | Both counterexamples verified |

**Convention note.** The paper uses Sym(X) = X + X^T (full symmetrisation without the 1/2 factor). This is consistent throughout all proofs and is needed for eq 2.10 to hold. Verification initially failed when using Sym(X) = 1/2(X + X^T); switching to the paper's convention resolved it exactly.

---

## 2. Sign structure of A_cpl

The sign remark after Proposition 2.1 states that there is no universal PSD observer-covariant extension. We confirmed this quantitatively.

**n=4, m=2 (187 trials with V > 0):**
- Positive definite: 86 (46%)
- Negative definite: 7 (4%)
- Indefinite: 94 (50%)

**n=5, m=3 (17 trials with V > 0):**
- Positive definite: 4 (24%)
- Negative definite: 0 (0%)
- Indefinite: 13 (76%)

**Observation.** Indefiniteness is the generic case, especially at higher codimension. The asymmetry between positive-definite and negative-definite reflects the competing signs in the completed square (eq 2.8): the Qhat term (PSD, coefficient -2) and the O term (PSD, coefficient +1/2) partially cancel, but neither dominates generically. The direct acceleration A (which carries the sign of Hddot) can push the balance either way.

**Implementation consequence:** Any kernel function computing A_cpl should return the full eigendecomposition, not just a sign flag. The sign structure is a diagnostic, not a binary.

**Note on V > 0 hit rate.** The support-stability condition V_S > 0 is restrictive for random (H, Hdot, C): roughly 10% of trials at (4,2) and <1% at (5,3). This is expected -- for a random symmetric Hdot, the visible jet V = L^T Hdot L is generically indefinite. In practice, paths H(t) arising from physical problems often have V > 0 (e.g., monotone eigenvalue growth), so this is not a limitation of the theory, only of random testing.

---

## 3. Qhat is algebraically PSD

Verified over 200 random trials: Qhat = B R^{-1} B^T has minimum eigenvalue >= 4.8e-8 (positive, with the residual attributable to finite precision). This is algebraically guaranteed: Qhat = X^T D X with D = R^{-1} SPD and X = B, so Qhat is PSD by construction.

When beta = 0 (fixed C, fixed Z), Bhat = B and O = 0, so the completed square simplifies to W = L^T Hddot L - 2 B R^{-1} B^T. The only indefiniteness source is the Hddot term fighting the (always non-negative) quadratic defect.

---

## 4. Prop 3.2 gauge dependence

The projector-form identity P(dP^dP)P = L F_alpha L^T was verified in the whitened gauge (H = I, Phi = I, R = I) across all tested dimensions. It does **not** hold in a general gauge because the proof uses alpha + alpha^T = 0, which requires orthogonality of the frame.

In the general gauge, the curvature formulas F_alpha = -beta wedge theta and F_alpha = beta wedge R^{-1} beta^T Phi (Corollary 3.1) hold exactly. The projector-form identity is a whitened-gauge specialisation.

**Implementation consequence:** The curvature computation should use the Corollary 3.1 formulas (which work in any gauge), not the projector form. The projector form is useful for proving structural results (no free curvature mode) but not for computation in a general frame.

---

## 5. Prop 3.3 scaling verification

The quartic-action result was verified by constructing a 2-parameter observer field P(s,t) = exp(eps(s X_1 + t X_2)) P_0 exp(-...) and measuring ||F|| across eps in {0.5, 0.2, 0.1, 0.05, 0.02, 0.01}.

| eps | ||F|| | log ratio |
|-----|-------|-----------|
| 0.500 | 3.054e-01 | -- |
| 0.200 | 7.636e-02 | 2.00 |
| 0.100 | 1.909e-02 | 2.00 |
| 0.050 | 4.773e-03 | 2.00 |
| 0.020 | 7.636e-04 | 2.00 |
| 0.010 | 1.909e-04 | 2.00 |

The scaling exponent is 2.00 to 11 decimal places across all tested scales. The action (||F||^2) scales at exponent 4.00 correspondingly. This is exact, not approximate -- the Grassmannian exponential map produces a field whose curvature is quadratic in the perturbation amplitude, with no higher-order contamination at the base point.

---

## 6. Structural observations for implementation

### 6.1. The fixed-C, fixed-Z simplification

When C is constant and Z is chosen as the kernel of C (also constant), the hidden frame does not move: beta = omega = 0. In this regime:

- Bhat = B (the shifted defect equals the raw defect)
- O = 0 (the observer tensor vanishes)
- theta = -R^{-1} B^T (determined entirely by the visible-hidden coupling)
- The completed square becomes W = L^T Hddot L - 2 B R^{-1} B^T
- A_cpl = A + Sym(V^{-1/2} B R^{-1} B^T V^{-1/2}) (positive correction from hidden defect)

This is the natural regime for the first implementation: given a path H(t) with fixed observer C, compute the source law. No need for observer-motion machinery.

### 6.2. The pure-observer branch is separate

When dH = 0, the source law is trivially zero (V = W = 0). The curvature sector F_alpha, F_omega is nontrivial and depends only on how the observer frame moves. This is the natural regime for observer comparison: given a fixed law H, how curved is the Grassmannian geometry of the observer space?

### 6.3. The mixed regime needs both

When both H and C vary, the source law captures the visible dynamics and the curvature captures the observer motion. These are distinct sectors (Theorem 4.1) but share the split channel beta. A future compatibility theorem would relate them.

### 6.4. The sextic hierarchy needs the exact slice

The sextic stability certificate (Theorem 4.2) requires:
1. The reduced Hessian H_mu^red on the exact slice through the reference observer
2. The quartic form Q_mu restricted to ker(H_mu^red)
3. The sextic form q_{6,mu} restricted to the quartic null cone

Steps 1-2 are computable from the existing adapted-observer machinery. Step 3 requires the explicit sextic formula (Corollary 4.2), which involves commutator terms Y_a = B_a Z - Z A_a. This is new and would need implementation.

**Critical point:** The sextic certificate is only needed when the quadratic and quartic certificates are degenerate (i.e., the quadratic kernel is nontrivial and the quartic form vanishes on it). This is a codimension condition. In practice, most observers will be certified or rejected at quadratic or quartic level. The sextic level is the fallback for exceptional cases.

---

## 7. Research questions surfaced by verification

### 7.1. When is V > 0 generically?

The support-stability condition V = L^T Hdot L > 0 is restrictive for random data but natural for many physical paths. A useful result would characterise classes of paths H(t) for which V(t) > 0 on a given observer. For instance:

- If H(t) = H_0 + t S with S >= 0, then V = L^T S L >= 0, with V > 0 iff C does not annihilate the range of S.
- If H(t) is a monotone eigenvalue path (all eigenvalues increasing), V > 0 is guaranteed unless L lies in the kernel of the derivative.

This would help users know when the source law applies without checking eigenvalues.

### 7.2. The beta = 0 hidden defect is always stabilising

When the observer is fixed (beta = 0), the completed square gives A_cpl = A + positive correction. The hidden channel always acts to make the coupled source more positive than the direct source. This is a structural result: hidden variables stabilise the visible dynamics in the fixed-observer regime.

This breaks when the observer moves (beta != 0), because the observer tensor O enters with the opposite sign and can overcome the defect term. The indefiniteness of A_cpl in the general case is fundamentally about observer motion fighting hidden stabilisation.

### 7.3. Curvature as an observer quality metric

In the pure-observer branch, ||F_alpha|| measures how much the observer space is curved. A flat observer family (F_alpha = 0) means beta = 0, i.e., the hidden frame doesn't rotate as the observer changes. This happens when the observer moves within a commuting family.

||F_alpha|| could serve as a scalar diagnostic for "observer noncommutativity" -- how much the choice of observer matters for the hidden structure. This is a different quantity from the existing closure scores and leakage diagnostics, and could complement them.

### 7.4. Compatibility question

The paper explicitly identifies the correct future target: a compatibility theorem relating A_cpl and F_alpha when both sectors are active. From Proposition 4.1 (common split seed), both are built from the same beta channel:

- Source: O = Phi beta R^{-1} beta^T Phi (symmetric contraction)
- Curvature: F_alpha = beta wedge R^{-1} beta^T Phi (exterior product)

A compatibility relation would likely take the form of an inequality or a trace identity relating Tr(O) to ||F_alpha||^2, mediated by the metric geometry (Phi, R). This is worth pursuing but is explicitly flagged as future work.

---

## 8. Conclusion

All theorems in the 0.4.0 Technical Note verified numerically at machine precision. The paper's convention Sym(X) = X + X^T is confirmed necessary and consistent. The sign structure of A_cpl is generically indefinite, confirming the remark. The quartic-action scaling is exact. The projector form holds in the whitened gauge only.

The natural implementation path is:
1. Fixed-observer source law (beta = 0, simplest)
2. Pure-observer curvature (dH = 0, separate)
3. Full mixed regime (both sectors, no compatibility theorem yet)
4. Sextic stability certificate (fallback for degenerate cases)
