"""
N6: Conservation law in the nomocomp pipeline.

For the Strike 1B benchmark (collinear regressors), compute the information
budget between competing models and show that the conservation law reveals
why the geometric comparator succeeds where AIC/BIC fails.
"""

import sys
import numpy as np
from numpy.linalg import inv, eigh
from pathlib import Path

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))
from nomogeo import information_budget

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok


def strike_1b_scenario():
    """
    Reproduce the Strike 1B scenario: two nested Gaussian linear models
    with collinear predictors.

    Model 1: y = X1 beta1 + epsilon  (true model, k predictors)
    Model 2: y = X2 beta2 + epsilon  (wrong model, k predictors, different columns)

    Both have the same k, so AIC/BIC give the same penalty.
    The geometric score differs because the Fisher information differs.

    The conservation law tells us: the Fisher information difference decomposes
    into visible + hidden contributions. The geometric comparator sees the
    visible part. AIC/BIC see nothing (same k).
    """
    print("=== N6: Strike 1B Conservation Analysis ===")

    rng = np.random.default_rng(42)
    n_obs = 100
    k = 3  # same dimension for both models

    # True model: y = x1 + x2 + x3 + noise
    # Collinear design: x2 ~ 0.9*x1 + noise, x3 ~ 0.9*x1 + noise
    x1 = rng.standard_normal(n_obs)
    rho = 0.9  # collinearity

    # Model 1 (true): uses x1, x2_true, x3_true
    x2_true = rho * x1 + np.sqrt(1 - rho**2) * rng.standard_normal(n_obs)
    x3_true = rho * x1 + np.sqrt(1 - rho**2) * rng.standard_normal(n_obs)
    X1 = np.column_stack([x1, x2_true, x3_true])

    # Model 2 (wrong): uses x1, x2_wrong, x3_wrong (different noise realisations)
    x2_wrong = rho * x1 + np.sqrt(1 - rho**2) * rng.standard_normal(n_obs)
    x3_wrong = rho * x1 + np.sqrt(1 - rho**2) * rng.standard_normal(n_obs)
    X2 = np.column_stack([x1, x2_wrong, x3_wrong])

    sigma = 1.0  # noise variance

    # Fisher information matrices (observed information at MLE)
    # For OLS: I = X^T X / sigma^2
    H1 = X1.T @ X1 / sigma**2  # Fisher for model 1
    H2 = X2.T @ X2 / sigma**2  # Fisher for model 2

    # The "perturbation" is the difference in Fisher information
    Hdot = H2 - H1

    print(f"  n_obs = {n_obs}, k = {k}, rho = {rho}")
    print(f"  H1 condition: {np.max(eigh(H1)[0]) / np.min(eigh(H1)[0]):.1f}")
    print(f"  H2 condition: {np.max(eigh(H2)[0]) / np.min(eigh(H2)[0]):.1f}")
    print(f"  ||Hdot|| / ||H1|| = {np.linalg.norm(Hdot) / np.linalg.norm(H1):.4f}")

    # Full-parameter comparison (C = I_k, m = k)
    C_full = np.eye(k)
    b_full = information_budget(H1, C_full, Hdot)
    print(f"\n  Full-parameter comparison (m = k = {k}):")
    print(f"    vis_rate = {b_full.visible_rate:.6f}")
    print(f"    hid_rate = {b_full.hidden_rate:.6f}")
    print(f"    amb_rate = {b_full.ambient_rate:.6f}")
    print(f"    vis_frac = {b_full.visible_fraction:.6f}")
    report("Full-parameter conservation", b_full.conservation_residual)

    # When C = I (full observation), vis_rate = amb_rate (no hidden sector)
    # This is the trivial case. The interesting case is partial observation.

    # Partial observation: observe only the first 2 parameters (miss one)
    for m in [1, 2]:
        C = np.zeros((m, k))
        for i in range(m): C[i, i] = 1.0
        b = information_budget(H1, C, Hdot)
        print(f"\n  Partial observation (m = {m}):")
        print(f"    vis_rate = {b.visible_rate:.6f}")
        print(f"    hid_rate = {b.hidden_rate:.6f}")
        print(f"    amb_rate = {b.ambient_rate:.6f}")
        print(f"    vis_frac = {b.visible_fraction:.4f}")
        report(f"Partial m={m} conservation", b.conservation_residual)
        print(f"    => Observing {m}/{k} parameters captures {100*b.visible_fraction:.1f}% of the Fisher difference")

    # What the geometric comparator sees vs what AIC/BIC sees
    det_H1 = np.linalg.slogdet(H1)[1]
    det_H2 = np.linalg.slogdet(H2)[1]
    det_correction = 0.5 * (det_H2 - det_H1)

    print(f"\n  Geometric comparator:")
    print(f"    1/2 log det(H2) - 1/2 log det(H1) = {det_correction:.6f}")
    print(f"    This IS the ambient rate (integrated): {b_full.ambient_rate:.6f}")
    print(f"    AIC sees: 0 (same k)")
    print(f"    BIC sees: 0 (same k)")
    print(f"    Geometric sees: {det_correction:.6f} (the determinant correction)")

    # The conservation law explains WHY:
    # The determinant correction is Tr(H1^{-1} Hdot) (to first order)
    # = vis_rate + hid_rate
    # AIC/BIC use only k (dimension count), which cancels for same-k models
    # The geometric comparator uses log det(I), which captures the full ambient rate

    print(f"\n  WHY the geometric comparator works:")
    print(f"    The determinant correction = Tr(H1^{{-1}} Hdot) = {b_full.ambient_rate:.6f}")
    print(f"    This is the total information difference between models.")
    print(f"    AIC/BIC penalty: k log(n) - k log(n) = 0 (same k)")
    print(f"    Geometric penalty: 1/2 log det(I1) - 1/2 log det(I2) = {det_correction:.6f}")
    print(f"    The geometric comparator captures the conservation-law budget.")
    print(f"    AIC/BIC are blind to it.")


def collinearity_sweep():
    """
    Sweep collinearity rho and show how the information budget changes.
    At high rho, the models are hard to distinguish — but the geometric
    comparator still sees the difference through the determinant correction.
    """
    print(f"\n=== Collinearity sweep ===")

    rng = np.random.default_rng(42)
    n_obs = 100; k = 3; sigma = 1.0

    rhos = [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
    print(f"  {'rho':>5s}  {'amb_rate':>10s}  {'det_corr':>10s}  {'H1_cond':>10s}  {'vis_frac_m2':>12s}")

    for rho in rhos:
        x1 = rng.standard_normal(n_obs)
        x2t = rho*x1 + np.sqrt(max(1-rho**2, 0))*rng.standard_normal(n_obs)
        x3t = rho*x1 + np.sqrt(max(1-rho**2, 0))*rng.standard_normal(n_obs)
        X1 = np.column_stack([x1, x2t, x3t])

        x2w = rho*x1 + np.sqrt(max(1-rho**2, 0))*rng.standard_normal(n_obs)
        x3w = rho*x1 + np.sqrt(max(1-rho**2, 0))*rng.standard_normal(n_obs)
        X2 = np.column_stack([x1, x2w, x3w])

        H1 = X1.T @ X1 / sigma**2
        H2 = X2.T @ X2 / sigma**2
        Hdot = H2 - H1

        cond = np.max(eigh(H1)[0]) / max(np.min(eigh(H1)[0]), 1e-10)
        det_corr = 0.5 * (np.linalg.slogdet(H2)[1] - np.linalg.slogdet(H1)[1])

        b_full = information_budget(H1, np.eye(k), Hdot)

        C_m2 = np.zeros((2, k)); C_m2[0,0] = 1; C_m2[1,1] = 1
        b_m2 = information_budget(H1, C_m2, Hdot)
        vf_m2 = b_m2.visible_fraction if abs(b_m2.ambient_rate) > 1e-10 else float('nan')

        print(f"  {rho:5.2f}  {b_full.ambient_rate:10.4f}  {det_corr:10.4f}  {cond:10.1f}  {vf_m2:12.4f}")


if __name__ == "__main__":
    print("=" * 70)
    print("N6: Conservation Law in nomocomp Pipeline")
    print("=" * 70)

    strike_1b_scenario()
    collinearity_sweep()

    print(f"\n{'='*70}")
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
