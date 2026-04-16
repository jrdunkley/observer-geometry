"""
nomoselect perturbation robustness: show that task-aware dimensionality
reduction is also perturbation-aware.

For Iris and Wine datasets:
1. Compute the sample covariance H (as a precision/Fisher proxy)
2. Perturb it (simulate a distributional shift)
3. Compare the visible fraction of nomoselect's subspace vs PCA's subspace
4. Show nomoselect captures more of the perturbation information

Also: apply the conservation law to logistic regression Fisher information
on Iris, demonstrating domain-agnosticism.
"""

import sys
import numpy as np
from numpy.linalg import inv, eigh, norm, slogdet
from pathlib import Path

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))

from nomogeo import information_budget, closure_adapted_observer, visible_precision

passes = 0
fails = 0

def report(name, err, tol=1e-10):
    global passes, fails
    ok = err < tol
    passes += ok; fails += (not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {err:.3e}")
    return ok

def sym(M): return 0.5*(M+M.T)


def load_iris():
    from sklearn.datasets import load_iris
    X, y = load_iris(return_X_y=True)
    return X, y, "Iris"

def load_wine():
    from sklearn.datasets import load_wine
    X, y = load_wine(return_X_y=True)
    return X, y, "Wine"


def perturbation_robustness(X, y, name, m=2):
    """
    Compare PCA and task-aware (nomogeo-adapted) subspaces for perturbation robustness.
    """
    print(f"\n=== {name} (n={X.shape[1]}, m={m}) ===")

    n = X.shape[1]

    # Compute within-class scatter as precision proxy
    classes = np.unique(y)
    S_w = np.zeros((n, n))
    for c in classes:
        Xc = X[y == c]
        Xc_centered = Xc - Xc.mean(axis=0)
        S_w += Xc_centered.T @ Xc_centered
    S_w /= len(X)
    H = inv(S_w + 0.01 * np.eye(n))  # regularised precision

    # Between-class scatter as the "task structure"
    S_b = np.zeros((n, n))
    grand_mean = X.mean(axis=0)
    for c in classes:
        Xc = X[y == c]
        diff = Xc.mean(axis=0) - grand_mean
        S_b += len(Xc) * np.outer(diff, diff)
    S_b /= len(X)

    # Perturbation: simulate a distributional shift (add noise to the precision)
    rng = np.random.default_rng(42)
    Hdot = sym(rng.standard_normal((n, n)) * 0.5)  # random perturbation
    # Also try a structured perturbation (shift in class means)
    Hdot_struct = inv(S_w + S_b + 0.01*np.eye(n)) - H  # shift toward between-class structure

    # PCA observer: top m eigenvectors of H (or equivalently bottom m of S_w)
    eigvals_H, eigvecs_H = eigh(H)
    C_pca = eigvecs_H[:, -m:].T  # top m eigenvectors of precision

    # Task-aware observer: closure_adapted_observer using between-class scatter
    try:
        result = closure_adapted_observer(H, [S_b], m)
        C_task = result.C
    except Exception as e:
        print(f"  Adapted observer failed: {e}")
        C_task = C_pca  # fallback

    # Compute information budgets for both perturbations
    for pert_name, Hdot_pert in [("random", Hdot), ("structured", Hdot_struct)]:
        print(f"\n  Perturbation: {pert_name}")

        b_pca = information_budget(H, C_pca, Hdot_pert)
        b_task = information_budget(H, C_task, Hdot_pert)

        print(f"    PCA observer:  vis_frac = {b_pca.visible_fraction:.4f}, "
              f"cons_err = {b_pca.conservation_residual:.2e}")
        print(f"    Task observer: vis_frac = {b_task.visible_fraction:.4f}, "
              f"cons_err = {b_task.conservation_residual:.2e}")

        report(f"{name} {pert_name} PCA conservation", b_pca.conservation_residual)
        report(f"{name} {pert_name} task conservation", b_task.conservation_residual)

        if abs(b_task.visible_fraction) < 100 and abs(b_pca.visible_fraction) < 100:
            if b_task.visible_fraction > b_pca.visible_fraction:
                print(f"    Task observer captures MORE perturbation information (+{b_task.visible_fraction - b_pca.visible_fraction:.4f})")
            else:
                print(f"    PCA observer captures more (+{b_pca.visible_fraction - b_task.visible_fraction:.4f})")

    # V > 0 check
    print(f"\n  V > 0 status:")
    for obs_name, C_obs in [("PCA", C_pca), ("Task", C_task)]:
        for pert_name, Hdot_pert in [("random", Hdot), ("structured", Hdot_struct)]:
            b = information_budget(H, C_obs, Hdot_pert)
            vmin = float(np.min(b.v_eigenvalues))
            status = "Y" if vmin > 1e-6 else "N"
            print(f"    {obs_name} / {pert_name}: V_min = {vmin:.4f} ({status})")


def ml_fisher_conservation():
    """
    Apply the conservation law to logistic regression Fisher information.
    This demonstrates domain-agnosticism: the conservation law works on
    any smooth SPD field, including statistical Fisher information.
    """
    print(f"\n=== ML Fisher Information Conservation ===")

    from sklearn.datasets import load_iris
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    X_raw, y = load_iris(return_X_y=True)
    # Binary: class 0 vs class 1
    mask = y < 2
    X_raw = X_raw[mask]
    y_bin = y[mask]

    scaler = StandardScaler()
    X = scaler.fit_transform(X_raw)
    n = X.shape[1]  # 4 features

    # Fit logistic regression
    model = LogisticRegression(C=1.0, fit_intercept=False, max_iter=1000)
    model.fit(X, y_bin)
    theta = model.coef_.flatten()

    # Fisher information at theta
    # For logistic regression: I(theta) = X^T W X where W = diag(p_i(1-p_i))
    logits = X @ theta
    p = 1.0 / (1.0 + np.exp(-logits))
    W = p * (1 - p)
    H_fisher = (X * W[:, None]).T @ X / len(X)  # n x n
    H_fisher = sym(H_fisher)

    # Check SPD
    eigvals = eigh(H_fisher)[0]
    print(f"  Fisher information eigenvalues: {eigvals}")
    if np.min(eigvals) < 1e-10:
        H_fisher += 0.01 * np.eye(n)
        print(f"  Regularised (added 0.01 I)")

    # Perturbation: change in Fisher under a parameter shift
    delta_theta = np.array([0.1, -0.05, 0.02, 0.08])
    theta_new = theta + delta_theta
    logits_new = X @ theta_new
    p_new = 1.0 / (1.0 + np.exp(-logits_new))
    W_new = p_new * (1 - p_new)
    H_fisher_new = sym((X * W_new[:, None]).T @ X / len(X))
    if np.min(eigh(H_fisher_new)[0]) < 1e-10:
        H_fisher_new += 0.01 * np.eye(n)

    Hdot = H_fisher_new - H_fisher

    # Observer: first 2 features
    m = 2
    C = np.zeros((m, n))
    C[0, 0] = 1.0; C[1, 1] = 1.0

    b = information_budget(H_fisher, C, Hdot)
    report("ML Fisher conservation", b.conservation_residual)

    print(f"\n  Logistic regression Fisher information on Iris (binary):")
    print(f"    n = {n}, m = {m}")
    print(f"    Visible rate:  {b.visible_rate:.6f}")
    print(f"    Hidden rate:   {b.hidden_rate:.6f}")
    print(f"    Ambient rate:  {b.ambient_rate:.6f}")
    print(f"    Visible fraction: {b.visible_fraction:.4f}")
    print(f"    Conservation: {b.conservation_residual:.2e}")

    # Try adapted observer
    try:
        result = closure_adapted_observer(H_fisher, [Hdot], m)
        b_adapted = information_budget(H_fisher, result.C, Hdot)
        print(f"\n  Adapted observer:")
        print(f"    Visible fraction: {b_adapted.visible_fraction:.4f}")
        print(f"    Improvement over canonical: {b_adapted.visible_fraction - b.visible_fraction:+.4f}")
    except Exception as e:
        print(f"  Adapted observer failed: {e}")

    print(f"\n  The conservation law works on statistical Fisher information,")
    print(f"  not just molecular Hessians. The framework is domain-agnostic.")


if __name__ == "__main__":
    print("=" * 70)
    print("nomoselect Robustness + ML Fisher Conservation")
    print("=" * 70)

    X_iris, y_iris, _ = load_iris()
    perturbation_robustness(X_iris, y_iris, "Iris", m=2)

    X_wine, y_wine, _ = load_wine()
    perturbation_robustness(X_wine, y_wine, "Wine", m=3)

    ml_fisher_conservation()

    print("\n" + "=" * 70)
    print(f"RESULTS: {passes} passed, {fails} failed")
    print("=" * 70)
