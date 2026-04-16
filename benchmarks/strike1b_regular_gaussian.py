"""
Strike 1B: Regular-Gaussian benchmark — AIC/BIC vs exact local evidence.

Proves that the omitted Laplace determinant correction can reverse model
ranking in a concrete nested Gaussian family fitted via statsmodels OLS.

Benchmark family
----------------
    n = 60, seed = 1
    x₁ ~ N(0,1),  z ~ N(0,1),  x₂ = 0.98·x₁ + √(1-0.98²)·z
    y = 1 + 2·x₁ + 0.1·x₂ + ε,   ε ~ N(0, 2²)

    M₁: intercept + x₁          (k₁ = 3: β₀, β₁, σ²)
    M₂: intercept + x₁ + x₂    (k₂ = 4: β₀, β₁, β₂, σ²)

Comparison scores (all in log-evidence units, higher = better)
--------------------------------------------------------------
    AIC-style:   log L - k  =  -n·S_vis - k
                 (from AIC = -2 log L + 2k; maximisation form = -AIC/2)
    BIC-style:   -n·S_vis - (k/2)·log(n)
    Exact local: -n·S_vis - (k/2)·log(n) - ½·log det(J) + (k/2)·log(2π)

where S_vis = NLL/n (per-sample negative log-likelihood at the MLE)
and J is the per-sample observed information (Hessian of S_vis).

The determinant correction is  Δ = -½·log det(J₂) + ½·log det(J₁)
                                    + (k₂/2)·log(2π) - (k₁/2)·log(2π)

Scope conditions:
    - Regular Gaussian nested family, fitted by OLS.
    - Flat improper prior on (β, σ²).  For nested models this is
      standard; the prior on shared parameters cancels in comparison.
    - "Exact local Laplace" = exact to quadratic order at the MLE.
      Exact for β given σ²; approximate for σ² (inverse-χ² posterior).
"""
from __future__ import annotations

import numpy as np
import statsmodels.api as sm
from dataclasses import dataclass


@dataclass
class ModelScore:
    """All comparison scores for one fitted model."""
    name: str
    k: int                    # active dimension (num params including σ²)
    nll: float                # total negative log-likelihood at MLE
    S_vis: float              # per-sample NLL = nll / n
    log_det_J: float          # log det of per-sample observed information
    aic_score: float          # -n·S_vis - k  = log L - k  (higher = better)
    bic_score: float          # -n·S_vis - (k/2)·log(n)
    exact_score: float        # -n·S_vis - (k/2)·log(n) - ½·log det(J) + (k/2)·log(2π)


def generate_data(n: int = 60, seed: int = 1, rho: float = 0.98):
    """Generate the benchmark dataset."""
    rng = np.random.default_rng(seed)
    x1 = rng.standard_normal(n)
    z = rng.standard_normal(n)
    x2 = rho * x1 + np.sqrt(1 - rho**2) * z
    eps = 2.0 * rng.standard_normal(n)
    y = 1.0 + 2.0 * x1 + 0.1 * x2 + eps
    return y, x1, x2


def fit_and_score(y, X, name: str, n: int) -> ModelScore:
    """Fit OLS and compute all comparison scores.

    The per-sample observed information J is the Hessian of S_vis(θ)
    with respect to the full parameter vector θ = (β, σ²), evaluated
    at the MLE.  For a Gaussian linear model with known structure,
    the MLE Hessian in the (β, σ²) chart is block-diagonal:

        H_ββ = X'X / (n σ²)        [k_β × k_β]
        H_σσ = 1 / (2 σ⁴)          [1 × 1]
        H_βσ = 0                    [off-diagonal]

    This is the per-sample Hessian (divided by n where needed).
    """
    model = sm.OLS(y, X)
    res = model.fit()

    k_beta = X.shape[1]          # number of regression coefficients
    k = k_beta + 1               # +1 for σ²
    sigma2_mle = np.sum(res.resid**2) / n  # MLE (not OLS) variance

    # Per-sample negative log-likelihood at MLE:
    # NLL = (n/2)·log(2π) + (n/2)·log(σ²) + n/2
    nll = 0.5 * n * np.log(2 * np.pi) + 0.5 * n * np.log(sigma2_mle) + 0.5 * n
    S_vis = nll / n

    # Per-sample observed information matrix in (β, σ²) chart.
    # H_ββ = X'X / (n · σ²_mle)
    XtX = X.T @ X
    H_bb = XtX / (n * sigma2_mle)

    # H_σσ = 1 / (2 σ⁴_mle)  — Hessian of S_vis w.r.t. σ²
    H_ss = 1.0 / (2.0 * sigma2_mle**2)

    # Full per-sample Hessian (block-diagonal)
    J = np.zeros((k, k))
    J[:k_beta, :k_beta] = H_bb
    J[k_beta, k_beta] = H_ss

    log_det_J = np.linalg.slogdet(J)[1]

    log_n = np.log(n)

    # Scores (all in log-evidence units, higher = better)
    # AIC = -2 log L + 2k  →  maximisation score = -AIC/2 = log L - k = -nll - k
    aic_score = -nll - k                    # log L - k
    bic_score = -nll - (k / 2.0) * log_n
    exact_score = -nll - (k / 2.0) * log_n - 0.5 * log_det_J + (k / 2.0) * np.log(2 * np.pi)

    return ModelScore(
        name=name, k=k, nll=nll, S_vis=S_vis,
        log_det_J=log_det_J,
        aic_score=aic_score,
        bic_score=bic_score,
        exact_score=exact_score,
    )


def run_benchmark(n: int = 60, seed: int = 1, rho: float = 0.98):
    """Run the full benchmark for one seed. Returns (M1_score, M2_score)."""
    y, x1, x2 = generate_data(n=n, seed=seed, rho=rho)

    # M1: intercept + x1
    X1 = sm.add_constant(np.column_stack([x1]))
    m1 = fit_and_score(y, X1, "M1", n)

    # M2: intercept + x1 + x2
    X2 = sm.add_constant(np.column_stack([x1, x2]))
    m2 = fit_and_score(y, X2, "M2", n)

    return m1, m2


def run_sweep(n_seeds: int = 200, n: int = 60, rho: float = 0.98):
    """Run the benchmark over many seeds. Returns list of (M1, M2) pairs."""
    results = []
    for seed in range(1, n_seeds + 1):
        m1, m2 = run_benchmark(n=n, seed=seed, rho=rho)
        results.append((m1, m2))
    return results


def print_single(m1: ModelScore, m2: ModelScore, n: int = 60):
    """Print the deterministic flagship case."""
    print("=" * 72)
    print("Strike 1B: Regular-Gaussian Flagship Benchmark")
    print("=" * 72)
    print(f"  n = {n}, generating model = M2 (y = 1 + 2x₁ + 0.1x₂ + ε)")
    print()

    for m in (m1, m2):
        print(f"  {m.name}:  k={m.k}  NLL={m.nll:.4f}  S_vis={m.S_vis:.6f}  "
              f"log det(J)={m.log_det_J:.4f}")

    print()
    print(f"  {'Score':<20s}  {'M1':>12s}  {'M2':>12s}  {'Δ(M2-M1)':>12s}  {'Prefers':>8s}")
    print(f"  {'-'*20}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*8}")

    for label, s1, s2 in [
        ("AIC-style", m1.aic_score, m2.aic_score),
        ("BIC-style", m1.bic_score, m2.bic_score),
        ("Exact local", m1.exact_score, m2.exact_score),
    ]:
        delta = s2 - s1
        pref = "M2" if delta > 0 else "M1"
        print(f"  {label:<20s}  {s1:>12.4f}  {s2:>12.4f}  {delta:>+12.4f}  {pref:>8s}")

    print()
    det_correction = (-0.5 * m2.log_det_J + 0.5 * m1.log_det_J
                      + (m2.k / 2.0) * np.log(2 * np.pi)
                      - (m1.k / 2.0) * np.log(2 * np.pi))
    print(f"  Determinant correction (Exact - BIC shift): {det_correction:+.4f}")
    print()


def print_sweep(results, n: int = 60):
    """Print the sweep summary."""
    n_seeds = len(results)
    aic_prefers_m2 = 0
    bic_prefers_m2 = 0
    exact_prefers_m2 = 0
    aic_wrong_exact_right = 0
    bic_wrong_exact_right = 0
    det_corrections = []

    for m1, m2 in results:
        aic_m2 = m2.aic_score > m1.aic_score
        bic_m2 = m2.bic_score > m1.bic_score
        exact_m2 = m2.exact_score > m1.exact_score

        if aic_m2:
            aic_prefers_m2 += 1
        if bic_m2:
            bic_prefers_m2 += 1
        if exact_m2:
            exact_prefers_m2 += 1
        if not aic_m2 and exact_m2:
            aic_wrong_exact_right += 1
        if not bic_m2 and exact_m2:
            bic_wrong_exact_right += 1

        det_corr = (-0.5 * m2.log_det_J + 0.5 * m1.log_det_J
                    + (m2.k / 2.0) * np.log(2 * np.pi)
                    - (m1.k / 2.0) * np.log(2 * np.pi))
        det_corrections.append(det_corr)

    det_corrections = np.array(det_corrections)

    print("=" * 72)
    print(f"Strike 1B: Sweep over {n_seeds} seeds  (n={n})")
    print("=" * 72)
    print()
    print(f"  {'Criterion':<25s}  {'Prefers M2':>12s}  {'Rate':>8s}")
    print(f"  {'-'*25}  {'-'*12}  {'-'*8}")
    print(f"  {'AIC-style':<25s}  {aic_prefers_m2:>12d}  {aic_prefers_m2/n_seeds:>8.1%}")
    print(f"  {'BIC-style':<25s}  {bic_prefers_m2:>12d}  {bic_prefers_m2/n_seeds:>8.1%}")
    print(f"  {'Exact local evidence':<25s}  {exact_prefers_m2:>12d}  {exact_prefers_m2/n_seeds:>8.1%}")
    print()
    print(f"  Reversals (criterion picks M1, exact picks M2):")
    print(f"    AIC → Exact reversal:  {aic_wrong_exact_right:>4d} / {n_seeds}  ({aic_wrong_exact_right/n_seeds:.1%})")
    print(f"    BIC → Exact reversal:  {bic_wrong_exact_right:>4d} / {n_seeds}  ({bic_wrong_exact_right/n_seeds:.1%})")
    print()
    print(f"  Determinant correction (Exact - BIC shift for M2 vs M1):")
    print(f"    mean   = {det_corrections.mean():+.4f}")
    print(f"    std    = {det_corrections.std():.4f}")
    print(f"    min    = {det_corrections.min():+.4f}")
    print(f"    max    = {det_corrections.max():+.4f}")
    print(f"    median = {np.median(det_corrections):+.4f}")
    print()

    return {
        "n_seeds": n_seeds,
        "aic_prefers_m2": aic_prefers_m2,
        "bic_prefers_m2": bic_prefers_m2,
        "exact_prefers_m2": exact_prefers_m2,
        "aic_wrong_exact_right": aic_wrong_exact_right,
        "bic_wrong_exact_right": bic_wrong_exact_right,
        "det_corrections": det_corrections,
    }


if __name__ == "__main__":
    # A. Deterministic flagship case
    m1, m2 = run_benchmark()
    print_single(m1, m2)

    # B. Sweep
    results = run_sweep(n_seeds=200)
    print_sweep(results)
