from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy import integrate, optimize, stats
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "audit" / "outputs"
OUT_JSON = OUT_DIR / "non_gaussian_stress_tests.json"
OUT_CSV = OUT_DIR / "non_gaussian_stress_tests.csv"


LogPDF = Callable[[float], float]


@dataclass(frozen=True)
class CaseResult:
    case: str
    quadratic_verdict: str
    full_law_verdict: str
    status: str
    metrics: dict[str, float | str]


def quad_exp(log_integrand: Callable[[float], float], epsabs: float = 1e-10) -> float:
    peak = _rough_peak(log_integrand)

    def shifted(x: float) -> float:
        return math.exp(log_integrand(x) - peak)

    value, _err = integrate.quad(shifted, -np.inf, np.inf, epsabs=epsabs, epsrel=1e-10, limit=300)
    return float(math.exp(peak) * value)


def _rough_peak(logf: Callable[[float], float]) -> float:
    xs = np.concatenate(
        [
            np.linspace(-12, 12, 481),
            np.array([-40.0, -25.0, 25.0, 40.0]),
        ]
    )
    return float(max(logf(float(x)) for x in xs))


def normalizer(log_unnormalized: LogPDF) -> float:
    return quad_exp(log_unnormalized)


def normalized_logpdf(log_unnormalized: LogPDF) -> LogPDF:
    log_z = math.log(normalizer(log_unnormalized))
    return lambda x: log_unnormalized(x) - log_z


def sym_kl(logp: LogPDF, logq: LogPDF) -> float:
    def p_int(x: float) -> float:
        lp = logp(x)
        return math.exp(lp) * (lp - logq(x))

    def q_int(x: float) -> float:
        lq = logq(x)
        return math.exp(lq) * (lq - logp(x))

    return float(_quad_real(p_int) + _quad_real(q_int))


def hellinger_sq(logp: LogPDF, logq: LogPDF) -> float:
    def affinity_int(x: float) -> float:
        return math.exp(0.5 * (logp(x) + logq(x)))

    affinity = _quad_real(affinity_int)
    return float(max(0.0, min(1.0, 1.0 - affinity)))


def _quad_real(f: Callable[[float], float]) -> float:
    value, _err = integrate.quad(f, -np.inf, np.inf, epsabs=1e-9, epsrel=1e-9, limit=300)
    return float(value)


def gaussian_logpdf(precision: float = 1.0) -> LogPDF:
    return lambda x: 0.5 * math.log(precision / (2.0 * math.pi)) - 0.5 * precision * x * x


def student_t_same_hessian_logpdf(nu: float, hessian: float = 1.0) -> LogPDF:
    scale = math.sqrt((nu + 1.0) / (nu * hessian))
    dist = stats.t(df=nu, loc=0.0, scale=scale)
    return lambda x: float(dist.logpdf(x))


def quartic_logpdf(hessian: float = 1.0, quartic: float = 0.0) -> LogPDF:
    return normalized_logpdf(lambda x: -0.5 * hessian * x * x - quartic * x**4)


def symmetric_tail_mixture_logpdf(eps: float, m: float, tail_sd: float, target_hessian: float = 1.0) -> tuple[LogPDF, float]:
    def hessian_for(s0: float) -> float:
        # For a symmetric mixture, f'(0)=0 and -d^2 log f(0) = -f''(0)/f(0).
        f0 = (1.0 - eps) * stats.norm.pdf(0.0, 0.0, s0) + eps * stats.norm.pdf(0.0, m, tail_sd)
        f2_central = (1.0 - eps) * stats.norm.pdf(0.0, 0.0, s0) * (-1.0 / s0**2)
        tail_pdf0 = stats.norm.pdf(0.0, m, tail_sd)
        f2_tail_one = tail_pdf0 * ((m * m / tail_sd**4) - (1.0 / tail_sd**2))
        f2 = f2_central + eps * f2_tail_one
        return float(-f2 / f0)

    root = optimize.brentq(lambda s: hessian_for(s) - target_hessian, 0.5, 2.0)

    def logpdf(x: float) -> float:
        logs = [
            math.log1p(-eps) + stats.norm.logpdf(x, 0.0, root),
            math.log(eps / 2.0) + stats.norm.logpdf(x, m, tail_sd),
            math.log(eps / 2.0) + stats.norm.logpdf(x, -m, tail_sd),
        ]
        return float(logsumexp(logs))

    return logpdf, float(root)


def truncated_left_gaussian_logpdf(bound: float) -> LogPDF:
    # Standard normal conditioned on x >= bound. Around x=0 with bound < 0,
    # the log-density Hessian is exactly -1, while the support is different.
    log_z = math.log(1.0 - stats.norm.cdf(bound))

    def logpdf(x: float) -> float:
        if x < bound:
            return -math.inf
        return float(stats.norm.logpdf(x) - log_z)

    return logpdf


def closure_cubic_marginal_logpdf(t: float, alpha: float) -> LogPDF:
    def inner_log_integral(x: float) -> float:
        # Integral of exp(-alpha z^4 - b z^2), with b allowed to be negative.
        b = 0.5 + t * x
        peak = 0.0 if b >= 0.0 else (b * b) / (4.0 * alpha)

        def shifted(z: float) -> float:
            return math.exp((-alpha * z**4 - b * z * z) - peak)

        value, _err = integrate.quad(shifted, -np.inf, np.inf, epsabs=1e-10, epsrel=1e-10, limit=300)
        if value <= 0.0:
            return -1.0e300
        return float(peak + math.log(value))

    log_unnorm = lambda x: -0.5 * x * x - alpha * x**4 + inner_log_integral(x)
    return normalized_logpdf(log_unnorm)


def case_false_collapse_tail() -> CaseResult:
    p = gaussian_logpdf(1.0)
    q = student_t_same_hessian_logpdf(nu=5.0, hessian=1.0)
    return CaseResult(
        case="same local Hessian, heavy-tail law",
        quadratic_verdict="collapse: H=1 for both laws, so Phi and all quadratic objects agree",
        full_law_verdict="visible 1D laws differ strongly in tails",
        status="pathology",
        metrics={
            "symmetric_KL": sym_kl(p, q),
            "hellinger_sq": hellinger_sq(p, q),
            "student_df": 5.0,
            "student_scale_for_H=1": math.sqrt(6.0 / 5.0),
        },
    )


def case_observer_ranking_reversal() -> CaseResult:
    x0 = gaussian_logpdf(1.0)
    x1 = gaussian_logpdf(1.08)
    y0 = quartic_logpdf(1.0, 0.0)
    y1 = quartic_logpdf(1.0, 0.35)
    skl_x = sym_kl(x0, x1)
    skl_y = sym_kl(y0, y1)
    return CaseResult(
        case="quadratic ranks x, full law ranks y",
        quadratic_verdict="x observer wins: Delta_H=(0.08, 0), y has matched Hessian",
        full_law_verdict="y observer wins by higher-order tail change",
        status="pathology",
        metrics={
            "sym_KL_x_precision_shift": skl_x,
            "sym_KL_y_quartic_tail_shift": skl_y,
            "ratio_y_over_x": skl_y / skl_x,
            "quadratic_delta_x": 0.08,
            "quadratic_delta_y": 0.0,
        },
    )


def case_quadratic_closure_higher_order_leakage() -> CaseResult:
    p0 = closure_cubic_marginal_logpdf(t=0.0, alpha=0.15)
    p1 = closure_cubic_marginal_logpdf(t=0.28, alpha=0.15)

    def mean(logp: LogPDF) -> float:
        return _quad_real(lambda x: x * math.exp(logp(x)))

    def second(logp: LogPDF) -> float:
        mu = mean(logp)
        return _quad_real(lambda x: (x - mu) ** 2 * math.exp(logp(x)))

    return CaseResult(
        case="cubic hidden-visible coupling invisible to Hessian",
        quadratic_verdict="no visible response and no closure leakage at order two: Hessian at the mode is unchanged",
        full_law_verdict="visible marginal of x changes through integration over hidden z",
        status="pathology",
        metrics={
            "sym_KL_visible_marginal_x": sym_kl(p0, p1),
            "hellinger_sq_visible_marginal_x": hellinger_sq(p0, p1),
            "mean_x_t0": mean(p0),
            "mean_x_t028": mean(p1),
            "var_x_t0": second(p0),
            "var_x_t028": second(p1),
        },
    )


def case_hidden_branch_same_hessian() -> CaseResult:
    base = gaussian_logpdf(1.0)
    mix, central_sd = symmetric_tail_mixture_logpdf(eps=0.05, m=5.0, tail_sd=0.35, target_hessian=1.0)
    branch_base = 2.0 * (1.0 - stats.norm.cdf(3.0))
    branch_mix = _quad_real(lambda x: (1.0 if abs(x) > 3.0 else 0.0) * math.exp(mix(x)))
    return CaseResult(
        case="matched Hessian with remote symmetric branch",
        quadratic_verdict="collapse near the central mode: H=1 after central width calibration",
        full_law_verdict="remote branch is visible to threshold/event observers",
        status="pathology",
        metrics={
            "sym_KL": sym_kl(base, mix),
            "branch_prob_abs_gt_3_base": branch_base,
            "branch_prob_abs_gt_3_mixture": branch_mix,
            "calibrated_central_sd": central_sd,
        },
    )


def case_support_boundary_same_hessian() -> CaseResult:
    bound = -2.0
    # KL(base || trunc) is infinite because base assigns mass to x < bound.
    # KL(trunc || base) is analytic for truncation of a normalized base law.
    kl_trunc_base = -math.log(1.0 - stats.norm.cdf(bound))
    return CaseResult(
        case="same interior Hessian, different support",
        quadratic_verdict="same local quadratic precision H=1 at x=0",
        full_law_verdict="support verdict differs: truncated law forbids an event with positive Gaussian probability",
        status="pathology",
        metrics={
            "KL_truncated_to_gaussian": kl_trunc_base,
            "KL_gaussian_to_truncated": "infinite",
            "gaussian_prob_x_lt_minus_2": stats.norm.cdf(bound),
            "truncated_prob_x_lt_minus_2": 0.0,
        },
    )


def case_student_t_scale_survival() -> CaseResult:
    base = student_t_same_hessian_logpdf(nu=7.0, hessian=1.0)
    strong = student_t_same_hessian_logpdf(nu=7.0, hessian=1.30)
    weak = student_t_same_hessian_logpdf(nu=7.0, hessian=1.08)
    skl_strong = sym_kl(base, strong)
    skl_weak = sym_kl(base, weak)
    return CaseResult(
        case="elliptic heavy-tail scale perturbation preserves ranking",
        quadratic_verdict="strong precision perturbation should outrank weak perturbation",
        full_law_verdict="Student-t visible KL preserves the same ordering",
        status="survival",
        metrics={
            "sym_KL_strong": skl_strong,
            "sym_KL_weak": skl_weak,
            "ratio_strong_over_weak": skl_strong / skl_weak,
            "nu": 7.0,
        },
    )


def case_small_quartic_survival() -> CaseResult:
    base = quartic_logpdf(1.0, 0.03)
    strong = quartic_logpdf(1.30, 0.03)
    weak = quartic_logpdf(1.08, 0.03)
    skl_strong = sym_kl(base, strong)
    skl_weak = sym_kl(base, weak)
    return CaseResult(
        case="small non-Gaussian quartic family preserves local ordering",
        quadratic_verdict="strong Hessian perturbation should outrank weak perturbation",
        full_law_verdict="full non-Gaussian KL preserves the same qualitative observer ranking",
        status="survival",
        metrics={
            "sym_KL_strong": skl_strong,
            "sym_KL_weak": skl_weak,
            "ratio_strong_over_weak": skl_strong / skl_weak,
            "quartic": 0.03,
        },
    )


def run() -> list[CaseResult]:
    return [
        case_false_collapse_tail(),
        case_observer_ranking_reversal(),
        case_quadratic_closure_higher_order_leakage(),
        case_hidden_branch_same_hessian(),
        case_support_boundary_same_hessian(),
        case_student_t_scale_survival(),
        case_small_quartic_survival(),
    ]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = run()
    serial = [r.__dict__ for r in results]
    OUT_JSON.write_text(json.dumps(serial, indent=2, sort_keys=True), encoding="utf-8")

    metric_keys = sorted({key for result in results for key in result.metrics})
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["case", "status", "quadratic_verdict", "full_law_verdict", *metric_keys],
        )
        writer.writeheader()
        for result in results:
            row = {
                "case": result.case,
                "status": result.status,
                "quadratic_verdict": result.quadratic_verdict,
                "full_law_verdict": result.full_law_verdict,
            }
            row.update(result.metrics)
            writer.writerow(row)
    print(OUT_JSON)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
