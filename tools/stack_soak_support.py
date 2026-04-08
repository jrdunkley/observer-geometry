from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for workspace_src in (ROOT / "nomodescent" / "src", ROOT / "evidence" / "src"):
    path_str = str(workspace_src)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from nomogeo import (
    clock,
    hidden_contraction,
    hidden_load,
    load_from_hidden_contraction,
    transport_hidden_load,
    visible_from_hidden_load,
    visible_geometry,
    visible_precision,
)
from nomodescent import (
    AssumptionEntry,
    AssumptionLedger,
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
    classify_relation,
    common_descent_test,
    minimal_refinement_search,
    staged_descent_check,
)
from evidence import (
    EvidenceBundle,
    ExtractionNote,
    ExtractionRecord,
    ObserverHypothesis,
    ProtocolObservation,
    SourceRef,
    assemble_problem_spec,
    encode_matrix_observation,
    infer_observer_candidates,
)


def execute_soak_task(task_kind: str, seed: int, variant: str = "default") -> dict[str, object]:
    if task_kind == "kernel":
        return _kernel_task(seed, variant)
    if task_kind == "descent":
        return _descent_task(seed, variant)
    if task_kind == "evidence":
        return _evidence_task(seed, variant)
    if task_kind == "example":
        return _example_task(seed, variant)
    raise ValueError(f"unknown task_kind '{task_kind}'")


def fail_soak_task(label: str, fail: bool = False) -> str:
    if fail:
        raise RuntimeError(f"intentional failure for {label}")
    return label


def _kernel_task(seed: int, variant: str) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    n = int(rng.integers(3, 7))
    m = int(rng.integers(1, n + 1))
    if variant == "boundary":
        h_eigs = np.geomspace(1e-6, 1.0, n)
        H = np.diag(h_eigs)
        C = _random_surjective(rng, m, n)
        phi = visible_precision(H, C)
        support_rank = max(1, min(n - 1, int(rng.integers(1, n + 1))))
        q = _orthogonal(rng, n)
        basis = q[:, :support_rank]
        values = np.geomspace(1e-4, 1.0, support_rank)
        T = basis @ np.diag(values) @ basis.T
        lambda_reduced = np.diag(np.geomspace(1e-8, 1e3, support_rank))
        X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
        load = hidden_load(T, X)
        zero_visible = visible_from_hidden_load(np.zeros((n, n), dtype=float), np.zeros((0, 0), dtype=float), lambda_representation="reduced")
        return {
            "task_kind": "kernel",
            "variant": variant,
            "seed": seed,
            "phi_condition": float(np.linalg.cond(phi)),
            "inverse_roundtrip_residual": float(np.linalg.norm(load.reduced_lambda - lambda_reduced, ord=np.inf)),
            "zero_support_residual": float(np.linalg.norm(zero_visible, ord=np.inf)),
            "large_load_clock": float(clock(lambda_reduced)),
        }

    H = _random_spd(rng, n)
    C = _random_surjective(rng, m, n)
    geometry = visible_geometry(H, C)
    support_rank = int(rng.integers(0, n + 1))
    q = _orthogonal(rng, n)
    basis = q[:, :support_rank]
    T = basis @ np.diag(np.linspace(1.0, 2.0, support_rank)) @ basis.T if support_rank else np.zeros((n, n), dtype=float)
    lambda_reduced = _random_psd(rng, support_rank, scale=0.05) if support_rank else np.zeros((0, 0), dtype=float)
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    load = hidden_load(T, X)

    chain_length = 24
    factor_total = np.eye(support_rank, dtype=float) if support_rank else np.zeros((0, 0), dtype=float)
    chain_clock = 0.0
    for _ in range(chain_length):
        load_i = _random_psd(rng, support_rank, scale=0.03) if support_rank else np.zeros((0, 0), dtype=float)
        if support_rank:
            factor_total = hidden_contraction(load_i) @ factor_total
            chain_clock += clock(load_i)
    composed = load_from_hidden_contraction(factor_total) if support_rank else np.zeros((0, 0), dtype=float)
    long_chain_clock = clock(composed) if support_rank else 0.0

    other = _random_psd(rng, support_rank, scale=0.05) if support_rank else np.zeros((0, 0), dtype=float)
    transported = transport_hidden_load(lambda_reduced, other) if support_rank else np.zeros((0, 0), dtype=float)
    transported_direct = load_from_hidden_contraction(hidden_contraction(other) @ hidden_contraction(lambda_reduced)) if support_rank else np.zeros((0, 0), dtype=float)
    projector_residual = max(
        float(np.linalg.norm(C @ geometry.lift - np.eye(m), ord=np.inf)),
        float(np.linalg.norm(C @ geometry.projector, ord=np.inf)),
        float(np.linalg.norm(geometry.projector @ geometry.projector - geometry.projector, ord=np.inf)),
        float(np.linalg.norm(geometry.projector.T @ H - H @ geometry.projector, ord=np.inf)),
    )
    return {
        "task_kind": "kernel",
        "variant": variant,
        "seed": seed,
        "projector_residual": projector_residual,
        "inverse_roundtrip_residual": float(np.linalg.norm(load.reduced_lambda - lambda_reduced, ord=np.inf)),
        "two_step_transport_residual": float(np.linalg.norm(transported - transported_direct, ord=np.inf)),
        "long_chain_clock_residual": abs(float(long_chain_clock - chain_clock)),
        "hidden_rank_gap": int(abs(load.rank - np.linalg.matrix_rank(T - X, tol=max(load.metadata.rank_tol, 1e-12)))),
    }


def _descent_task(seed: int, variant: str) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    if variant == "factorisation":
        fine = ObserverSpec(name="fine", matrix=_random_surjective(rng, 3, 5))
        D = rng.normal(size=(2, 3))
        coarse = ObserverSpec(name="coarse", matrix=D @ fine.matrix)
        result = classify_relation(coarse, fine)
        return {
            "task_kind": "descent",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "a_residual": float(result.residuals["a_through_b_residual"]),
            "b_residual": float(result.residuals["b_through_a_residual"]),
        }
    if variant == "tower":
        H = _random_spd(rng, 5)
        first = ObserverSpec(name="first", matrix=_random_surjective(rng, 4, 5))
        second = ObserverSpec(name="second", matrix=_random_surjective(rng, 2, 4))
        result = staged_descent_check(H, [first, second])
        return {
            "task_kind": "descent",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "tower_residual": float(result.residuals["tower_residual"]),
        }
    if variant == "completion_exact":
        sigma = _random_spd(rng, 3)
        problem = ProblemSpec(
            name=f"completion_exact_{seed}",
            latent_dim=3,
            observers=(
                ObserverSpec(name="full", matrix=np.eye(3)),
                ObserverSpec(name="panel", matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            ),
            evidence=(
                VisibleEvidenceSpec(name="sigma_full", observer="full", kind="covariance", matrix=sigma),
                VisibleEvidenceSpec(name="sigma_panel", observer="panel", kind="covariance", matrix=sigma[:2, :2]),
            ),
            assumptions=AssumptionLedger(entries=(AssumptionEntry(label="gaussian", statement="Gaussian common descent", exact=True),)),
            goals=(GoalSpec(kind="common_completion"),),
        )
        result = common_descent_test(problem)
        return {
            "task_kind": "descent",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "linear_residual": float(result.residuals["linear_residual"]),
            "psd_margin": float(result.residuals["psd_margin"]),
        }
    if variant == "completion_incompatible":
        problem = ProblemSpec(
            name=f"completion_bad_{seed}",
            latent_dim=2,
            observers=(ObserverSpec(name="obs", matrix=[[1.0, 0.0]]),),
            evidence=(
                VisibleEvidenceSpec(name="obs_1", observer="obs", kind="covariance", matrix=[[1.0]]),
                VisibleEvidenceSpec(name="obs_2", observer="obs", kind="covariance", matrix=[[2.0]]),
            ),
            goals=(GoalSpec(kind="common_completion"),),
        )
        result = common_descent_test(problem)
        return {
            "task_kind": "descent",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "linear_residual": float(result.residuals["linear_residual"]),
            "certificate_count": len(result.certificates),
        }
    if variant == "refinement":
        H = _random_spd(rng, 4)
        target_a = ObserverSpec(name="target_a", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        target_b = ObserverSpec(name="target_b", matrix=[[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        candidates = [
            ObserverSpec(name="panel_3", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
            ObserverSpec(name="full_4", matrix=np.eye(4)),
        ]
        result = minimal_refinement_search(H, [target_a, target_b], candidates)
        return {
            "task_kind": "descent",
            "variant": variant,
            "seed": seed,
            "winner": result.winner,
            "score_panel_3": float(result.scores["panel_3"]),
            "score_full_4": float(result.scores["full_4"]),
        }
    raise ValueError(f"unknown descent variant '{variant}'")


def _evidence_task(seed: int, variant: str) -> dict[str, object]:
    if variant == "underdetermined":
        bundle = _ambiguous_bundle(seed)
        result = assemble_problem_spec(bundle)
        return {
            "task_kind": "evidence",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "required_decisions": len(result.required_human_decisions),
            "unresolved_count": len(result.audit.unresolved_items),
        }
    if variant == "assembled":
        bundle = _ambiguous_bundle(seed)
        suggestion = infer_observer_candidates(bundle, "protocol")
        result = assemble_problem_spec(bundle, observer_selection={"protocol": suggestion.ranked_candidates[0]})
        return {
            "task_kind": "evidence",
            "variant": variant,
            "seed": seed,
            "classification": result.classification,
            "selected_observer": suggestion.ranked_candidates[0],
            "exact_flag": result.exact,
        }
    if variant == "suggestion":
        bundle = _ambiguous_bundle(seed)
        result = infer_observer_candidates(bundle, "protocol")
        return {
            "task_kind": "evidence",
            "variant": variant,
            "seed": seed,
            "top_candidate": result.ranked_candidates[0],
            "score_gap": float(result.audit.residuals["score_gap_top2"]),
            "candidate_count": len(result.ranked_candidates),
        }
    raise ValueError(f"unknown evidence variant '{variant}'")


def _example_task(seed: int, variant: str) -> dict[str, object]:
    if variant == "entanglement":
        from examples.entanglement_hidden_load.run_main import entanglement_hidden_load, two_mode_squeezed_thermal_covariance

        r = 0.2 + 0.01 * (seed % 50)
        summary = entanglement_hidden_load(two_mode_squeezed_thermal_covariance(r))
        return {
            "task_kind": "example",
            "variant": variant,
            "seed": seed,
            "tau_residual": abs(float(summary["tau"] - 2.0 * summary["mutual_information"])),
            "tau": float(summary["tau"]),
        }
    if variant == "bell":
        from examples.bell_common_gluing.run_main import classify_family, linear_completion_residual

        sample = ("compatible", 0.0, 0.6, True) if seed % 3 == 0 else ("variance_only", 0.22, 0.6, False) if seed % 3 == 1 else ("correlator_only", 0.0, 0.82, True)
        name, delta, rho, allow_approx = sample
        family = {
            "00": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
            "01": np.array([[1.0 + delta, rho * np.sqrt(1.0 + delta)], [rho * np.sqrt(1.0 + delta), 1.0]], dtype=float),
            "10": np.array([[1.0, rho], [rho, 1.0]], dtype=float),
            "11": np.array([[1.0, -rho], [-rho, 1.0]], dtype=float),
        }
        bundle = EvidenceBundle(
            name=f"bell_example_{name}_{seed}",
            latent_dim=4,
            protocol_observations=tuple(
                ProtocolObservation(
                    name=key,
                    facts=(f"A{key[0]}", f"B{key[1]}"),
                    candidate_family=(key,),
                    source_ref=SourceRef(source="bell_example", location=f"context {key}"),
                    extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
                )
                for key in ("00", "01", "10", "11")
            ),
            observer_hypotheses=tuple(
                ObserverHypothesis(
                    name=key,
                    matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]] if key == "00" else [[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]] if key == "01" else [[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]] if key == "10" else [[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
                    protocol_name=key,
                    features=(f"A{key[0]}", f"B{key[1]}"),
                    source_ref=SourceRef(source="bell_example", location="observer family"),
                    extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
                )
                for key in ("00", "01", "10", "11")
            ),
            matrix_observations=tuple(
                encode_matrix_observation(
                    name=f"sigma_{key}",
                    matrix_role="visible_object",
                    matrix_kind="covariance",
                    matrix=tuple(tuple(float(entry) for entry in row) for row in matrix),
                    observer_name=key,
                    source="bell_example",
                    location=f"context {key}",
                    quote="Pairwise covariance used in example Bell sweep.",
                )
                for key, matrix in family.items()
            ),
        )
        assembly = assemble_problem_spec(bundle)
        if assembly.problem_spec is None:
            raise RuntimeError("Bell example soak bundle should assemble")
        result = common_descent_test(assembly.problem_spec, allow_approximate_psd_search=allow_approx, psd_search_grid_size=31)
        return {
            "task_kind": "example",
            "variant": variant,
            "seed": seed,
            "sample": name,
            "phase_class": int(classify_family(delta, rho)),
            "downstream_classification": result.classification,
            "linear_residual": float(
                linear_completion_residual(
                    {
                        (0, 0): family["00"],
                        (0, 1): family["01"],
                        (1, 0): family["10"],
                        (1, 1): family["11"],
                    }
                )
            ),
        }
    if variant == "arrow":
        from examples.arrow_rank_deficiency.run_main import latent_precision, observer_ceiling, observers

        H = latent_precision()
        clocks = []
        for name, observer in observers().items():
            phi = visible_precision(H, observer)
            load = hidden_load(observer_ceiling(name, H), phi, support_mode="ambient")
            clocks.append(float(load.clock))
        return {
            "task_kind": "example",
            "variant": variant,
            "seed": seed,
            "clock_order_ok": bool(clocks[0] >= clocks[1] >= clocks[2]),
            "clock_gap": float(clocks[0] - clocks[2]),
        }
    raise ValueError(f"unknown example variant '{variant}'")


def _ambiguous_bundle(seed: int) -> EvidenceBundle:
    rng = np.random.default_rng(seed)
    cov = np.array([[1.0 + 0.01 * (seed % 5)]], dtype=float)
    protocol = ProtocolObservation(
        name="protocol",
        facts=("measure_x1", "measure_x2"),
        candidate_family=("direct", "aggregate"),
        source_ref=SourceRef(source="soak_bundle", location="protocol section"),
        extraction=ExtractionRecord(extraction_mode="exact_extraction", epistemic_status="exact", authoritative=True),
    )
    hypotheses = (
        ObserverHypothesis(
            name="direct",
            matrix=[[1.0, 0.0]],
            protocol_name="protocol",
            features=("measure_x1",),
            source_ref=SourceRef(source="soak_bundle", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol encoded as direct x1 measurement",),
        ),
        ObserverHypothesis(
            name="aggregate",
            matrix=[[1.0, 1.0]],
            protocol_name="protocol",
            features=("measure_x1", "measure_x2"),
            source_ref=SourceRef(source="soak_bundle", location="candidate family"),
            extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=True),
            assumption_statements=("protocol encoded as x1+x2 aggregate measurement",),
        ),
    )
    matrix = encode_matrix_observation(
        name="cov_protocol",
        matrix_role="visible_object",
        matrix_kind="covariance",
        matrix=((float(cov[0, 0]),),),
        observer_name="protocol",
        source="soak_bundle",
        location="result table",
        quote="Visible covariance under the unresolved protocol observer.",
    )
    note = ExtractionNote(
        name="model_choice",
        statement="observer matrix is inferred from protocol wording and remains unresolved without explicit selection",
        load_bearing=True,
        source_ref=SourceRef(source="soak_bundle"),
        extraction=ExtractionRecord(extraction_mode="encoded_inference", epistemic_status="inferred", authoritative=False),
    )
    return EvidenceBundle(
        name=f"ambiguous_bundle_{seed}",
        latent_dim=2,
        matrix_observations=(matrix,),
        protocol_observations=(protocol,),
        observer_hypotheses=hypotheses,
        notes=(note,),
        description=f"Soak bundle seeded with {seed}",
    )


def _random_spd(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return a.T @ a + n * np.eye(n, dtype=float)


def _random_surjective(rng: np.random.Generator, m: int, n: int) -> np.ndarray:
    while True:
        matrix = rng.normal(size=(m, n))
        if np.linalg.matrix_rank(matrix) == m:
            return matrix


def _random_psd(rng: np.random.Generator, n: int, scale: float = 1.0) -> np.ndarray:
    if n == 0:
        return np.zeros((0, 0), dtype=float)
    factor = rng.normal(size=(n, n))
    return scale * (factor @ factor.T)


def _orthogonal(rng: np.random.Generator, n: int) -> np.ndarray:
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return q

