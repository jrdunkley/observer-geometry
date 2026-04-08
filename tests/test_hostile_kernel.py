from __future__ import annotations

import numpy as np
import pytest

from nomogeo import (
    batch_map,
    clock,
    dv_bridge,
    gaussian_data_processing_contraction,
    hidden_load,
    local_visible_calculus,
    transport_hidden_load,
    visible_from_hidden_load,
    visible_precision,
)
from nomogeo.batch import ThreadPoolExecutor
from nomogeo.exceptions import BatchTaskError, InputValidationError, SupportError

from .helpers import add_pair, fail_on_negative, random_psd_on_subspace, random_spd, random_surjective


def test_domain_rejection_and_useful_exceptions() -> None:
    H = np.eye(2)
    non_surjective = np.array([[1.0, 0.0], [2.0, 0.0]])
    with pytest.raises(InputValidationError, match="surjective"):
        visible_precision(H, non_surjective)

    with pytest.raises(InputValidationError, match="positive definite"):
        visible_precision(np.array([[1.0, 0.0], [0.0, 0.0]]), np.array([[1.0, 0.0]]))

    with pytest.raises(InputValidationError, match="incompatible latent dimension"):
        local_visible_calculus(np.eye(2), np.array([[1.0, 0.0, 0.0]]), np.eye(2))

    T = np.array([[1.0, 0.0], [0.0, 0.0]])
    X_bad = np.array([[0.5, 0.2], [0.2, 0.1]])
    with pytest.raises(SupportError, match="preserve the support"):
        hidden_load(T, X_bad)

    with pytest.raises(InputValidationError, match="positive semidefinite"):
        clock(np.diag([-0.1, 0.2]))

    with pytest.raises(InputValidationError, match="positive semidefinite"):
        transport_hidden_load(np.diag([-0.1, 0.2]), np.diag([0.2, 0.3]))

    with pytest.raises(InputValidationError, match="ambiguous"):
        visible_from_hidden_load(np.array([[2.0, 0.7], [0.7, 1.5]]), np.diag([0.2, 0.5]))

    with pytest.raises(InputValidationError, match="support_mode"):
        hidden_load(np.eye(2), np.eye(2), support_mode="bad")

    with pytest.raises(InputValidationError, match="support_mode"):
        visible_from_hidden_load(np.eye(2), np.eye(2), support_mode="bad", lambda_representation="ambient")

    with pytest.raises(InputValidationError, match="backend"):
        batch_map(add_pair, [(1, 2)], backend="bad")

    with pytest.raises(BatchTaskError, match="batch task 0 failed"):
        batch_map(add_pair, [{"x": 1}], backend="serial")


def test_near_boundary_stability_and_scale_behaviour() -> None:
    H = np.diag([1e-4, 1.0, 3.0])
    C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    phi = visible_precision(H, C)
    assert np.all(np.isfinite(phi))
    assert np.allclose(phi, phi.T, atol=1e-12, rtol=1e-12)

    T = np.diag([1.0, 1e-5, 0.0])
    lambda_small = np.diag([1e-7, 5e-7])
    X = visible_from_hidden_load(T, lambda_small, lambda_representation="reduced")
    load = hidden_load(T, X)
    assert np.allclose(load.reduced_lambda, lambda_small, atol=1e-9, rtol=1e-3)

    H_base = random_spd(np.random.default_rng(901), 4)
    C_base = random_surjective(np.random.default_rng(902), 2, 4)
    phi_base = visible_precision(H_base, C_base)
    assert np.allclose(visible_precision(1e-6 * H_base, C_base), 1e-6 * phi_base, atol=1e-12, rtol=1e-8)
    assert np.allclose(visible_precision(1e6 * H_base, C_base), 1e6 * phi_base, atol=1e-3, rtol=1e-8)

    T_base, _basis = random_psd_on_subspace(np.random.default_rng(903), 5, 3)
    lambda_base = np.diag([0.2, 0.4, 0.8])
    X_base = visible_from_hidden_load(T_base, lambda_base, lambda_representation="reduced")
    scaled = hidden_load(1e-6 * T_base, 1e-6 * X_base).reduced_lambda
    assert np.allclose(scaled, lambda_base, atol=1e-9, rtol=1e-9)


def test_public_functions_do_not_mutate_inputs() -> None:
    rng = np.random.default_rng(904)
    H = random_spd(rng, 4)
    C = random_surjective(rng, 2, 4)
    raw_delta = rng.normal(size=(4, 4))
    Delta = 0.5 * (raw_delta + raw_delta.T)
    T, _basis = random_psd_on_subspace(rng, 4, 2)
    lambda_reduced = np.diag([0.2, 0.5])
    X = visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    Jhat = np.array([[0.0, 1.0, 0.0, -0.5], [-1.0, 0.0, 0.3, 0.0], [0.0, -0.3, 0.0, 0.2], [0.5, 0.0, -0.2, 0.0]])

    copies = [arr.copy() for arr in (H, C, Delta, T, lambda_reduced, X, Jhat)]
    visible_precision(H, C)
    local_visible_calculus(H, C, Delta)
    hidden_load(T, X)
    visible_from_hidden_load(T, lambda_reduced, lambda_representation="reduced")
    transport_hidden_load(lambda_reduced, lambda_reduced)
    clock(lambda_reduced)
    dv_bridge(H, Jhat)
    gaussian_data_processing_contraction(H, H + np.eye(4), np.eye(4), np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0]]))

    for original, copy in zip((H, C, Delta, T, lambda_reduced, X, Jhat), copies, strict=True):
        assert np.array_equal(original, copy)


def test_batch_determinism_failure_propagation_and_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    rng = np.random.default_rng(905)
    tasks = []
    for _ in range(6):
        H = random_spd(rng, 4)
        C = random_surjective(rng, 2, 4)
        tasks.append({"H": H, "C": C})

    from .helpers import visible_precision_trace_task

    serial = batch_map(visible_precision_trace_task, tasks, backend="serial")
    thread = batch_map(visible_precision_trace_task, tasks, backend="thread", max_workers=2)
    process = batch_map(visible_precision_trace_task, tasks, backend="process", max_workers=2)

    assert serial == thread == process

    failing_tasks = [1, -1, 2]
    for backend in ("serial", "thread", "process"):
        with pytest.raises(BatchTaskError) as excinfo:
            batch_map(fail_on_negative, failing_tasks, backend=backend, max_workers=2)
        assert excinfo.value.index == 1
        assert excinfo.value.task == -1
        assert "batch task 1 failed" in str(excinfo.value)

    class BrokenProcessPool:
        def __init__(self, *args: object, **kwargs: object) -> None:
            raise PermissionError("blocked process pool")

    monkeypatch.setattr("nomogeo.batch.ProcessPoolExecutor", BrokenProcessPool)
    fallback = batch_map(visible_precision_trace_task, tasks, backend="process", max_workers=2)
    assert fallback == serial
