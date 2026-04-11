import numpy as np
import pytest

from nomogeo import (
    affine_hidden_branch_reversal,
    guarded_fibre_dominance,
    staged_affine_hidden_elimination,
    variable_precision_affine_hidden_reduction,
)
from nomogeo.exceptions import InputValidationError


def _normalize_log_integral(logw: np.ndarray, x: np.ndarray) -> float:
    shift = float(np.max(logw))
    z = float(np.trapezoid(np.exp(logw - shift), x))
    return shift + float(np.log(z))


def test_variable_precision_affine_hidden_reduction_matches_numeric_integral() -> None:
    v = np.linspace(-0.7, 0.7, 81)
    h = np.linspace(-9.0, 9.0, 12001)
    A = 0.08 * v**4 + 0.03 * v
    J = (0.2 + 0.45 * np.sin(1.3 * v) + 0.1 * v**2)[:, None]
    D_scalar = 1.4 * np.exp(0.8 * v - 0.2 * v**2)
    D = D_scalar[:, None, None]

    result = variable_precision_affine_hidden_reduction(A, J, D)

    numeric = []
    for idx, v_i in enumerate(v):
        action = A[idx] + 0.5 * D_scalar[idx] * h**2 + J[idx, 0] * h
        numeric.append(-_normalize_log_integral(-action, h))
    numeric = np.asarray(numeric)
    centered_numeric = numeric - np.mean(numeric)
    centered_formula = result.visible_action - np.mean(result.visible_action)

    assert np.max(np.abs(centered_numeric - centered_formula)) <= 1e-12
    assert np.allclose(result.variational_action, A - 0.5 * J[:, 0] ** 2 / D_scalar)
    assert np.allclose(result.fibre_volume, 0.5 * np.log(D_scalar))
    assert "not a generic non-Gaussian" in result.metadata.notes[-1]


def test_variable_precision_affine_hidden_reduction_detects_logdet_branch_flip() -> None:
    A = np.array([0.0, 0.0])
    J = np.zeros((2, 1))
    D = np.array([[[0.1]], [[10.0]]])

    result = variable_precision_affine_hidden_reduction(A, J, D)

    assert np.allclose(result.variational_action, [0.0, 0.0])
    assert result.visible_action[0] < result.visible_action[1]
    assert result.visible_action[0] - result.visible_action[1] == pytest.approx(0.5 * np.log(0.1 / 10.0))


def test_affine_hidden_branch_reversal_detects_fibre_volume_flip() -> None:
    variational = np.array([0.0, 0.2])
    fibre = np.array([0.5, -0.1])

    result = affine_hidden_branch_reversal(variational, fibre)

    assert result.variational_winners == (0,)
    assert result.visible_winners == (1,)
    assert result.reversal
    assert not result.preserved
    assert result.branch_reversal_matrix[0, 1] == 1.0
    assert "not a generic" in result.metadata.notes[-1]


def test_guarded_fibre_dominance_refuses_small_denominator_ratio() -> None:
    fibre = np.array([-1.0, 1.0])
    flat_variational = np.array([2.0, 2.0])

    undefined = guarded_fibre_dominance(flat_variational, fibre, denominator_floor=1e-6)

    assert undefined.fibre_centered_norm > 0.0
    assert undefined.variational_centered_norm == pytest.approx(0.0)
    assert not undefined.ratio_defined
    assert undefined.ratio is None

    defined = guarded_fibre_dominance(np.array([0.0, 2.0]), fibre, denominator_floor=1e-6)
    assert defined.ratio_defined
    assert defined.ratio == pytest.approx(1.0)


def test_staged_affine_hidden_elimination_matches_one_step_all_orders() -> None:
    rng = np.random.default_rng(3102)
    raw = rng.normal(size=(5, 5))
    D = raw.T @ raw + 1.3 * np.eye(5)
    J = rng.normal(size=5)
    A = 0.4
    one_step = variable_precision_affine_hidden_reduction(A, J, D).visible_action.item()

    for order in ([0, 1, 2, 3, 4], [4, 1, 3, 0, 2], [2, 2]):
        if len(set(order)) != len(order):
            with pytest.raises(InputValidationError):
                staged_affine_hidden_elimination(A, J, D, order)
            continue
        labels = list(range(5))
        current_A = A
        current_J = J.copy()
        current_D = D.copy()
        for label in order:
            idx = labels.index(label)
            stage = staged_affine_hidden_elimination(current_A, current_J, current_D, [idx])
            current_A = stage.action
            current_J = stage.coupling
            current_D = stage.hidden_precision
            labels = [item for item in labels if item != label]
        assert stage.visible_action == pytest.approx(one_step)
        assert current_J.shape == (0,)
        assert current_D.shape == (0, 0)


def test_affine_hidden_reduction_rejects_bad_inputs() -> None:
    with pytest.raises(InputValidationError):
        variable_precision_affine_hidden_reduction(np.zeros(2), np.zeros((3, 1)), np.ones((2, 1, 1)))
    with pytest.raises(InputValidationError):
        variable_precision_affine_hidden_reduction(np.zeros(2), np.zeros((2, 1)), -np.ones((2, 1, 1)))
    with pytest.raises(InputValidationError):
        staged_affine_hidden_elimination(0.0, np.ones((2, 1)), np.eye(2), [0])
    with pytest.raises(InputValidationError):
        staged_affine_hidden_elimination(0.0, np.ones(2), np.eye(2), [])
    with pytest.raises(InputValidationError):
        staged_affine_hidden_elimination(0.0, np.ones(2), np.eye(2), [2])
    with pytest.raises(InputValidationError):
        affine_hidden_branch_reversal(np.zeros(2), np.zeros(3))
    with pytest.raises(InputValidationError):
        guarded_fibre_dominance(np.zeros(2), np.zeros(2), norm="bad")
