"""Tests for tower_affine_hidden_elimination (affine.py 0.3.3 addition)."""
import numpy as np
import pytest

from nomogeo import (
    staged_affine_hidden_elimination,
    variable_precision_affine_hidden_reduction,
)
from nomogeo.affine import tower_affine_hidden_elimination
from nomogeo.regime_types import TowerEliminationRecord, TowerStageRecord
from nomogeo.exceptions import InputValidationError


def _random_spd(d: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(d, d))
    return raw.T @ raw + 0.5 * np.eye(d)


# ── Basic operation ──────────────────────────────────────────────────


def test_single_stage_matches_staged_elimination() -> None:
    """One tower stage should exactly match staged_affine_hidden_elimination."""
    D = _random_spd(4, seed=10)
    J = np.array([0.3, -0.1, 0.7, 0.5])
    A = 1.2

    tower = tower_affine_hidden_elimination(
        A, J, D, stages=[("all", [0, 1, 2, 3])]
    )
    staged = staged_affine_hidden_elimination(A, J, D, [0, 1, 2, 3])

    assert tower.final_action == pytest.approx(staged.visible_action)
    assert tower.final_coupling.shape == (0,)
    assert tower.final_precision.shape == (0, 0)
    assert len(tower.stages) == 1
    assert tower.stages[0].block_name == "all"
    assert tower.stages[0].hidden_dim == 4


def test_two_stage_matches_one_shot() -> None:
    """Eliminating in two stages should match eliminating all at once."""
    D = _random_spd(6, seed=20)
    J = np.array([0.1, -0.2, 0.3, -0.4, 0.5, -0.6])
    A = 0.5

    # One-shot via variable_precision
    one_shot = variable_precision_affine_hidden_reduction(A, J, D)
    one_shot_action = float(one_shot.visible_action)

    # Two-stage tower: eliminate indices [0,1,2] then [0,1,2] of remainder
    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("layer_1", [0, 1, 2]), ("layer_2", [0, 1, 2])],
    )

    assert tower.final_action == pytest.approx(one_shot_action)
    assert tower.final_coupling.shape == (0,)
    assert tower.final_precision.shape == (0, 0)
    assert len(tower.stages) == 2
    assert tower.stages[0].block_name == "layer_1"
    assert tower.stages[0].hidden_dim == 3
    assert tower.stages[1].block_name == "layer_2"
    assert tower.stages[1].hidden_dim == 3


def test_three_stage_matches_one_shot() -> None:
    """Three stages eliminating one variable at a time."""
    D = _random_spd(3, seed=30)
    J = np.array([1.0, 2.0, 3.0])
    A = 0.0

    one_shot = variable_precision_affine_hidden_reduction(A, J, D)
    one_shot_action = float(one_shot.visible_action)

    # Eliminate one variable at a time — always index 0 of remainder.
    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("var_0", [0]), ("var_1", [0]), ("var_2", [0])],
    )

    assert tower.final_action == pytest.approx(one_shot_action)
    assert len(tower.stages) == 3
    for s in tower.stages:
        assert s.hidden_dim == 1


def test_partial_elimination_leaves_residual() -> None:
    """Eliminating fewer than all variables leaves a residual."""
    D = _random_spd(5, seed=40)
    J = np.ones(5)
    A = 1.0

    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("first_two", [0, 1])],
    )

    assert tower.final_coupling.shape == (3,)
    assert tower.final_precision.shape == (3, 3)
    # Residual precision should be SPD (Schur complement of SPD is SPD).
    eigvals = np.linalg.eigvalsh(tower.final_precision)
    assert np.all(eigvals > 0)


# ── Log-determinant accumulation ─────────────────────────────────────


def test_total_half_log_det_equals_sum_of_increments() -> None:
    """total_half_log_det should equal the sum of per-stage half_log_det_increments."""
    D = _random_spd(6, seed=50)
    J = np.zeros(6)
    A = 0.0

    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("block_a", [0, 1]), ("block_b", [0, 1]), ("block_c", [0, 1])],
    )

    sum_inc = sum(s.half_log_det_increment for s in tower.stages)
    assert tower.total_half_log_det == pytest.approx(sum_inc)


def test_total_half_log_det_matches_half_logdet() -> None:
    """When all variables are eliminated, total_half_log_det = (1/2) logdet(D)."""
    D = _random_spd(4, seed=60)
    J = np.zeros(4)
    A = 0.0

    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("all", [0, 1, 2, 3])],
    )

    expected_log_det = 0.5 * np.linalg.slogdet(D)[1]
    assert tower.total_half_log_det == pytest.approx(expected_log_det)


# ── Order invariance ─────────────────────────────────────────────────


def test_elimination_order_invariance() -> None:
    """Schur complement composition is order-invariant: different orderings
    should produce the same final action."""
    D = _random_spd(4, seed=70)
    J = np.array([0.5, -0.3, 0.2, 0.1])
    A = 0.7

    # Order 1: [0,1] then [0,1]
    t1 = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("first", [0, 1]), ("second", [0, 1])],
    )

    # Order 2: [2,3] then [0,1] — but indices are relative to current,
    # so this eliminates original indices 2,3 first, then 0,1.
    t2 = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("last", [2, 3]), ("first", [0, 1])],
    )

    assert t1.final_action == pytest.approx(t2.final_action)
    assert t1.total_half_log_det == pytest.approx(t2.total_half_log_det)


# ── Input validation ─────────────────────────────────────────────────


def test_rejects_empty_stages() -> None:
    D = np.eye(2)
    J = np.zeros(2)
    with pytest.raises(InputValidationError, match="at least one"):
        tower_affine_hidden_elimination(0.0, J, D, stages=[])


def test_rejects_exhausted_hidden_space() -> None:
    """Trying to eliminate when no hidden coordinates remain."""
    D = np.eye(2)
    J = np.zeros(2)
    with pytest.raises(InputValidationError, match="no hidden"):
        tower_affine_hidden_elimination(
            0.0, J, D,
            stages=[("all", [0, 1]), ("ghost", [0])],
        )


def test_rejects_non_scalar_A() -> None:
    with pytest.raises(InputValidationError, match="scalar"):
        tower_affine_hidden_elimination(
            np.array([0.0, 1.0]), np.zeros(2), np.eye(2),
            stages=[("x", [0])],
        )


def test_rejects_bad_coupling_shape() -> None:
    with pytest.raises(InputValidationError):
        tower_affine_hidden_elimination(
            0.0, np.zeros((2, 2)), np.eye(2),
            stages=[("x", [0])],
        )


# ── Stage record metadata ───────────────────────────────────────────


def test_stage_records_have_correct_structure() -> None:
    """Check that TowerStageRecord fields are populated correctly."""
    D = _random_spd(4, seed=80)
    J = np.ones(4)
    A = 0.0

    tower = tower_affine_hidden_elimination(
        A, J, D,
        stages=[("block_a", [0, 1]), ("block_b", [0, 1])],
    )

    assert isinstance(tower, TowerEliminationRecord)
    for s in tower.stages:
        assert isinstance(s, TowerStageRecord)
        assert s.half_log_det_increment > 0  # SPD implies positive log-det
        assert s.schur_condition_number is not None
        assert s.schur_condition_number >= 1.0
