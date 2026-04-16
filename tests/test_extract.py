"""Tests for the geometry extraction layer."""
from __future__ import annotations

import numpy as np
import pytest
from numpy.linalg import eigh, inv
from numpy.testing import assert_allclose

from nomogeo.extract import extract_supervised, extract_covariance, SupervisedGeometry
from nomogeo.source import information_budget
from nomogeo.exceptions import InputValidationError


def _iris():
    from sklearn.datasets import load_iris
    return load_iris(return_X_y=True)


def _wine():
    from sklearn.datasets import load_wine
    return load_wine(return_X_y=True)


class TestExtractSupervised:

    def test_returns_supervised_geometry(self):
        X, y = _iris()
        result = extract_supervised(X, y, task="fisher")
        assert isinstance(result, SupervisedGeometry)
        assert result.H.shape == (4, 4)
        assert result.Hdot.shape == (4, 4)
        assert result.n_features == 4
        assert result.n_classes == 3

    def test_H_is_spd(self):
        X, y = _iris()
        result = extract_supervised(X, y)
        eigvals = eigh(result.H)[0]
        assert np.min(eigvals) > 0

    def test_H_is_symmetric(self):
        X, y = _iris()
        result = extract_supervised(X, y)
        assert_allclose(result.H, result.H.T, atol=1e-14)

    def test_Hdot_is_symmetric(self):
        X, y = _iris()
        result = extract_supervised(X, y)
        assert_allclose(result.Hdot, result.Hdot.T, atol=1e-14)

    def test_Hdot_is_psd(self):
        """Between-class scatter should be PSD."""
        X, y = _iris()
        result = extract_supervised(X, y)
        eigvals = eigh(result.Hdot)[0]
        assert np.min(eigvals) > -1e-10

    def test_conservation_with_budget(self):
        """Extracted geometry should satisfy conservation law."""
        X, y = _iris()
        geom = extract_supervised(X, y)
        m = 2
        C = np.zeros((m, geom.n_features))
        C[0, 0] = 1.0; C[1, 1] = 1.0
        b = information_budget(geom.H, C, geom.Hdot)
        assert b.conservation_residual < 1e-12

    @pytest.mark.parametrize("task", ["fisher", "equal_weight", "minority"])
    def test_all_task_types(self, task):
        X, y = _iris()
        result = extract_supervised(X, y, task=task)
        assert result.H.shape == (4, 4)
        assert result.Hdot.shape == (4, 4)

    def test_class_scatters_count(self):
        X, y = _iris()
        result = extract_supervised(X, y)
        assert len(result.class_scatters) == 3

    def test_wine_dataset(self):
        X, y = _wine()
        result = extract_supervised(X, y)
        assert result.n_features == 13
        assert result.n_classes == 3
        eigvals = eigh(result.H)[0]
        assert np.min(eigvals) > 0

    def test_invalid_task_raises(self):
        X, y = _iris()
        with pytest.raises(InputValidationError):
            extract_supervised(X, y, task="invalid")

    def test_single_class_raises(self):
        X = np.random.randn(20, 3)
        y = np.zeros(20)
        with pytest.raises(InputValidationError):
            extract_supervised(X, y)

    def test_1d_X_raises(self):
        with pytest.raises(InputValidationError):
            extract_supervised(np.array([1, 2, 3]), np.array([0, 1, 0]))

    def test_adapted_observer_on_extracted_geometry(self):
        """Full pipeline: extract → adapted observer → budget."""
        from nomogeo import closure_adapted_observer
        X, y = _iris()
        geom = extract_supervised(X, y)
        m = 2
        result = closure_adapted_observer(geom.H, [geom.Hdot], m)
        b = information_budget(geom.H, result.C, geom.Hdot)
        assert b.conservation_residual < 1e-12
        # Adapted should capture more than canonical
        C_canon = np.zeros((m, geom.n_features))
        C_canon[0, 0] = 1.0; C_canon[1, 1] = 1.0
        b_canon = information_budget(geom.H, C_canon, geom.Hdot)
        # Not guaranteed to beat canonical in vis_frac (depends on geometry)
        # but conservation should hold for both
        assert b_canon.conservation_residual < 1e-12


class TestExtractCovariance:

    def test_basic(self):
        cov = np.array([[1.0, 0.5], [0.5, 1.0]])
        pert = np.array([[0.1, 0.0], [0.0, -0.1]])
        H, Hdot = extract_covariance(cov, pert)
        assert H.shape == (2, 2)
        assert Hdot.shape == (2, 2)
        # H should be SPD
        assert np.min(eigh(H)[0]) > 0
        # H and Hdot should be symmetric
        assert_allclose(H, H.T, atol=1e-14)
        assert_allclose(Hdot, Hdot.T, atol=1e-14)

    def test_conservation(self):
        rng = np.random.default_rng(42)
        n = 5
        A = rng.standard_normal((n, n))
        cov = 0.5 * (A @ A.T + (A @ A.T).T) + np.eye(n)
        pert = 0.5 * (rng.standard_normal((n, n)) + rng.standard_normal((n, n)).T) * 0.1
        H, Hdot = extract_covariance(cov, pert)
        m = 2
        C = np.zeros((m, n)); C[0, 0] = 1.0; C[1, 1] = 1.0
        b = information_budget(H, C, Hdot)
        assert b.conservation_residual < 1e-11

    def test_zero_perturbation(self):
        cov = np.eye(3)
        pert = np.zeros((3, 3))
        H, Hdot = extract_covariance(cov, pert)
        assert_allclose(Hdot, np.zeros((3, 3)), atol=1e-14)
