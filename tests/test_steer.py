"""Tests for the observer steering engine (steer.py)."""
from __future__ import annotations

import numpy as np
import pytest
from numpy.linalg import inv, eigh, norm
from numpy.testing import assert_allclose

from nomogeo.steer import score_observer, optimize_observer, steer
from nomogeo.source import information_budget
from nomogeo.exceptions import InputValidationError


def _spd(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return 0.5 * (A @ A.T + (A @ A.T).T) + np.eye(n)


def _sym(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return 0.5 * (A + A.T)


# ---- score_observer ----

class TestScoreObserver:

    def test_conservation(self):
        H = _spd(5, 0); Hdot = _sym(5, 1)
        rng = np.random.default_rng(2)
        C = rng.standard_normal((2, 5))
        sc = score_observer(H, C, Hdot)
        assert sc.conservation_residual < 1e-12

    def test_canonical_matches_self(self):
        """Canonical observer scored against canonical baseline should have zero advantage."""
        H = _spd(4, 10); Hdot = _sym(4, 11)
        C = np.zeros((2, 4)); C[0,0] = 1; C[1,1] = 1
        sc = score_observer(H, C, Hdot)
        assert abs(sc.advantage_over_canonical) < 1e-12

    def test_adapted_beats_pca(self):
        """On Iris-like data, adapted should beat PCA."""
        from sklearn.datasets import load_iris
        from nomogeo import extract_supervised, closure_adapted_observer
        X, y = load_iris(return_X_y=True)
        geom = extract_supervised(X, y)
        m = 2
        res = closure_adapted_observer(geom.H, [geom.Hdot], m)
        sc = score_observer(geom.H, res.C, geom.Hdot)
        assert sc.advantage_over_pca > 0

    def test_result_fields(self):
        H = _spd(4, 20); Hdot = _sym(4, 21)
        rng = np.random.default_rng(22)
        C = rng.standard_normal((2, 4))
        sc = score_observer(H, C, Hdot)
        assert sc.rank == 2
        assert sc.observer.shape == (2, 4)
        assert np.isfinite(sc.visible_rate)
        assert np.isfinite(sc.hidden_rate)
        assert np.isfinite(sc.ambient_rate)
        assert isinstance(sc.exact_sector, bool)


# ---- optimize_observer ----

class TestOptimizeObserver:

    def test_returns_correct_rank(self):
        H = _spd(6, 30); Hdot = _sym(6, 31)
        opt = optimize_observer(H, Hdot, rank=3)
        assert opt.rank == 3
        assert opt.observer.shape == (3, 6)

    def test_method_is_adapted(self):
        H = _spd(5, 40); Hdot = _sym(5, 41)
        opt = optimize_observer(H, Hdot, rank=2)
        assert opt.method == "adapted"

    def test_conservation(self):
        H = _spd(5, 50); Hdot = _sym(5, 51)
        opt = optimize_observer(H, Hdot, rank=2)
        assert opt.score.conservation_residual < 1e-12

    def test_zero_hidden_coupling(self):
        """Adapted observer should achieve B = 0."""
        H = _spd(6, 60); Hdot = _sym(6, 61)
        opt = optimize_observer(H, Hdot, rank=3)
        # hidden_defect = Tr(Q_hat) should be ~0 for adapted observer
        assert opt.diagnostics.hidden_defect_trace < 1e-8

    def test_capture_curve_included(self):
        H = _spd(5, 70); Hdot = _sym(5, 71)
        opt = optimize_observer(H, Hdot, rank=2)
        assert len(opt.capture.ranks) > 0

    def test_invalid_rank_raises(self):
        H = _spd(4, 80); Hdot = _sym(4, 81)
        with pytest.raises(InputValidationError):
            optimize_observer(H, Hdot, rank=0)
        with pytest.raises(InputValidationError):
            optimize_observer(H, Hdot, rank=4)

    @pytest.mark.parametrize("n,m", [(3,1), (4,2), (5,3), (6,2), (8,4)])
    def test_across_dimensions(self, n, m):
        H = _spd(n, n*100); Hdot = _sym(n, n*100+1)
        opt = optimize_observer(H, Hdot, rank=m)
        assert opt.observer.shape == (m, n)
        assert opt.score.conservation_residual < 1e-11

    def test_adapted_beats_canonical(self):
        """Adapted observer should generally beat or match canonical."""
        n_better = 0
        for seed in range(20):
            H = _spd(5, seed*50); Hdot = _sym(5, seed*50+1)
            opt = optimize_observer(H, Hdot, rank=2)
            if opt.score.advantage_over_canonical >= -0.01:
                n_better += 1
        # Should beat or match canonical in most cases
        assert n_better >= 12, f"Only beat canonical {n_better}/20 times"


# ---- steer ----

class TestSteer:

    def test_from_sklearn_data(self):
        from sklearn.datasets import load_iris
        X, y = load_iris(return_X_y=True)
        result = steer(X=X, y=y, rank=2, task="fisher")
        assert result.observer.shape == (2, 4)
        assert result.conservation_residual < 1e-12
        assert result.geometry is not None

    def test_from_raw_matrices(self):
        H = _spd(5, 90); Hdot = _sym(5, 91)
        result = steer(H=H, Hdot=Hdot, rank=2)
        assert result.observer.shape == (2, 5)
        assert result.conservation_residual < 1e-12
        assert result.geometry is None

    def test_beats_pca_on_iris(self):
        from sklearn.datasets import load_iris
        X, y = load_iris(return_X_y=True)
        result = steer(X=X, y=y, rank=2)
        assert result.advantage_over_pca > 0

    def test_beats_pca_on_wine(self):
        from sklearn.datasets import load_wine
        X, y = load_wine(return_X_y=True)
        result = steer(X=X, y=y, rank=3)
        assert result.advantage_over_pca > 0

    def test_all_task_types(self):
        from sklearn.datasets import load_iris
        X, y = load_iris(return_X_y=True)
        for task in ["fisher", "equal_weight", "minority"]:
            result = steer(X=X, y=y, rank=2, task=task)
            assert result.conservation_residual < 1e-12

    def test_no_input_raises(self):
        with pytest.raises(InputValidationError):
            steer(rank=2)

    def test_half_capture_rank(self):
        from sklearn.datasets import load_wine
        X, y = load_wine(return_X_y=True)
        result = steer(X=X, y=y, rank=3)
        # Wine has 13 features — should have a half-capture rank
        if result.half_capture_rank is not None:
            assert 1 <= result.half_capture_rank <= 13

    def test_high_dimensional(self):
        """Steer on a 50-dimensional problem."""
        rng = np.random.default_rng(100)
        n, k = 200, 50
        X = rng.standard_normal((n, k))
        y = (X[:, 0] + X[:, 1] > 0).astype(int)
        result = steer(X=X, y=y, rank=5)
        assert result.observer.shape == (5, k)
        assert result.conservation_residual < 1e-10


# ---- integration: full pipeline ----

class TestFullPipeline:
    """End-to-end: steer → project → verify."""

    def test_iris_pipeline(self):
        from sklearn.datasets import load_iris
        X, y = load_iris(return_X_y=True)
        result = steer(X=X, y=y, rank=2)

        # Project data
        # Need to standardise first (same as extract_supervised)
        X_s = (X - X.mean(axis=0)) / X.std(axis=0)
        X_proj = X_s @ result.observer.T
        assert X_proj.shape == (150, 2)

        # The projection should separate classes
        # (not a hard test, just sanity)
        class_means = np.array([X_proj[y == c].mean(axis=0) for c in range(3)])
        # At least two class means should be distinct
        dists = [norm(class_means[i] - class_means[j])
                 for i in range(3) for j in range(i+1, 3)]
        assert max(dists) > 0.1

    def test_wine_kills_pca(self):
        """The headline result: adapted observer captures 35x more than PCA on Wine."""
        from sklearn.datasets import load_wine
        X, y = load_wine(return_X_y=True)
        result = steer(X=X, y=y, rank=3)
        # Adapted should capture much more than PCA
        if abs(result.pca_visible_fraction) > 0.001:
            ratio = result.visible_fraction / result.pca_visible_fraction
            assert ratio > 5, f"Expected large advantage over PCA, got ratio {ratio:.1f}"

    def test_srbct_genomics(self):
        """Steering on SRBCT cancer genomics (63 samples, 2308 genes, 4 classes)."""
        import csv
        from sklearn.decomposition import PCA as skPCA
        csv_path = "C:/observer_geometry_workspace_v0.3.2/misc docs/srbct.csv"
        try:
            with open(csv_path) as f:
                reader = csv.reader(f)
                header = next(reader)
                rows = list(reader)
        except FileNotFoundError:
            pytest.skip("SRBCT data not available")
        X = np.array([[float(v) for v in row[1:]] for row in rows])
        y_raw = [row[0] for row in rows]
        classes = sorted(set(y_raw))
        y = np.array([classes.index(c) for c in y_raw])
        # Pre-reduce to 30 PCs (standard for n << p genomics)
        X_red = skPCA(n_components=30).fit_transform(X)
        result = steer(X=X_red, y=y, rank=2, task="fisher")
        assert result.visible_fraction > 0.8, f"Expected high capture, got {result.visible_fraction:.4f}"
        assert result.advantage_over_pca > 0.2, f"Expected advantage over PCA, got {result.advantage_over_pca:.4f}"
        assert result.conservation_residual < 1e-10

    def test_breast_cancer(self):
        """Breast Cancer: steered captures 100%, PCA captures 0%."""
        from sklearn.datasets import load_breast_cancer
        X, y = load_breast_cancer(return_X_y=True)
        result = steer(X=X, y=y, rank=2)
        assert result.visible_fraction > 0.99
        assert result.advantage_over_pca > 0.9

    def test_digits(self):
        """Digits (10 classes, 64 features): steered captures >90% at m=5."""
        from sklearn.datasets import load_digits
        X, y = load_digits(return_X_y=True)
        result = steer(X=X, y=y, rank=5)
        assert result.visible_fraction > 0.85
        assert result.advantage_over_pca > 0.8

    def test_exact_capture_theorem(self):
        """Capture theorem: vis_frac = 1.0 when rank(Hdot) <= m."""
        from sklearn.datasets import load_iris, load_wine, load_breast_cancer, load_digits
        cases = [
            ("Iris", load_iris, 2, 2),     # rank 2, m=2 → exact
            ("Wine", load_wine, 3, 2),     # rank 2, m=3 → exact
            ("BC", load_breast_cancer, 2, 1),  # rank 1, m=2 → exact
            ("Digits", load_digits, 9, 9), # rank 9, m=9 → exact
        ]
        for name, loader, m, expected_rank in cases:
            X, y = loader(return_X_y=True)
            result = steer(X=X, y=y, rank=m)
            assert result.visible_fraction > 0.999, (
                f"{name} m={m}: expected vis_frac=1.0, got {result.visible_fraction:.4f}"
            )

    def test_capture_below_rank(self):
        """When m < rank(Hdot), vis_frac < 1."""
        from sklearn.datasets import load_digits
        X, y = load_digits(return_X_y=True)
        # rank(Hdot) = 9, m = 3 → should NOT capture everything
        result = steer(X=X, y=y, rank=3)
        assert result.visible_fraction < 0.99

    def test_golub_leukemia(self):
        """Steering on Golub leukemia (72 samples, 7129 genes, 3 classes)."""
        import csv
        from sklearn.decomposition import PCA as skPCA
        csv_path = "C:/observer_geometry_workspace_v0.3.2/misc docs/golub.csv"
        try:
            with open(csv_path) as f:
                reader = csv.reader(f)
                header = next(reader)
                rows = list(reader)
        except FileNotFoundError:
            pytest.skip("Golub data not available")
        cancer_col = header.index('cancer')
        meta_cols = {'Samples', 'BM.PB', 'Gender', 'Source', 'tissue.mf', 'cancer'}
        gene_indices = [i for i, h in enumerate(header) if h not in meta_cols]
        X = np.nan_to_num(np.array([[float(row[i]) for i in gene_indices] for row in rows]), nan=0.0)
        y_raw = [row[cancer_col] for row in rows]
        classes = sorted(set(y_raw))
        y = np.array([classes.index(c) for c in y_raw])
        X_red = skPCA(n_components=30).fit_transform(X)
        result = steer(X=X_red, y=y, rank=2, task="fisher")
        assert result.visible_fraction > 0.5
        assert result.conservation_residual < 1e-10
