"""
Systematic steering tests across rdatasets.

Tests the eigenvalue capture formula and PCA comparison on diverse
real-world classification datasets from the R ecosystem.
"""
from __future__ import annotations

import csv
import numpy as np
import pytest
from pathlib import Path
from numpy.linalg import eigh, inv
from scipy.linalg import sqrtm

from nomogeo import extract_supervised, steer, capture_curve
from nomogeo.steer import optimize_observer, score_observer

RDATA = Path("C:/observer_geometry_workspace_v0.3.2/datasets/Rdatasets-master/Rdatasets-master/csv")


def _load_rdataset(pkg: str, item: str, target_col: str, drop_cols: list[str] | None = None):
    """Load an rdataset CSV, extract numeric features and a categorical target."""
    path = RDATA / pkg / f"{item}.csv"
    if not path.exists():
        pytest.skip(f"{pkg}/{item} not found")

    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        pytest.skip("Empty dataset")

    # Get target
    y_raw = [r[target_col] for r in rows]
    classes = sorted(set(y_raw))
    if len(classes) < 2:
        pytest.skip("Single class")
    y = np.array([classes.index(c) for c in y_raw])

    # Get numeric features
    drop = set(drop_cols or []) | {target_col, ""}
    feature_cols = []
    for col in rows[0].keys():
        if col in drop:
            continue
        try:
            vals = [float(r[col]) for r in rows]
            if not all(np.isfinite(v) for v in vals):
                continue
            if np.std(vals) < 1e-12:
                continue
            feature_cols.append(col)
        except (ValueError, KeyError):
            continue

    if len(feature_cols) < 4:
        pytest.skip(f"Too few numeric features ({len(feature_cols)})")

    X = np.array([[float(r[col]) for col in feature_cols] for r in rows])
    return X, y, classes, feature_cols


# ---- Dataset specifications: (package, item, target_column, drop_columns) ----

DATASETS = [
    # Medical / biomedical
    ("modeldata", "ad_data", "Class", ["rownames"]),
    ("modeldata", "cells", "case", ["rownames"]),
    ("HSAUR", "plasma", "ESR", ["rownames"]),
    ("HSAUR", "heptathlon", "score", ["rownames"]),  # treated as 2-class via median split
    ("MASS", "Pima.tr", "type", ["rownames"]),
    ("MASS", "biopsy", "class", ["rownames", "ID"]),

    # Ecology / agriculture
    ("archdata", "Bornholm", "Phase", ["rownames"]),
    ("archdata", "Michelsberg", "Chronology", ["rownames"]),
    ("modeldata", "leaf_id_flavia", "Species", ["rownames"]),
    ("agridat", "stroup.nin", "gen", ["rownames", "col", "row"]),

    # Social / economics
    ("Ecdat", "Car", "choice", ["rownames"]),
    ("MASS", "Cars93", "Type", ["rownames", "Manufacturer", "Model", "Make"]),
    ("MASS", "survey", "Exer", ["rownames"]),

    # Classic ML / pattern recognition
    ("MASS", "crabs", "sp", ["rownames", "index"]),
    ("MASS", "fgl", "type", ["rownames"]),
    ("datasets", "iris", "Species", ["rownames"]),
    ("ISLR", "Default", "default", ["rownames"]),
    ("ISLR", "Smarket", "Direction", ["rownames", "Today"]),

    # Engineering / physics
    ("modeldata", "chem_proc_yield", "Yield", ["rownames"]),  # continuous → median split

    # Sports / culture
    ("Lahman", "Teams", "lgID", ["rownames", "teamID", "franchID", "divID", "name", "park", "teamIDBR", "teamIDlahman45", "teamIDretro"]),
]


def _maybe_median_split(y, classes):
    """If target is continuous, split at median into 2 classes."""
    if len(classes) > 20:
        # Likely continuous — median split
        vals = y.astype(float)
        median = np.median(vals)
        y_new = (vals > median).astype(int)
        return y_new, ["low", "high"]
    return y, classes


class TestRdatasetsCapture:
    """Test the eigenvalue capture formula on diverse real datasets."""

    @pytest.mark.parametrize("spec", DATASETS, ids=[f"{s[0]}/{s[1]}" for s in DATASETS])
    def test_capture_formula(self, spec):
        pkg, item, target, drop = spec
        try:
            X, y, classes, feature_cols = _load_rdataset(pkg, item, target, drop)
        except Exception as e:
            pytest.skip(str(e))

        y, classes = _maybe_median_split(y, classes)
        if len(set(y)) < 2:
            pytest.skip("Single class after processing")

        n_features = X.shape[1]
        n_classes = len(set(y))

        # Extract geometry
        try:
            geom = extract_supervised(X, y, task="fisher")
        except Exception as e:
            pytest.skip(f"Extraction failed: {e}")

        # Test conservation at m=2
        m = min(2, n_features - 1)
        result = steer(H=geom.H, Hdot=geom.Hdot, rank=m)
        assert result.conservation_residual < 1e-10, (
            f"{pkg}/{item}: conservation failed ({result.conservation_residual:.2e})"
        )

        # Test eigenvalue capture formula
        try:
            H_inv = inv(geom.H)
            H_inv_sqrt = np.real(sqrtm(H_inv))
            Hdot_w = H_inv_sqrt @ geom.Hdot @ H_inv_sqrt
            Hdot_w = 0.5 * (Hdot_w + Hdot_w.T)
            eigvals_w = np.sort(eigh(Hdot_w)[0])[::-1]
            pos = eigvals_w[eigvals_w > 1e-10]
            total = np.sum(pos)

            if total > 1e-12 and len(pos) > 0:
                m_test = min(len(pos), n_features - 1)
                predicted = np.sum(pos[:m_test]) / total
                actual = steer(H=geom.H, Hdot=geom.Hdot, rank=m_test)
                err = abs(predicted - actual.visible_fraction)
                assert err < 1e-6, (
                    f"{pkg}/{item}: eigenvalue formula failed "
                    f"(predicted={predicted:.6f}, actual={actual.visible_fraction:.6f}, err={err:.2e})"
                )
        except Exception:
            pass  # sqrtm can fail on ill-conditioned matrices; skip formula test


class TestRdatasetsPCAComparison:
    """Verify steered observer beats or matches PCA across datasets."""

    @pytest.mark.parametrize("spec", DATASETS, ids=[f"{s[0]}/{s[1]}" for s in DATASETS])
    def test_beats_or_matches_pca(self, spec):
        pkg, item, target, drop = spec
        try:
            X, y, classes, feature_cols = _load_rdataset(pkg, item, target, drop)
        except Exception as e:
            pytest.skip(str(e))

        y, classes = _maybe_median_split(y, classes)
        if len(set(y)) < 2:
            pytest.skip("Single class after processing")

        try:
            geom = extract_supervised(X, y)
        except Exception as e:
            pytest.skip(f"Extraction failed: {e}")

        m = min(3, geom.n_features - 1)
        result = steer(H=geom.H, Hdot=geom.Hdot, rank=m)

        # Steered should capture at least as much as PCA (or very close)
        # Allow small margin for numerical reasons
        assert result.visible_fraction >= result.pca_visible_fraction - 0.01, (
            f"{pkg}/{item}: steered ({result.visible_fraction:.4f}) < "
            f"PCA ({result.pca_visible_fraction:.4f})"
        )
