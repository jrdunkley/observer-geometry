from __future__ import annotations

import numpy as np
import pytest

from nomodescent import (
    ObserverSpec,
    classify_qd_relation,
    common_refinement_test,
    false_collapse_diagnostic,
)
from nomodescent.exceptions import DescentInputError


def test_common_refinement_rejects_mismatched_latent_dimensions() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 0.0]])
    observer_b = ObserverSpec(name="b", matrix=[[1.0, 0.0, 0.0]])
    with pytest.raises(DescentInputError, match="same latent dimension"):
        common_refinement_test((observer_a, observer_b))


def test_false_collapse_requires_multiple_summaries() -> None:
    observer = ObserverSpec(name="obs", matrix=[[1.0]])
    problem = type(
        "ProblemStub",
        (),
        {},
    )
    with pytest.raises(DescentInputError, match="at least two coarse summaries"):
        false_collapse_diagnostic(
            fine_problem=problem,  # type: ignore[arg-type]
            coarse_summaries={"only": np.array([1.0])},
            coarse_summary_label="single",
        )


def test_false_collapse_rejects_incompatible_summary_shapes() -> None:
    observer_full = ObserverSpec(name="full", matrix=np.eye(2))
    from nomodescent import GoalSpec, ProblemSpec, VisibleEvidenceSpec

    problem = ProblemSpec(
        name="toy",
        latent_dim=2,
        observers=(observer_full,),
        evidence=(VisibleEvidenceSpec(name="cov", observer="full", kind="covariance", matrix=[[1.0, 0.0], [0.0, 1.0]]),),
        goals=(GoalSpec(kind="common_completion"),),
    )
    with pytest.raises(DescentInputError, match="same shape"):
        false_collapse_diagnostic(
            fine_problem=problem,
            coarse_summaries={"a": np.array([1.0, 2.0]), "b": np.eye(2)},
            coarse_summary_label="bad",
        )


def test_nearly_factorising_pair_stays_non_nested_below_tolerance() -> None:
    observer_a = ObserverSpec(name="a", matrix=[[1.0, 1.0e-6]])
    observer_b = ObserverSpec(name="b", matrix=[[1.0, 0.0]])
    result = classify_qd_relation(observer_a, observer_b, tol=1e-12)
    assert result.classification == "non_nested_observers"
