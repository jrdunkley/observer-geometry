from __future__ import annotations

from pathlib import Path
import shutil

import pytest

from nomodescent import (
    AssumptionEntry,
    AssumptionLedger,
    GoalSpec,
    ObserverSpec,
    ProblemSpec,
    VisibleEvidenceSpec,
)
from nomodescent.exceptions import DescentInputError


def test_problem_spec_roundtrip_json() -> None:
    problem = ProblemSpec(
        name="roundtrip",
        latent_dim=2,
        observers=(
            ObserverSpec(name="full", matrix=[[1.0, 0.0], [0.0, 1.0]]),
            ObserverSpec(name="coarse", matrix=[[1.0, 0.0]]),
        ),
        evidence=(VisibleEvidenceSpec(name="coarse_cov", observer="coarse", kind="covariance", matrix=[[1.2]]),),
        assumptions=AssumptionLedger(entries=(AssumptionEntry(label="a1", statement="Gaussian regime", exact=True),)),
        goals=(GoalSpec(kind="common_completion"),),
    )
    tmp_dir = Path("tests") / "_tmp_roundtrip"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)
    path = tmp_dir / "problem.json"
    problem.to_json(path)
    loaded = ProblemSpec.from_json(path)
    assert loaded.to_dict() == problem.to_dict()
    shutil.rmtree(tmp_dir)


def test_problem_spec_rejects_unknown_observer_reference() -> None:
    with pytest.raises(DescentInputError):
        ProblemSpec(
            name="bad",
            latent_dim=2,
            observers=(ObserverSpec(name="full", matrix=[[1.0, 0.0], [0.0, 1.0]]),),
            evidence=(VisibleEvidenceSpec(name="e1", observer="missing", kind="covariance", matrix=[[1.0]]),),
        )


def test_problem_spec_rejects_nonsurjective_observer() -> None:
    with pytest.raises(DescentInputError):
        ObserverSpec(name="bad", matrix=[[1.0, 0.0], [2.0, 0.0]])


def test_visible_evidence_rejects_invalid_kind() -> None:
    with pytest.raises(DescentInputError):
        VisibleEvidenceSpec(name="e", observer="obs", kind="spectrum", matrix=[[1.0]])
