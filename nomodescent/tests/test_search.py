from __future__ import annotations

import numpy as np

from nomodescent import ObserverSpec, minimal_refinement_search


def test_minimal_refinement_search_is_deterministic() -> None:
    H = np.array(
        [
            [3.2, 0.6, 0.4, 0.0],
            [0.6, 2.7, 0.3, 0.1],
            [0.4, 0.3, 2.1, 0.0],
            [0.0, 0.1, 0.0, 1.8],
        ]
    )
    targets = [
        ObserverSpec(name="a", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
        ObserverSpec(name="b", matrix=[[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
    ]
    candidates = [
        ObserverSpec(name="panel", matrix=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
        ObserverSpec(name="full", matrix=np.eye(4)),
    ]
    first = minimal_refinement_search(H, targets, candidates, objective="min_dimension_then_logdet")
    second = minimal_refinement_search(H, targets, candidates, objective="min_dimension_then_logdet")
    assert first.winner == "panel"
    assert first.scores == second.scores
    assert first.audit.falsification_route
