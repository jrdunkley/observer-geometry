from __future__ import annotations

import json

import numpy as np

from nomogeo import declared_ladder_dimension_cost_intervals


def main() -> None:
    names = ["small_clean", "medium", "large"]
    scores = np.array([4.0, 5.8, 6.4])
    dimensions = np.array([1.0, 2.0, 4.0])
    result = declared_ladder_dimension_cost_intervals(scores, dimensions)
    rows = []
    for idx, name in enumerate(names):
        upper = result.interval_upper[idx]
        rows.append(
            {
                "candidate": name,
                "score": float(result.scores[idx]),
                "dimension": float(result.dimensions[idx]),
                "wins_for_cost_lower": float(result.interval_lower[idx]),
                "wins_for_cost_upper": None if np.isinf(upper) else float(upper),
                "interval_nonempty": bool(result.interval_nonempty[idx]),
            }
        )
    print(json.dumps({"declared_ladder_dimension_cost": rows}, indent=2))


if __name__ == "__main__":
    main()
