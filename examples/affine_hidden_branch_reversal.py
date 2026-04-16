from __future__ import annotations

import json

import numpy as np

from nomogeo import affine_hidden_branch_reversal, guarded_fibre_dominance


def main() -> None:
    variational = np.array([0.0, 0.2])
    fibre = np.array([0.5, -0.1])
    reversal = affine_hidden_branch_reversal(variational, fibre)
    flat_guard = guarded_fibre_dominance(np.array([2.0, 2.0]), np.array([-1.0, 1.0]), denominator_floor=1e-6)
    live_guard = guarded_fibre_dominance(np.array([0.0, 2.0]), np.array([-1.0, 1.0]), denominator_floor=1e-6)
    print(
        json.dumps(
            {
                "affine_hidden_branch_reversal": {
                    "variational_winners": list(reversal.variational_winners),
                    "visible_winners": list(reversal.visible_winners),
                    "reversal": reversal.reversal,
                    "visible_action": reversal.visible_action.tolist(),
                },
                "guarded_fibre_dominance": {
                    "flat_variational_ratio_defined": flat_guard.ratio_defined,
                    "live_variational_ratio_defined": live_guard.ratio_defined,
                    "live_ratio": live_guard.ratio,
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
