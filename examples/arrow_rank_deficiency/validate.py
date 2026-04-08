from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from nomogeo import hidden_load, visible_precision

from examples.arrow_rank_deficiency.run_main import latent_precision, observer_ceiling, observers

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h = latent_precision()
    observer_map = observers()

    plurality = observer_map["plurality"]
    approval = observer_map["approval"]
    ranked = observer_map["ranked_choice"]

    phi_plurality = visible_precision(h, plurality)
    phi_approval = visible_precision(h, approval)
    phi_ranked = visible_precision(h, ranked)

    load_plurality = hidden_load(observer_ceiling("plurality", h), phi_plurality, support_mode="ambient")
    load_approval = hidden_load(observer_ceiling("approval", h), phi_approval, support_mode="ambient")
    load_ranked = hidden_load(observer_ceiling("ranked_choice", h), phi_ranked, support_mode="ambient")

    stage_ranked_to_approval = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    stage_approval_to_plurality = np.array([[1.0, 0.0]])
    composition_residual = float(
        np.linalg.norm(
            visible_precision(visible_precision(phi_ranked, stage_ranked_to_approval), stage_approval_to_plurality)
            - phi_plurality,
            ord=np.inf,
        )
    )

    report = {
        "composition_residual_abs": composition_residual,
        "clock_plurality": float(load_plurality.clock),
        "clock_approval": float(load_approval.clock),
        "clock_ranked_choice": float(load_ranked.clock),
        "clock_ordering_plurality_ge_approval_ge_ranked": bool(
            load_plurality.clock >= load_approval.clock >= load_ranked.clock
        ),
        "plurality_clock_minus_direct_gap_abs": abs(
            float(load_plurality.clock - (np.linalg.slogdet(h[:1, :1])[1] - np.linalg.slogdet(phi_plurality)[1]))
        ),
        "approval_clock_minus_direct_gap_abs": abs(
            float(load_approval.clock - (np.linalg.slogdet(h[:2, :2])[1] - np.linalg.slogdet(phi_approval)[1]))
        ),
        "ranked_clock_minus_direct_gap_abs": abs(
            float(load_ranked.clock - (np.linalg.slogdet(h)[1] - np.linalg.slogdet(phi_ranked)[1]))
        ),
        "rank_deficiencies": {
            "plurality": h.shape[0] - plurality.shape[0],
            "approval": h.shape[0] - approval.shape[0],
            "ranked_choice": h.shape[0] - ranked.shape[0],
        },
        "hidden_ranks": {
            "plurality": int(load_plurality.rank),
            "approval": int(load_approval.rank),
            "ranked_choice": int(load_ranked.rank),
        },
    }
    (OUT / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    thresholds = {
        "composition_residual_abs_max": 1e-10,
        "clock_minus_direct_gap_abs_max": 1e-10,
        "ranked_choice_clock_abs_max": 1e-10,
    }
    audit = {
        "example": "arrow_rank_deficiency",
        "claim": (
            "In this synthetic observer-geometry example, plurality, approval, and ranked choice erase different latent "
            "preference directions, and the hidden-load clocks quantify that loss."
        ),
        "headline_quantities": report,
        "thresholds": thresholds,
        "passes": {
            "composition": bool(report["composition_residual_abs"] <= thresholds["composition_residual_abs_max"]),
            "determinant_gap_cross_check": bool(
                max(
                    report["plurality_clock_minus_direct_gap_abs"],
                    report["approval_clock_minus_direct_gap_abs"],
                    report["ranked_clock_minus_direct_gap_abs"],
                )
                <= thresholds["clock_minus_direct_gap_abs_max"]
            ),
            "clock_ordering": bool(report["clock_ordering_plurality_ge_approval_ge_ranked"]),
            "ranked_choice_hidden_load_vanishes": bool(
                abs(report["clock_ranked_choice"]) <= thresholds["ranked_choice_clock_abs_max"]
            ),
        },
        "scripts": {
            "main": "python -m examples.arrow_rank_deficiency.run_main",
            "validate": "python -m examples.arrow_rank_deficiency.validate",
            "audit": "python -m examples.arrow_rank_deficiency.run_all",
        },
    }
    audit["passes"]["all"] = bool(all(audit["passes"].values()))
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
